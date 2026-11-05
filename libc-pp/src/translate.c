#include <errno.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "allocator.h"
#include "ast.h"
#include "diag.h"
#include "parser.tab.h"
#include "scope.h"
#include "translate.h"

/* printf into the arena, so that the string lasts as long as the tree. */
static char *arena_printf (const char *format, ...)
{
  va_list ap;
  char *s;
  int n;

  va_start (ap, format);
  n = vsnprintf (NULL, 0, format, ap);
  va_end (ap);
  s = allocate ((size_t)n + 1);
  va_start (ap, format);
  vsnprintf (s, (size_t)n + 1, format, ap);
  va_end (ap);
  return s;
}

/* see if a token of kind 'kind' is a child of the subtree 'specs'. */
static int has_token (const struct ast *specs, int kind, const char *text)
{
  for (const struct ast *n = specs->child; n; n = n->next)
    if (n->kind == AST_TOKEN && ast_token (n)->kind == kind
        && (!text || strcmp (ast_token (n)->text, text) == 0))
      return 1;
  return 0;
}

/* A plain name, as in `scalar a;` or `scalar a __attribute__ ((unused));`:
 * not a pointer, array or function, and with no initializer yet. */
static int plain_name (const struct ast *d)
{
  if (d->kind == AST_ATTRIBUTED_DECLARATOR)
    d = d->child;
  return d->kind == AST_DECL_NAME;
}

static int is_scalar_specs (const struct ast *specs)
{
  return has_token (specs, T_TYPE_NAME, "scalar");
}

/* The specifiers of the declaration or parameter where d, a declarator, is a
 * plain name, or NULL. */
static const struct ast *declaring_specs (const struct ast *d)
{
  if (d->kind != AST_DECL_NAME)
    return NULL;
  for (const struct ast *n = d->parent; n; n = n->parent) {
    if (n->kind == AST_DECLARATION || n->kind == AST_PARAMETER)
      return n->child;
    if (n->kind != AST_ATTRIBUTED_DECLARATOR && n->kind != AST_INIT_DECLARATOR
        && n->kind != AST_INIT_DECLARATOR_LIST)
      return NULL;
  }
  return NULL;
}

void check_scalar_assignment (const struct scope_table *scopes,
                              const struct ast *target, const struct ast *op)
{
  const struct token *name, *at = ast_token (op);
  const struct var *v;
  const struct ast *specs;

  while (target->kind == AST_PAREN_EXPRESSION)
    target = target->child->next;
  if (target->kind != AST_TOKEN)
    return;
  name = ast_token (target);
  if (name->kind != T_IDENTIFIER && name->kind != T_FOREACH)
    return;
  v = scope_find (scopes, name->text);
  if (!v || v->kind != SCOPE_VARIABLE || !v->decl)
    return;
  specs = declaring_specs (v->decl);
  if (!specs || !is_scalar_specs (specs))
    return;
  diag (at->file, at->line, DIAG_WARNING, "'%s' changes scalar '%s', whose "
        "value cwrap gives it", at->text, name->text);
  at = ast_first_token (v->decl);
  diag (at->file, at->line, DIAG_NOTE, "'%s' is declared here", name->text);
}

static void warn_scalar_initializer (const struct ast *init)
{
  const struct token *t = ast_first_token (init);
  const char *name = ast_declared_name (init->child);

  diag (t->file, t->line, DIAG_WARNING, "scalar '%s' has a value of its own, "
        "so cwrap does not number it", name);
}

/* Parses text, external declarations that may use the type scalar, for code
 * that the translator adds. Returns the first declaration, chained to the rest
 * by next. */
static struct ast *parse_snippet (const char *text)
{
  static const char prelude[] = "typedef int scalar;";
  size_t n = strlen (prelude) + strlen (text) + 1;
  char *source = allocate (n);
  struct ast *root, *first, *d;
  FILE *in;
  int rejected;

  snprintf (source, n, "%s%s", prelude, text);
  in = fmemopen (source, n - 1, "r");
  if (!in) {
    diag_tool (DIAG_ERROR, "fmemopen: %s", strerror (errno));
    exit (1);
  }
  root = parse_file (in, &rejected);
  fclose (in);
  if (!root || root->child->next->kind == AST_TOKEN)
    abort ();
  first = root->child->next;              /* past the prelude */
  for (d = first; d->next->kind != AST_TOKEN; d = d->next)
    d->parent = NULL;
  d->parent = d->next = NULL;             /* without the end of input */
  return first;
}

/* Puts the chain of declarations from parse_snippet () before pos. */
static void insert_all_before (struct ast *pos, struct ast *decls)
{
  while (decls) {
    struct ast *next = decls->next;

    ast_insert_before (pos, decls);
    decls = next;
  }
}

/* The first token under n, to change its trivia. */
static struct token *first_token (struct ast *n)
{
  while (n->kind != AST_TOKEN)
    n = n->child;
  return &((struct ast_leaf *)n)->tok;
}

/* The name of the index of a foreach loop over the variable %s. */
#define FOREACH_INDEX "__i_%s"

/* Gives every token under n the place of `at`, for the AST dump. */
static void set_location (struct ast *n, const struct token *at)
{
  if (n->kind == AST_TOKEN) {
    struct token *t = &((struct ast_leaf *)n)->tok;

    t->file = at->file;
    t->line = at->line;
    return;
  }
  for (n = n->child; n; n = n->next)
    set_location (n, at);
}

/* The trivia of the tokens from first to last, in order. */
static const char *joined_trivia (struct ast *const *tokens, int n)
{
  size_t len = 1;
  char *s;

  for (int i = 0; i < n; i++)
    len += strlen (ast_token (tokens[i])->trivia);
  s = allocate (len);
  for (int i = 0; i < n; i++)
    strcat (s, ast_token (tokens[i])->trivia);
  return s;
}

struct ast *foreach_to_for (struct ast *kw, struct ast *lparen, struct ast *var,
                            struct ast *in, struct ast *coll, struct ast *rparen,
                            struct ast *body, int is_cache)
{
  const char *v = ast_token (var)->text, *c = ast_token (coll)->text;
  struct ast *header[] = { var, in, coll, rparen };
  const char *index = arena_printf (FOREACH_INDEX, v);
  const char *start = "0";
  struct ast *fn, *loop, *block, *n;
  const char *trivia;

  /* A dirty cache is updated once, before the first element is read. */
  if (is_cache)
    start = arena_printf ("%s.update && %s.is_dirty ? (%s.is_dirty = 0, "
                          "%s.update (&%s), 0) : 0", c, c, c, c, c);
  fn = parse_snippet (arena_printf ("void f (void) { for (int %s = %s; %s < %s.n; "
                                    "++%s) { __auto_type %s = %s.a[%s]; } }",
                                    index, start, index, c, index, v, c, index));
  /* void f (void) { for ... } : the for statement is in the body's items. */
  loop = fn->child->next->next->child->next->child;
  loop->parent = loop->next = NULL;
  set_location (loop, ast_token (kw));

  /* The header keeps its lines: `for` and `(` take the trivia of `foreach` and
   * `(`, and the `)` takes all of the trivia of `v in V)` if it spans lines. */
  first_token (loop)->trivia = ast_token (kw)->trivia;
  first_token (loop->child->next)->trivia = ast_token (lparen)->trivia;
  trivia = joined_trivia (header, 4);
  for (n = loop->child; n; n = n->next)
    if (n->kind == AST_TOKEN && strcmp (ast_token (n)->text, ")") == 0)
      first_token (n)->trivia = strchr (trivia, '\n') ? trivia
                                                       : ast_token (rparen)->trivia;

  /* { __auto_type v = V.a[i]; } gets the body after its declaration. */
  block = loop->last;
  if (body)
    ast_add (block->child->next, body);
  return loop;
}

/* Frees a local scalar when its scope ends, as scalar.h defines it: cleanup
 * passes the variable's address. */
#define SCALAR_CLEANUP "__cwrap_scalar_free"

/* For the names of the slots of static scalars in functions. */
static int local_slot;

/* The first item of the block in `void f (void) { text }`. */
static struct ast *parse_block_item (const char *text)
{
  struct ast *fn = parse_snippet (arena_printf ("void f (void) { %s }", text));
  struct ast *item = fn->child->next->next->child->next->child;

  item->parent = item->next = NULL;
  return item;
}

/* Puts n just after d, among the children of d's parent. */
static void insert_after (struct ast *d, struct ast *n)
{
  if (d->next) {
    ast_insert_before (d->next, n);
    return;
  }
  n->next = NULL;
  ast_add (d->parent, n);
}

void local_scalars (struct ast *decl)
{
  struct ast *specs, *list;
  int is_static;

  if (decl->kind == AST_EXTENSION_DECLARATION)
    decl = decl->child->next;
  if (decl->kind != AST_DECLARATION)
    return;
  specs = decl->child;
  list = specs->next;
  if (list->kind != AST_INIT_DECLARATOR_LIST || !is_scalar_specs (specs)
      || has_token (specs, T_TYPEDEF, NULL) || has_token (specs, T_EXTERN, NULL)
      || has_token (specs, T_REGISTER, NULL) || has_token (specs, T_THREAD_LOCAL, NULL))
    return;
  is_static = has_token (specs, T_STATIC, NULL);
  for (struct ast *d = list->child, *next; d; d = next) {
    const char *name = ast_declared_name (d);

    next = d->next;
    if (!plain_name (d))
      continue;
    if (is_static) {
      /* static scalar s, *const slot ... = &s; */
      struct ast *more = parse_snippet (arena_printf (
        "static scalar __v, *const __cwrap_local_slot_%d "
        "__attribute__ ((section (\"scalars\"), used)) = &%s;", local_slot++, name));
      struct ast *comma = more->child->next->child->next, *slot = comma->next;

      insert_after (d, comma);
      insert_after (comma, slot);
    } else {
      /* scalar x __attribute__ ((cleanup (...))) = scalar_new (); */
      struct ast *item = parse_block_item (arena_printf (
        "scalar __v __attribute__ ((cleanup (%s))) = scalar_new ();",
        SCALAR_CLEANUP));
      struct ast *init = item->child->next->child;
      struct ast *placeholder = init->child->child;

      ast_replace (d, init);
      ast_replace (placeholder, d);
    }
  }
}

/* Each `scalar a;` at file scope, headers included, gets a slot in the section
 * scalars, on the same line so that no line moves:
 *
 *   scalar a; static scalar *const __cwrap_scalar_slot_0 ... = &a;
 *
 * The linker joins the slots of every file into one array, which scalar.c
 * numbers at startup. extern and typedef
 * declarations, pointers, arrays and functions get no slot, and neither do
 * names with an initializer of their own, which get a warning. */
static void register_scalars (struct ast *root)
{
  int slot = 0;

  for (struct ast *n = root->child, *next; n; n = next) {
    struct ast *decl = n, *specs, *list;

    next = n->next;
    if (decl->kind == AST_EXTENSION_DECLARATION)
      decl = decl->child->next;
    if (decl->kind != AST_DECLARATION)
      continue;
    specs = decl->child;
    list = specs->next;
    if (list->kind != AST_INIT_DECLARATOR_LIST
        || !is_scalar_specs (specs)
        || has_token (specs, T_TYPEDEF, NULL) || has_token (specs, T_EXTERN, NULL))
      continue;
    for (struct ast *d = list->child; d; d = d->next) {
      struct ast *reg;

      if (d->kind == AST_INIT_DECLARATOR && plain_name (d->child))
        warn_scalar_initializer (d);
      if (!plain_name (d))
        continue;
      reg = parse_snippet (arena_printf ("static scalar *const __cwrap_scalar_slot_%d "
                                         "__attribute__ ((section (\"scalars\"), used)) = &%s;",
                                         slot++, ast_declared_name (d)));
      first_token (reg)->trivia = " ";
      insert_all_before (next, reg);
    }
  }
}


/* Caches */

/* Frees the array of a local cache when its scope ends: see <cwrap/cache.h>. */
#define CACHE_CLEANUP "__cwrap_cache_free"

/* The element types of cache (T) seen in this file, by struct tag. */
struct cache_type {
  const char *tag;        /* cache_int */
  const char *type;       /* T as written, tokens one space apart */
  const char *typedefs;   /* `typedef int N;` for each typedef name in T, to
                           * parse the definition apart from the file */
  int defined;
  struct cache_type *next;
};
static struct cache_type *cache_types;

static struct cache_type *find_cache_type (const char *tag)
{
  for (struct cache_type *ct = cache_types; ct; ct = ct->next)
    if (strcmp (ct->tag, tag) == 0)
      return ct;
  return NULL;
}

/* Writes what cache_type needs of the tokens under n: T's text, the tag, in
 * which `*` is p and any other punctuator its code, and the typedef names. */
static void describe_type (const struct ast *n, FILE *text, FILE *tag,
                           FILE *typedefs)
{
  if (n->kind != AST_TOKEN) {
    for (n = n->child; n; n = n->next)
      describe_type (n, text, tag, typedefs);
    return;
  }
  const struct token *t = ast_token (n);
  unsigned char c = (unsigned char)t->text[0];

  fprintf (text, "%s%s", ftell (text) ? " " : "", t->text);
  if (c == '_' || c == '$' || c >= 0x80 || (c >= 'A' && c <= 'Z')
      || (c >= 'a' && c <= 'z') || (c >= '0' && c <= '9'))
    fprintf (tag, "_%s", t->text);
  else if (strcmp (t->text, "*") == 0)
    fprintf (tag, "_p");
  else
    for (const char *p = t->text; *p; p++)
      fprintf (tag, "_x%02x", (unsigned char)*p);
  if (t->kind == T_TYPE_NAME)
    fprintf (typedefs, "typedef int %s;", t->text);
}

/* The trivia of every token under n, in order. */
static void subtree_trivia (const struct ast *n, FILE *out)
{
  if (n->kind == AST_TOKEN) {
    fputs (ast_token (n)->trivia, out);
    return;
  }
  for (n = n->child; n; n = n->next)
    subtree_trivia (n, out);
}

/* What a FILE from open_memstream () holds, now in the arena. */
static char *memstream_text (FILE *f, char **buf, size_t *len)
{
  char *s;

  if (fclose (f) != 0)
    abort ();
  s = allocate (*len + 1);
  memcpy (s, *buf, *len);
  s[*len] = '\0';
  free (*buf);
  return s;
}

#define MEMSTREAM(name)                                                 \
  char *name##_buf = NULL;                                              \
  size_t name##_len = 0;                                                \
  FILE *name = open_memstream (&name##_buf, &name##_len)
#define MEMSTREAM_TEXT(name) memstream_text (name, &name##_buf, &name##_len)

struct ast *cache_specifier (struct ast *kw, struct ast *lparen, struct ast *type,
                             struct ast *rparen, int *rejected)
{
  MEMSTREAM (text);
  MEMSTREAM (tag);
  MEMSTREAM (typedefs);
  MEMSTREAM (trivia);
  struct ast *parts[] = { lparen, type, rparen };
  const char *type_text, *tag_text, *between;
  struct cache_type *ct;
  struct ast *spec;

  if (!text || !tag || !typedefs || !trivia)
    abort ();
  fputs ("cache", tag);
  describe_type (type, text, tag, typedefs);
  for (int i = 0; i < 3; i++)
    subtree_trivia (parts[i], trivia);
  type_text = MEMSTREAM_TEXT (text);
  tag_text = MEMSTREAM_TEXT (tag);
  between = MEMSTREAM_TEXT (trivia);

  ct = find_cache_type (tag_text);
  if (!ct) {
    ct = allocate (sizeof (*ct));
    ct->tag = tag_text;
    ct->type = type_text;
    ct->typedefs = MEMSTREAM_TEXT (typedefs);
    ct->defined = 0;
    ct->next = cache_types;
    cache_types = ct;
  } else {
    MEMSTREAM_TEXT (typedefs);
    if (strcmp (ct->type, type_text) != 0) {
      const struct token *t = ast_token (kw);

      diag (t->file, t->line, DIAG_ERROR, "cache (%s) and cache (%s) would "
            "both be struct %s", type_text, ct->type, tag_text);
      (*rejected)++;
    }
  }

  /* struct cache_T *__p; : its specifier. `struct` takes the trivia of `cache`,
   * and the tag that of the rest, if it spans lines. */
  spec = parse_snippet (arena_printf ("struct %s *__p;", tag_text))->child->child;
  spec->parent = spec->next = NULL;
  set_location (spec, ast_token (kw));
  first_token (spec)->trivia = ast_token (kw)->trivia;
  first_token (spec->last)->trivia = strchr (between, '\n') ? between : " ";
  return ast_node (AST_CACHE_SPECIFIER, &spec, 1);
}

static int is_cache_specs (const struct ast *specs)
{
  for (const struct ast *n = specs->child; n; n = n->next)
    if (n->kind == AST_CACHE_SPECIFIER)
      return 1;
  return 0;
}

const struct ast *cache_variable (const struct scope_table *scopes,
                                  const struct ast *expr)
{
  const struct token *name;
  const struct var *v;
  const struct ast *specs;

  while (expr->kind == AST_PAREN_EXPRESSION)
    expr = expr->child->next;
  if (expr->kind != AST_TOKEN)
    return NULL;
  name = ast_token (expr);
  if (name->kind != T_IDENTIFIER && name->kind != T_FOREACH && name->kind != T_CACHE)
    return NULL;
  v = scope_find (scopes, name->text);
  if (!v || v->kind != SCOPE_VARIABLE || !v->decl)
    return NULL;
  specs = declaring_specs (v->decl);
  return specs && is_cache_specs (specs) ? v->decl : NULL;
}

static void note_declared (const struct ast *decl)
{
  const struct token *t = ast_first_token (decl);

  diag (t->file, t->line, DIAG_NOTE, "'%s' is declared here",
        ast_declared_name (decl));
}

int check_cache_copy (const struct scope_table *scopes, const struct ast *expr,
                      const struct ast *at)
{
  const struct ast *decl = cache_variable (scopes, expr);
  const struct token *t = ast_token (at);

  if (!decl)
    return 0;
  diag (t->file, t->line, DIAG_ERROR, "'%s' copies cache '%s', and two caches "
        "must not share an array", t->text, ast_declared_name (decl));
  note_declared (decl);
  return 1;
}

void check_cache_arguments (const struct scope_table *scopes,
                            const struct ast *args)
{
  for (const struct ast *n = args->child; n; n = n->next) {
    const struct ast *decl = cache_variable (scopes, n);
    const struct token *t;

    if (!decl)
      continue;
    t = ast_first_token (n);
    diag (t->file, t->line, DIAG_WARNING, "cache '%s' is passed by value, and "
          "the copy shares its array; pass &%s", ast_declared_name (decl),
          ast_declared_name (decl));
  }
}

int check_cache_declaration (const struct ast *decl)
{
  const struct token *t;

  if (!is_cache_specs (decl->child) || !has_token (decl->child, T_TYPEDEF, NULL))
    return 0;
  t = ast_first_token (decl);
  diag (t->file, t->line, DIAG_ERROR, "a typedef of a cache type: write cache "
        "(T) where it is used");
  return 1;
}

int local_caches (struct ast *decl)
{
  struct ast *specs, *list;
  int errors = 0;

  if (decl->kind == AST_EXTENSION_DECLARATION)
    decl = decl->child->next;
  if (decl->kind != AST_DECLARATION)
    return 0;
  specs = decl->child;
  list = specs->next;
  /* Static storage starts at zero and lasts, and register has no address for
   * the cleanup to take. */
  if (list->kind != AST_INIT_DECLARATOR_LIST || !is_cache_specs (specs)
      || has_token (specs, T_TYPEDEF, NULL) || has_token (specs, T_EXTERN, NULL)
      || has_token (specs, T_STATIC, NULL) || has_token (specs, T_THREAD_LOCAL, NULL)
      || has_token (specs, T_REGISTER, NULL))
    return 0;
  for (struct ast *d = list->child, *next; d; d = next) {
    struct ast *target = d->kind == AST_INIT_DECLARATOR ? d->child : d;

    next = d->next;
    if (!plain_name (target))
      continue;
    if (d->kind == AST_INIT_DECLARATOR) {
      /* x = { ... } becomes x __attribute__ ((cleanup (...))) = { ... } */
      struct ast *item, *attributed;

      if (d->last->kind != AST_INITIALIZER_LIST_EXPR) {
        const struct token *t = ast_first_token (d->last);

        diag (t->file, t->line, DIAG_ERROR, "cache '%s' is initialized with a "
              "value, and two caches must not share an array; use { ... } or "
              "cache_push ()", ast_declared_name (target));
        errors++;
        continue;
      }
      item = parse_block_item (arena_printf (
        "int __v __attribute__ ((cleanup (%s)));", CACHE_CLEANUP));
      attributed = item->child->next->child;
      ast_replace (target, attributed);
      ast_replace (attributed->child, target);
    } else {
      /* x becomes x __attribute__ ((cleanup (...))) = { 0 } */
      struct ast *item = parse_block_item (arena_printf (
        "int __v __attribute__ ((cleanup (%s))) = { 0 };", CACHE_CLEANUP));
      struct ast *init = item->child->next->child;

      ast_replace (d, init);
      ast_replace (init->child->child, d);
    }
  }
  return errors;
}

/* Puts the definition of each struct cache_T used under n, inner ones first,
 * before ext, the external declaration n is in, if this file has none yet.
 * The definition goes on ext's line, after its trivia, so no line moves. */
static void define_caches_in (struct ast *n, struct ast *ext)
{
  struct cache_type *ct;
  struct ast *def;

  if (n->kind == AST_TOKEN)
    return;
  for (struct ast *k = n->child; k; k = k->next)
    define_caches_in (k, ext);
  if (n->kind != AST_CACHE_SPECIFIER)
    return;
  ct = find_cache_type (ast_token (n->child->last)->text);
  if (!ct || ct->defined)
    return;
  ct->defined = 1;
  def = parse_snippet (arena_printf ("%sstruct %s { __typeof__ (%s) *a; "
                                     "int n, cap, is_dirty; "
                                     "void (*update) (struct %s *self); };",
                                     ct->typedefs, ct->tag, ct->type, ct->tag));
  while (def->next)                       /* past the typedefs */
    def = def->next;
  set_location (def, ast_first_token (ext));
  first_token (def)->trivia = first_token (ext)->trivia;
  first_token (ext)->trivia = " ";
  ast_insert_before (ext, def);
}

static void define_caches (struct ast *root)
{
  for (struct ast *n = root->child; n; n = n->next)
    define_caches_in (n, n);
}

/* The translation itself: rewrite the tree here. */
static void transform (struct ast *root)
{
  register_scalars (root);
  define_caches (root);
}

/* gcc gets the input as it came. It either accepts code our grammar does not
 * cover, or reports the error better than we could. */
static int copy_file (FILE *in, FILE *out)
{
  char buf[1 << 16];
  size_t n;

  rewind (in);
  while ((n = fread (buf, 1, sizeof (buf), in)) > 0)
    if (fwrite (buf, 1, n, out) != n)
      return 1;
  return ferror (in);
}

int translate_file (const char *in_path, const char *out_path)
{
  FILE *in, *out;
  struct ast *root;
  int failed, parsed, rejected;

  in = fopen (in_path, "r");
  if (!in) {
    diag_tool (DIAG_ERROR, "%s: %s", in_path, strerror (errno));
    return TRANSLATE_FAILED;
  }
  cache_types = NULL;
  root = parse_file (in, &rejected);
  if (rejected) {
    fclose (in);
    free_arena ();
    return TRANSLATE_FAILED;
  }
  parsed = root != NULL;
  if (root) {
    transform (root);
    if (getenv ("CWRAP_DUMP_AST"))
      ast_dump (root, stderr);
  }

  out = fopen (out_path, "w");
  if (!out) {
    diag_tool (DIAG_ERROR, "%s: %s", out_path, strerror (errno));
    fclose (in);
    free_arena ();
    return TRANSLATE_FAILED;
  }
  if (root)
    ast_emit (root, out);
  failed = root ? 0 : copy_file (in, out);
  fclose (in);
  free_arena ();
  failed |= ferror (out);
  if (fclose (out) != 0 || failed) {
    diag_tool (DIAG_ERROR, "%s: write failed", out_path);
    return TRANSLATE_FAILED;
  }
  return parsed ? TRANSLATE_OK : TRANSLATE_SKIPPED;
}
