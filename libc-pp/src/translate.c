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

/* The euclidean dimension from -dimension=, or 0: then scalar, vector and
 * tensor are ordinary names, and nothing here changes them. */
static int dimension;

void translate_set_dimension (int d)
{
  dimension = d;
}

/* The scalars in a tensor, at most */
#define MAX_PARTS 9

/* The type of cwrap's that specs name: scalar, vector or tensor, or NULL. */
static const char *scalar_type (const struct ast *specs)
{
  static const char *const types[] = { "scalar", "vector", "tensor" };

  for (int i = 0; i < (dimension >= 2 ? 3 : dimension); i++)
    if (has_token (specs, T_TYPE_NAME, types[i]))
      return types[i];
  return NULL;
}

/* The scalars in a variable of the type that specs name, as the members to
 * append to its name: "" for a scalar, ".x" ".y" for a vector, ".x.x" ".x.y"
 * ".y.x" ".y.y" for a tensor in two dimensions. Returns how many: 0 if the
 * type is not one of cwrap's. */
static int scalar_parts (const struct ast *specs, const char *parts[MAX_PARTS])
{
  static const char *const axes[] = { "x", "y", "z" };
  const char *type = scalar_type (specs);
  int n = 0;

  if (!type)
    return 0;
  if (strcmp (type, "scalar") == 0)
    parts[n++] = "";
  else if (strcmp (type, "vector") == 0)
    for (int i = 0; i < dimension; i++)
      parts[n++] = arena_printf (".%s", axes[i]);
  else
    for (int i = 0; i < dimension; i++)
      for (int j = 0; j < dimension; j++)
        parts[n++] = arena_printf (".%s.%s", axes[i], axes[j]);
  return n;
}

static int is_scalar_specs (const struct ast *specs)
{
  return scalar_type (specs) != NULL;
}

/* The tokens under n, with no space between them, as in `v.x`. */
static const char *expression_text (const struct ast *n)
{
  if (n->kind == AST_TOKEN)
    return ast_token (n)->text;
  const char *s = "";
  for (n = n->child; n; n = n->next)
    s = arena_printf ("%s%s", s, expression_text (n));
  return s;
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

/* target may also be a member of a vector or tensor, as in t.x or t.x.y. */
void check_scalar_assignment (const struct scope_table *scopes,
                              const struct ast *target, const struct ast *op)
{
  static const char *const by_depth[] = { "scalar", "vector", "tensor" };
  const struct token *name, *at = ast_token (op);
  const struct ast *root = target, *specs;
  const struct var *v;
  const char *type;
  int members = 0, depth;

  /* (t.x).y: down to t, counting the members */
  for (;;) {
    while (root->kind == AST_PAREN_EXPRESSION)
      root = root->child->next;
    if (root->kind != AST_MEMBER || strcmp (ast_token (root->child->next)->text, ".") != 0)
      break;
    root = root->child;
    members++;
  }
  if (root->kind != AST_TOKEN)
    return;
  name = ast_token (root);
  if (name->kind != T_IDENTIFIER && name->kind != T_FOREACH && name->kind != T_CACHE)
    return;
  v = scope_find (scopes, name->text);
  if (!v || v->kind != SCOPE_VARIABLE || !v->decl)
    return;
  specs = declaring_specs (v->decl);
  if (!specs || !(type = scalar_type (specs)))
    return;
  depth = strcmp (type, "scalar") == 0 ? 0 : strcmp (type, "vector") == 0 ? 1 : 2;
  if (members > depth)
    return;
  diag (at->file, at->line, DIAG_WARNING, "'%s' changes %s '%s', whose value "
        "cwrap gives it", at->text, by_depth[depth - members],
        expression_text (target));
  at = ast_first_token (v->decl);
  diag (at->file, at->line, DIAG_NOTE, "'%s' is declared here", name->text);
}

static void warn_scalar_initializer (const struct ast *init, const char *type)
{
  const struct token *t = ast_first_token (init);
  const char *name = ast_declared_name (init->child);

  diag (t->file, t->line, DIAG_WARNING, "%s '%s' has a value of its own, so "
        "cwrap does not number it", type, name);
}

/* Parses text, external declarations that may use the type scalar, for code
 * that the translator adds. Returns the first declaration, chained to the rest
 * by next. */
static struct ast *parse_snippet (const char *text)
{
  static const char prelude[] = "typedef int scalar; typedef int vector; "
                                 "typedef int tensor;";
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
  if (!root)
    abort ();
  first = root->child->next->next->next;  /* past the three typedefs */
  if (first->kind == AST_TOKEN)
    abort ();
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

/* Code that cwrap writes as text has no declarations behind its names: each
 * name spelled `name` under n gets the declaration of the token it was
 * written from, so that a @loop body renames it along with the rest. */
static void bind_spelled (struct ast *n, const struct token *from)
{
  if (n->kind == AST_TOKEN) {
    struct token *t = &((struct ast_leaf *)n)->tok;

    if (!t->binding && strcmp (t->text, from->text) == 0)
      t->binding = from->binding;
    return;
  }
  for (n = n->child; n; n = n->next)
    bind_spelled (n, from);
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

  bind_spelled (loop, ast_token (var));
  bind_spelled (loop, ast_token (coll));

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

struct ast *local_scalars (struct ast *decl)
{
  struct ast *item = decl, *specs, *list, *with_slots;
  const char *parts[MAX_PARTS], *type, *slots = "";
  int nparts, is_static;

  if (decl->kind == AST_EXTENSION_DECLARATION)
    decl = decl->child->next;
  if (decl->kind != AST_DECLARATION)
    return item;
  specs = decl->child;
  list = specs->next;
  if (list->kind != AST_INIT_DECLARATOR_LIST || !(type = scalar_type (specs))
      || has_token (specs, T_TYPEDEF, NULL) || has_token (specs, T_EXTERN, NULL)
      || has_token (specs, T_REGISTER, NULL) || has_token (specs, T_THREAD_LOCAL, NULL))
    return item;
  nparts = scalar_parts (specs, parts);
  is_static = has_token (specs, T_STATIC, NULL);
  for (struct ast *d = list->child, *next; d; d = next) {
    const char *name = ast_declared_name (d);

    next = d->next;
    if (!plain_name (d))
      continue;
    if (is_static) {
      /* static vector s; static scalar *const slot ... = &s.x; ... */
      for (int k = 0; k < nparts; k++)
        slots = arena_printf ("%s static scalar *const __cwrap_local_slot_%d "
                              "__attribute__ ((section (\"scalars\"), used)) = &%s%s;",
                              slots, local_slot++, name, parts[k]);
    } else {
      /* scalar x __attribute__ ((cleanup (...))) = scalar_new (); and the
       * same with vector_new () or tensor_new () */
      struct ast *item = parse_block_item (arena_printf (
        "%s __v __attribute__ ((cleanup (__cwrap_%s_free))) = %s_new ();",
        type, type, type));
      struct ast *init = item->child->next->child;
      struct ast *placeholder = init->child->child;

      ast_replace (d, init);
      ast_replace (placeholder, d);
    }
  }
  if (!*slots)
    return item;
  /* The slots follow the declaration, on its line, as items of their own. */
  with_slots = ast_new (AST_BLOCK_ITEMS);
  ast_add (with_slots, item);
  for (struct ast *n = parse_snippet (slots), *next; n; n = next) {
    next = n->next;
    n->next = NULL;
    for (struct ast *d = list->child; d; d = d->next)
      if (plain_name (d))
        bind_spelled (n, ast_first_token (d));
    ast_add (with_slots, n);
  }
  return with_slots;
}

/* Each `scalar a;` at file scope, headers included, gets a slot in the section
 * scalars, on the same line so that no line moves:
 *
 *   scalar a; static scalar *const __cwrap_scalar_slot_0 ... = &a;
 *
 * A vector or tensor gets a slot for each of its scalars, &t.x.x, &t.x.y, ...
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
      const char *parts[MAX_PARTS], *slots = "";
      int nparts = scalar_parts (specs, parts);

      if (d->kind == AST_INIT_DECLARATOR && plain_name (d->child))
        warn_scalar_initializer (d, scalar_type (specs));
      if (!plain_name (d))
        continue;
      /* one slot for each scalar: &t.x.x, &t.x.y, ... for a tensor */
      for (int k = 0; k < nparts; k++)
        slots = arena_printf ("%s static scalar *const __cwrap_scalar_slot_%d "
                              "__attribute__ ((section (\"scalars\"), used)) = &%s%s;",
                              slots, slot++, ast_declared_name (d), parts[k]);
      insert_all_before (next, parse_snippet (slots));
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

/* Block macros: @loop */

/* For the names of hidden locals: each expansion has its own. */
static int loop_expansions;

static int count_newlines (const char *s)
{
  int n = 0;

  for (; *s; s++)
    n += *s == '\n';
  return n;
}

/* Puts a linemarker in front of t's trivia, so that t lands on its own line
 * of its own file wherever it is emitted, and the tokens after it with it.
 * Trivia that has a linemarker of its own needs none. */
/* Does trivia hold a linemarker, `# 12 "file"`? A #pragma line does not
 * count: it sets no line. */
static int has_linemarker (const char *tr)
{
  for (const char *p = tr; p; p = strchr (p, '\n')) {
    if (*p == '\n')
      p++;
    p += strspn (p, " \t");
    if (*p == '#') {
      p += 1 + strspn (p + 1, " \t");
      if (*p >= '0' && *p <= '9')
        return 1;
    }
  }
  return 0;
}

static void mark_token (struct token *t)
{
  const char *tr = t->trivia;
  int line;

  if (has_linemarker (tr))
    return;
  line = t->line - count_newlines (tr);
  t->trivia = arena_printf ("\n# %d \"%s\"\n%s", line > 0 ? line : 1,
                            t->file, tr);
}

/* The last token under n. */
static struct token *last_token (struct ast *n)
{
  while (n->kind != AST_TOKEN)
    n = n->last;
  return &((struct ast_leaf *)n)->tok;
}

/* The first token after n, in the order of emission, or NULL. */
static struct token *token_after (struct ast *n)
{
  while (n && !n->next)
    n = n->parent;
  return n ? first_token (n->next) : NULL;
}

static void each_leaf (struct ast *n, void (*f) (struct ast *, void *), void *arg)
{
  if (n->kind == AST_TOKEN) {
    f (n, arg);
    return;
  }
  for (struct ast *k = n->child; k; k = k->next)
    each_leaf (k, f, arg);
}

/* The trivia of leaf, for the stub of a definition: line breaks and
 * linemarkers, so that no line moves, but no other directive. A #pragma in a
 * body, as from _Pragma, belongs to the code where the loop is used, and has
 * no place where it is defined. */
static void add_trivia (struct ast *leaf, void *out)
{
  const char *p = ast_token (leaf)->trivia;

  while (*p) {
    const char *eol = strchr (p, '\n');
    size_t n = eol ? (size_t)(eol - p) + 1 : strlen (p);
    const char *q = p + strspn (p, " \t");

    const char *after = q + 1 + strspn (q + 1, " \t");

    if (*q == '#' && !(*after >= '0' && *after <= '9'))
      fputc ('\n', out);         /* a directive other than a linemarker */
    else
      fwrite (p, 1, n, out);
    p += n;
  }
}

struct ast *loop_stub (const struct ast *source)
{
  MEMSTREAM (trivia);
  struct token t = *ast_first_token (source);

  if (!trivia)
    abort ();
  each_leaf ((struct ast *)source, add_trivia, trivia);
  t.kind = T_IDENTIFIER;
  t.text = "";
  t.trivia = MEMSTREAM_TEXT (trivia);
  t.binding = NULL;
  struct ast *leaf = ast_leaf (&t);
  return ast_node (AST_LOOP_DEFINITION, &leaf, 1);
}

static void add_text (struct ast *leaf, void *out)
{
  fprintf (out, "%s ", ast_token (leaf)->text);
}

int loop_same (const struct loop_def *a, const struct loop_def *b)
{
  MEMSTREAM (ta);
  MEMSTREAM (tb);

  if (!ta || !tb)
    abort ();
  each_leaf (a->source, add_text, ta);
  each_leaf (b->source, add_text, tb);
  return strcmp (MEMSTREAM_TEXT (ta), MEMSTREAM_TEXT (tb)) == 0;
}

/* The parameters of def, in order, into ps; returns how many. (void) is none. */
static int loop_params (const struct loop_def *def, struct ast **ps, int max)
{
  struct ast *list = def->params;
  int n = 0;

  if (!list)
    return 0;
  if (list->kind == AST_PARAMETER_TYPE_LIST)
    list = list->child;
  for (struct ast *p = list->child; p; p = p->next)
    if (p->kind == AST_PARAMETER && n < max)
      ps[n++] = p;
  if (n == 1 && !ps[0]->child->next && has_token (ps[0]->child, T_VOID, NULL))
    return 0;
  return n;
}

#define MAX_LOOP_PARAMS 64

int loop_check_params (const struct loop_def *def)
{
  struct ast *ps[MAX_LOOP_PARAMS];
  int n = loop_params (def, ps, MAX_LOOP_PARAMS), errors = 0;

  if (def->params && def->params->kind == AST_PARAMETER_TYPE_LIST) {
    diag (def->at->file, def->at->line, DIAG_ERROR, "loop '%s' takes a fixed "
          "number of arguments: no '...'", def->name);
    errors++;
  }
  for (int i = 0; i < n; i++)
    if (!ast_declared_name (ps[i])) {
      const struct token *t = ast_first_token (ps[i]);

      diag (t->file, t->line, DIAG_ERROR, "parameter %d of loop '%s' has no "
            "name", i + 1, def->name);
      errors++;
    }
  return errors;
}

static int injected (const struct loop_def *def, const char *name)
{
  for (int i = 0; i < def->ninject; i++)
    if (strcmp (def->inject[i].name, name) == 0)
      return 1;
  return 0;
}

/* Is leaf, a name, one that may stand for a variable or typedef here? Not a
 * member, designator, label, tag or attribute, which share the spelling of a
 * name but not its meaning. */
static int names_ordinary (const struct ast *leaf)
{
  const struct ast *p = leaf->parent;

  switch (p->kind) {
  case AST_MEMBER:
    return leaf != p->last;
  case AST_DESIGNATOR: case AST_OLD_DESIGNATOR: case AST_OFFSETOF_MEMBER:
  case AST_NAMED_LABEL: case AST_GOTO_STATEMENT: case AST_LABEL_LIST:
  case AST_LABEL_ADDRESS: case AST_STRUCT_SPECIFIER: case AST_ENUM_SPECIFIER:
  case AST_ATTRIBUTE_ITEM: case AST_ASM_LABELS: case AST_ASM_OPERAND:
    return 0;
  case AST_DECL_NAME:
    /* a struct member's declarator */
    for (p = p->parent; p; p = p->parent) {
      if (p->kind == AST_STRUCT_DECLARATOR)
        return 0;
      if (p->kind != AST_DECLARATOR && p->kind != AST_PAREN_DECLARATOR
          && p->kind != AST_SUFFIXED_DECLARATOR
          && p->kind != AST_ATTRIBUTED_DECLARATOR)
        return 1;
    }
    return 1;
  default:
    return 1;
  }
}

struct rename {
  const struct loop_def *def;
  int expansion;
  const struct var *const *given;
};

/* The body's own names, where declared and where used, get a name of their
 * own for this expansion, unless the block is to see them. */
static void rename_leaf (struct ast *leaf, void *arg)
{
  const struct rename *r = arg;
  struct token *t = &((struct ast_leaf *)leaf)->tok;

  /* a name the block sees: the declaration the block saw stands for it */
  for (int i = 0; t->binding && i < r->def->ninject; i++)
    if (t->binding == r->def->inject[i].var) {
      t->binding = r->given[i];
      return;
    }
  if ((t->kind != T_IDENTIFIER && t->kind != T_TYPE_NAME)
      || !t->binding || t->binding->loop != r->def->serial
      || injected (r->def, t->text) || !names_ordinary (leaf))
    return;
  t->text = arena_printf ("__cwrap_%s_%s_%d", r->def->name, t->text,
                          r->expansion);
}

/* One space before each token, at the place of `at`: code that cwrap adds
 * on the line of the use. */
static void flatten_leaf (struct ast *leaf, void *at)
{
  struct token *t = &((struct ast_leaf *)leaf)->tok;

  t->trivia = " ";
  t->file = ((const struct token *)at)->file;
  t->line = ((const struct token *)at)->line;
}

static struct ast *new_token (int kind, const char *text, const struct token *at)
{
  struct token t = *at;

  t.kind = kind;
  t.text = text;
  t.trivia = " ";
  t.binding = NULL;
  return ast_leaf (&t);
}

/* The DECL_NAME token under a declarator. */
static struct token *declared_token (struct ast *d)
{
  if (d->kind == AST_DECL_NAME)
    return &((struct ast_leaf *)d->child)->tok;
  for (d = d->child; d; d = d->next)
    if (d->kind != AST_TOKEN) {
      struct token *t = declared_token (d);

      if (t)
        return t;
    }
  return NULL;
}

/* specs declarator = value; from a parameter of the loop, renamed for this
 * expansion. If name is not NULL, the declarator declares it instead. */
static struct ast *bind_parameter (struct ast *param, const char *name,
                                   struct ast *value, const struct rename *r,
                                   const struct token *at)
{
  struct ast *specs = ast_copy (param->child);
  struct ast *declarator = ast_copy (param->child->next);
  struct ast *init, *list, *kids[3];

  each_leaf (specs, rename_leaf, (void *)r);
  each_leaf (declarator, rename_leaf, (void *)r);
  if (name)
    declared_token (declarator)->text = name;
  kids[0] = declarator;
  kids[1] = new_token ('=', "=", at);
  kids[2] = value;
  init = ast_node (AST_INIT_DECLARATOR, kids, 3);
  list = ast_node (AST_INIT_DECLARATOR_LIST, &init, 1);
  kids[0] = specs;
  kids[1] = list;
  kids[2] = new_token (';', ";", at);
  each_leaf (specs, flatten_leaf, (void *)at);
  each_leaf (declarator, flatten_leaf, (void *)at);
  return ast_node (AST_DECLARATION, kids, 3);
}

static struct ast *find_kind (struct ast *n, enum ast_kind kind)
{
  if (n->kind == kind)
    return n;
  if (n->kind == AST_TOKEN)
    return NULL;
  for (n = n->child; n; n = n->next) {
    struct ast *found = find_kind (n, kind);

    if (found)
      return found;
  }
  return NULL;
}

struct ast *loop_expand (const struct loop_def *def, struct ast *name,
                         struct ast *args, struct ast *block,
                         const struct var *const *given, int *rejected)
{
  const struct token *at = ast_token (name);
  struct ast *ps[MAX_LOOP_PARAMS], *vs[MAX_LOOP_PARAMS];
  int np = loop_params (def, ps, MAX_LOOP_PARAMS), nv = 0;
  struct rename r = { def, ++loop_expansions, given };
  struct ast *items, *body, *placeholder, *open, *close, *kids[3];
  struct token *after, *end;

  for (struct ast *a = args ? args->child : NULL; a; a = a->next)
    if (!(a->kind == AST_TOKEN && strcmp (ast_token (a)->text, ",") == 0)
        && nv < MAX_LOOP_PARAMS)
      vs[nv++] = a;
  if (nv != np) {
    diag (at->file, at->line, DIAG_ERROR, "loop '%s' takes %d argument%s, not "
          "%d", def->name, np, np == 1 ? "" : "s", nv);
    diag (def->at->file, def->at->line, DIAG_NOTE, "'%s' is defined here",
          def->name);
    (*rejected)++;
    return block;
  }
  for (int i = 0; i < nv; i++) {
    vs[i]->parent = vs[i]->next = NULL;
    each_leaf (vs[i], flatten_leaf, (void *)at);
  }

  /* The arguments are all evaluated before a parameter the block sees comes
   * into scope, as `T p = p;` would read the new p: those go through a
   * hidden local first. */
  items = ast_new (AST_BLOCK_ITEMS);
  for (int i = 0; i < np; i++) {
    const char *pname = ast_declared_name (ps[i]);

    if (injected (def, pname)) {
      const char *tmp = arena_printf ("__cwrap_%s_arg%d_%d", def->name, i,
                                      r.expansion);
      ast_add (items, bind_parameter (ps[i], tmp, vs[i], &r, at));
      vs[i] = new_token (T_IDENTIFIER, tmp, at);
    } else {
      ast_add (items, bind_parameter (ps[i], NULL, vs[i], &r, at));
      vs[i] = NULL;
    }
  }
  for (int i = 0; i < np; i++)
    if (vs[i])
      ast_add (items, bind_parameter (ps[i], NULL, vs[i], &r, at));

  /* The body, on the lines of the definition, and the block on its own. */
  body = ast_copy (def->body);
  each_leaf (body, rename_leaf, &r);
  placeholder = find_kind (body, AST_LOOP_BODY);
  after = token_after (placeholder);
  if (after)
    mark_token (after);
  mark_token (first_token (body));
  mark_token (first_token (block));
  end = last_token (block);
  ast_replace (placeholder, block);
  ast_add (items, body);

  /* { at the place of the name, } on the line where the block ends, so the
   * code after it keeps its lines. */
  open = new_token ('{', "{", at);
  first_token (open)->trivia = at->trivia;
  close = new_token ('}', "}", end);
  first_token (close)->trivia = arena_printf ("\n# %d \"%s\"\n", end->line,
                                              end->file);
  kids[0] = open;
  kids[1] = items;
  kids[2] = close;
  return ast_node (AST_COMPOUND_STATEMENT, kids, 3);
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
