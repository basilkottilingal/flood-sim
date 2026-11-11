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
  int hide;               /* the names given to the block are renamed too */
};

/* The body's own names, where declared and where used, get a name of their
 * own for this expansion, unless the block is to see them. */
static void rename_leaf (struct ast *leaf, void *arg)
{
  const struct rename *r = arg;
  struct token *t = &((struct ast_leaf *)leaf)->tok;
  int given = 0;

  /* a name the block sees: the declaration the block saw stands for it. A
   * hidden argument of a loop used in the body (`__cwrap_..._arg0_1`) keeps
   * its name either way. */
  for (int i = 0; t->binding && i < r->def->ninject; i++)
    if (t->binding == r->def->inject[i].var) {
      if (r->hide && strcmp (t->text, r->def->inject[i].name) == 0) {
        given = 1;
        break;
      }
      if (r->given)
        t->binding = r->given[i];
      return;
    }
  if ((t->kind != T_IDENTIFIER && t->kind != T_TYPE_NAME)
      || !t->binding || (!given && t->binding->loop != r->def->serial)
      || (!r->hide && injected (r->def, t->text)) || !names_ordinary (leaf))
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

/* With hide, the block sees none of the loop's names: they are renamed like
 * the rest, to __cwrap_<loop>_<name>_<expansion>. */
static struct ast *expand_loop (const struct loop_def *def, struct ast *name,
                                struct ast *args, struct ast *block,
                                const struct var *const *given, int hide,
                                int *rejected)
{
  const struct token *at = ast_token (name);
  struct ast *ps[MAX_LOOP_PARAMS], *vs[MAX_LOOP_PARAMS];
  int np = loop_params (def, ps, MAX_LOOP_PARAMS), nv = 0;
  struct rename r = { def, ++loop_expansions, given, hide };
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

struct ast *loop_expand (const struct loop_def *def, struct ast *name,
                         struct ast *args, struct ast *block,
                         const struct var *const *given, int *rejected)
{
  return expand_loop (def, name, args, block, given, 0, rejected);
}

/* A global scalar is never freed */

static const struct ast *declarator_specs (const struct ast *d);
static int declarator_indirections (const struct ast *d);

/* The declaration of the global scalar, vector or tensor that e names, or
 * of which it names a member (t.x.y), or NULL: one at file scope, static or
 * extern. *root is then the name, *member whether e is a member. */
static const struct ast *global_scalar_decl (const struct ast *e,
                                             const struct ast **specs_out,
                                             int *member)
{
  const struct ast *specs, *decl, *top = e;

  for (;;) {
    if (e->kind == AST_PAREN_EXPRESSION)
      e = e->child->next;
    else if (e->kind == AST_MEMBER && strcmp (ast_token (e->child->next)->text, ".") == 0)
      e = e->child;
    else
      break;
  }
  if (e->kind != AST_TOKEN || !ast_token (e)->binding
      || ast_token (e)->binding->kind != SCOPE_VARIABLE
      || !(decl = ast_token (e)->binding->decl)
      || !(specs = declarator_specs (decl)) || !scalar_type (specs)
      || declarator_indirections (decl) != 0)
    return NULL;
  {
    const struct ast *d = specs->parent;  /* the declaration */

    if (d->parent && d->parent->kind == AST_EXTENSION_DECLARATION)
      d = d->parent;
    if (!has_token (specs, T_STATIC, NULL) && !has_token (specs, T_EXTERN, NULL)
        && !(d->parent && d->parent->kind == AST_TRANSLATION_UNIT))
      return NULL;
  }
  *specs_out = specs;
  *member = e != top;
  return decl;
}

int check_global_free (const struct ast *fn, const struct ast *args)
{
  const struct ast *specs, *decl;
  const struct token *t;
  int member;

  if (!dimension || fn->kind != AST_TOKEN || !args || !args->child)
    return 0;
  t = ast_token (fn);
  if (t->kind != T_IDENTIFIER || (strcmp (t->text, "scalar_free") != 0
                                  && strcmp (t->text, "vector_free") != 0
                                  && strcmp (t->text, "tensor_free") != 0))
    return 0;
  if (!(decl = global_scalar_decl (args->child, &specs, &member)))
    return 0;
  diag (t->file, t->line, DIAG_ERROR, "%s () of '%s', %s global %s: it is "
        "numbered once for the whole run and never freed", t->text,
        expression_text (args->child), member ? "a member of a" : "a",
        scalar_type (specs));
  note_declared (decl);
  return 1;
}

/* Stencils: foreach () and s[...] */

int stencil_enabled (void)
{
  return dimension > 0;
}

/* The pointers and arrays a declarator wraps around its name: how many
 * subscripts or `*` go from the variable to its type. -1 for a function. */
static int declarator_indirections (const struct ast *d)
{
  int n = 0;

  switch (d->kind) {
  case AST_DECL_NAME:
    return 0;
  case AST_POINTER:
    for (const struct ast *k = d->child; k; k = k->next)
      if (k->kind == AST_TOKEN)
        n++;                              /* '*' */
      else if (k->kind == AST_POINTER)
        n += declarator_indirections (k);
    return n;
  case AST_DECLARATOR: case AST_PAREN_DECLARATOR: case AST_SUFFIXED_DECLARATOR:
  case AST_ATTRIBUTED_DECLARATOR:
    for (const struct ast *k = d->child; k; k = k->next) {
      int m;

      switch (k->kind) {
      case AST_ARRAY_SUFFIX:
        n++;
        break;
      case AST_FUNCTION_SUFFIX:
        return -1;
      case AST_DECL_NAME: case AST_POINTER: case AST_DECLARATOR:
      case AST_PAREN_DECLARATOR: case AST_SUFFIXED_DECLARATOR:
        if ((m = declarator_indirections (k)) < 0)
          return -1;
        n += m;
        break;
      default:
        break;
      }
    }
    return n;
  default:
    return -1;
  }
}

/* The specifiers of the declaration or parameter that declares d, whatever
 * its declarator, or NULL. */
static const struct ast *declarator_specs (const struct ast *d)
{
  for (const struct ast *n = d->parent; n; n = n->parent) {
    if (n->kind == AST_DECLARATION || n->kind == AST_PARAMETER)
      return n->child;
    if (n->kind != AST_ATTRIBUTED_DECLARATOR && n->kind != AST_INIT_DECLARATOR
        && n->kind != AST_INIT_DECLARATOR_LIST)
      return NULL;
  }
  return NULL;
}

/* What e is, of cwrap's types: *depth 0 for a scalar, 1 a vector, 2 a tensor,
 * and *ind the subscripts or `*` still to go, for an array of them or a
 * pointer. Returns 0 if e is none of these, or not an expression cwrap can
 * follow: a name, (e), e.x, e->x, e[i] and *e. */
static int scalar_shape (const struct ast *e, int *depth, int *ind)
{
  const struct token *t;
  const struct ast *specs;
  const char *type;

  switch (e->kind) {
  case AST_TOKEN:
    t = ast_token (e);
    if (t->kind != T_IDENTIFIER && t->kind != T_FOREACH && t->kind != T_CACHE)
      return 0;
    if (!t->binding || t->binding->kind != SCOPE_VARIABLE || !t->binding->decl
        || !(specs = declarator_specs (t->binding->decl))
        || !(type = scalar_type (specs)) || has_token (specs, T_TYPEDEF, NULL)
        || (*ind = declarator_indirections (t->binding->decl)) < 0)
      return 0;
    *depth = strcmp (type, "scalar") == 0 ? 0 : strcmp (type, "vector") == 0 ? 1 : 2;
    return 1;
  case AST_PAREN_EXPRESSION:
    return scalar_shape (e->child->next, depth, ind);
  case AST_MEMBER:
    if (!scalar_shape (e->child, depth, ind) || *depth == 0)
      return 0;
    if (strcmp (ast_token (e->child->next)->text, "->") == 0) {
      if (*ind != 1)
        return 0;
      *ind = 0;
    } else if (*ind != 0)
      return 0;
    --*depth;
    return 1;
  case AST_SUBSCRIPT:
    if (!scalar_shape (e->child, depth, ind) || *ind == 0)
      return 0;
    --*ind;
    return 1;
  case AST_UNARY_OP:
    if (strcmp (ast_token (e->child)->text, "*") != 0
        || !scalar_shape (e->child->next, depth, ind) || *ind == 0)
      return 0;
    --*ind;
    return 1;
  default:
    return 0;
  }
}

static int is_scalar_expression (const struct ast *e)
{
  int depth, ind;

  return scalar_shape (e, &depth, &ind) && depth == 0 && ind == 0;
}

/* The offsets of s[a, b, c], a comma expression, into os; returns how many. */
static int stencil_offsets (struct ast *e, struct ast **os, int max)
{
  int n = 0;

  if (!e)
    return 0;
  if (e->kind == AST_COMMA)
    n = stencil_offsets (e->child, os, max), e = e->last;
  if (n < max)
    os[n] = e;
  return n + 1;
}

static void stencil_error (int *rejected, const struct ast *at,
                           const char *format, ...)
{
  const struct token *t = ast_first_token (at);
  char *msg;
  va_list ap;
  int n;

  va_start (ap, format);
  n = vsnprintf (NULL, 0, format, ap);
  va_end (ap);
  msg = allocate ((size_t)n + 1);
  va_start (ap, format);
  vsnprintf (msg, (size_t)n + 1, format, ap);
  va_end (ap);
  diag (t->file, t->line, DIAG_ERROR, "%s", msg);
  (*rejected)++;
}

struct stencil_loop *stencil_loop_begin (struct stencil_loop *outer,
                                         struct ast *kw, int *rejected)
{
  struct stencil_loop *loop = allocate (sizeof (*loop));

  loop->name = ast_token (kw)->text;
  if (outer) {
    stencil_error (rejected, kw, "%s () inside %s (): stencil loops do not "
                   "nest", loop->name, outer->name);
    diag (outer->at->file, outer->at->line, DIAG_NOTE, "the outer %s () is "
          "here", outer->name);
  }
  loop->outer = outer;
  loop->at = ast_token (kw);
  loop->errors = *rejected;
  return loop;
}

void stencil_mover_begin (struct stencil_loop *loop, struct ast *kw,
                          int *rejected)
{
  const char *name = ast_token (kw)->text;

  if (!loop) {
    stencil_error (rejected, kw, "%s () outside a stencil loop, whose point "
                   "it moves", name);
    return;
  }
  if (loop->in_mover)
    stencil_error (rejected, kw, "%s () inside %s (): a point is moved once",
                   name, loop->mover_name);
  else
    loop->mover_name = name;
  loop->in_mover++;
}

struct ast *stencil_mover_end (struct stencil_loop *loop, struct ast *kw,
                               struct ast *open, struct ast *radius,
                               struct ast *close, struct ast *body)
{
  struct ast *node = ast_new (AST_STENCIL_MOVER);

  ast_add (ast_add (ast_add (ast_add (ast_add (node, kw), open), radius), close),
           body);
  if (!loop)
    return node;
  loop->in_mover--;
  if (loop->nmover == loop->mcap) {
    struct ast **m = allocate ((size_t)(loop->mcap = 2 * loop->mcap + 4) * sizeof (*m));

    if (loop->nmover)
      memcpy (m, loop->mover, (size_t)loop->nmover * sizeof (*m));
    loop->mover = m;
  }
  loop->mover[loop->nmover++] = node;
  return node;
}

struct ast *stencil_coarse (struct stencil_loop *loop, struct ast *fn,
                            struct ast *open, struct ast *args,
                            struct ast *close, int *rejected)
{
  struct ast *node, *scalar, *offsets = NULL, *a;
  const struct token *t;
  int n = 0;

  if (!loop || fn->kind != AST_TOKEN)
    return NULL;
  t = ast_token (fn);
  if (t->kind != T_IDENTIFIER || strcmp (t->text, "coarse") != 0 || t->binding
      || !(scalar = args->child) || !is_scalar_expression (scalar))
    return NULL;
  /* coarse (s, i, j) as s(i, j): the scalar, then the offsets, a comma
   * expression, as in s[i, j] */
  for (a = scalar->next; a; a = a->next)
    n += !(a->kind == AST_TOKEN && strcmp (ast_token (a)->text, ",") == 0);
  if (n > dimension)
    stencil_error (rejected, fn, "'%s' has %d offsets, in %d dimension%s",
                   expression_text (args), n, dimension,
                   dimension == 1 ? "" : "s");
  for (a = scalar->next; a; ) {
    struct ast *next = a->next;

    if (a->kind != AST_TOKEN || strcmp (ast_token (a)->text, ",") != 0) {
      a->parent = a->next = NULL;
      if (!offsets)
        offsets = a;
      else {
        struct ast *kids[3] = { offsets, new_token (',', ",", ast_first_token (a)), a };

        offsets = ast_node (AST_COMMA, kids, 3);
      }
    }
    a = next;
  }
  scalar->parent = scalar->next = NULL;
  first_token (scalar)->trivia = arena_printf ("%s%s", t->trivia,
                                               first_token (scalar)->trivia);
  node = ast_new (AST_STENCIL_COARSE);
  ast_add (ast_add (ast_add (ast_add (node, scalar), open), offsets), close);
  if (loop->naccess == loop->cap) {
    struct ast **m = allocate ((size_t)(loop->cap = 2 * loop->cap + 8) * sizeof (*m));

    if (loop->naccess)
      memcpy (m, loop->access, (size_t)loop->naccess * sizeof (*m));
    loop->access = m;
  }
  loop->access[loop->naccess++] = node;
  return node;
}

/* How an access reads in a message: s[1,0], or coarse (s, 1) */
static const char *access_text (const struct ast *access)
{
  if (access->kind != AST_STENCIL_COARSE)
    return expression_text (access);
  const struct ast *off = access->child->next->next;

  return off == access->last
    ? arena_printf ("coarse (%s)", expression_text (access->child))
    : arena_printf ("coarse (%s, %s)", expression_text (access->child),
                    expression_text (off));
}

struct ast *stencil_access (struct stencil_loop *loop, struct ast *base,
                            struct ast *open, struct ast *offsets,
                            struct ast *close, int *rejected)
{
  struct ast *os[4], *node;
  int n;

  if (!is_scalar_expression (base)) {
    if (offsets)
      return NULL;
    stencil_error (rejected, base, "'%s[]', but '%s' is not a scalar",
                   expression_text (base), expression_text (base));
    node = ast_new (AST_SUBSCRIPT);
    ast_add (ast_add (ast_add (node, base), open), close);
    return node;
  }
  node = ast_new (AST_STENCIL_ACCESS);
  ast_add (ast_add (ast_add (ast_add (node, base), open), offsets), close);
  if (!loop) {
    stencil_error (rejected, base, "'%s[...]', a cell of scalar '%s', outside "
                   "foreach ()", expression_text (base), expression_text (base));
    return node;
  }
  n = stencil_offsets (offsets, os, 4);
  if (n > dimension)
    stencil_error (rejected, base, "'%s' has %d offsets, in %d dimension%s",
                   expression_text (node), n, dimension, dimension == 1 ? "" : "s");
  if (loop->naccess == loop->cap) {
    struct ast **a = allocate ((size_t)(loop->cap = 2 * loop->cap + 8) * sizeof (*a));

    if (loop->naccess)
      memcpy (a, loop->access, (size_t)loop->naccess * sizeof (*a));
    loop->access = a;
  }
  loop->access[loop->naccess++] = node;
  return node;
}

/* Is n under root? */
static int is_under (const struct ast *n, const struct ast *root)
{
  for (; n; n = n->parent)
    if (n == root)
      return 1;
  return 0;
}

/* The variable whose storage e changes, in e = ..., e++, &e: through (e),
 * members and subscripts down to its name. NULL if it is not a name. */
static const struct var *changed_variable (const struct ast *e)
{
  for (;;) {
    if (e->kind == AST_PAREN_EXPRESSION)
      e = e->child->next;
    else if (e->kind == AST_MEMBER || e->kind == AST_SUBSCRIPT)
      e = e->child;
    else
      break;
  }
  return e->kind == AST_TOKEN ? ast_token (e)->binding : NULL;
}

/* The operand that n changes, if n is an assignment, ++, -- or &; else NULL. */
static const struct ast *changed_operand (const struct ast *n)
{
  switch (n->kind) {
  case AST_ASSIGNMENT: case AST_POSTFIX_OP:
    return n->child;
  case AST_PREFIX_OP:
    return n->child->next;
  case AST_UNARY_OP:
    return strcmp (ast_token (n->child)->text, "&") == 0 ? n->child->next : NULL;
  default:
    return NULL;
  }
}

/* Does anything under n change v? */
static int changes_variable (const struct ast *n, const struct var *v)
{
  const struct ast *target;

  if (n->kind == AST_TOKEN)
    return 0;
  if ((target = changed_operand (n)) && changed_variable (target) == v)
    return 1;
  for (n = n->child; n; n = n->next)
    if (changes_variable (n, v))
      return 1;
  return 0;
}

/* Can the block of scalar e be found once per leaf, before the cells? Only if
 * e is made of names declared outside the body, which it does not change, and
 * of constants: s, u.x, list[2], list[n]. */
static int stencil_invariant (const struct ast *e, const struct ast *body)
{
  if (e->kind == AST_TOKEN) {
    const struct token *t = ast_token (e);

    if (t->kind == T_CONSTANT)
      return 1;
    if (t->kind != T_IDENTIFIER && t->kind != T_FOREACH && t->kind != T_CACHE)
      return strchr ("()[].*+-", t->text[0]) != NULL
             || strcmp (t->text, "->") == 0;
    return t->binding && t->binding->kind == SCOPE_VARIABLE && t->binding->decl
           && !is_under (t->binding->decl, body)
           && !changes_variable (body, t->binding);
  }
  switch (e->kind) {
  case AST_MEMBER:                        /* the member's name is no variable */
    return stencil_invariant (e->child, body);
  case AST_PAREN_EXPRESSION: case AST_SUBSCRIPT: case AST_UNARY_OP: case AST_ADDITIVE: case AST_MULTIPLICATIVE:
    break;
  default:
    return 0;
  }
  for (e = e->child; e; e = e->next)
    if (!stencil_invariant (e, body))
      return 0;
  return 1;
}

/* Does the access write its cell: s[] = ..., s[] += ..., s[]++, &s[]? */
static int stencil_writes (const struct ast *access)
{
  const struct ast *n = access, *p;

  while ((p = n->parent) && p->kind == AST_PAREN_EXPRESSION)
    n = p;
  return p && changed_operand (p) == n;
}

/* The first leaf under n spelled text, or NULL. */
static struct ast *find_leaf (struct ast *n, const char *text)
{
  if (n->kind == AST_TOKEN)
    return strcmp (ast_token (n)->text, text) == 0 ? n : NULL;
  for (n = n->child; n; n = n->next) {
    struct ast *found = find_leaf (n, text);

    if (found)
      return found;
  }
  return NULL;
}

/* Puts e where the placeholder leaf `name` is under n. */
static void fill_placeholder (struct ast *n, const char *name, struct ast *e)
{
  struct ast *leaf = find_leaf (n, name);

  if (!leaf)
    abort ();
  e->parent = e->next = NULL;
  ast_replace (leaf, e);
}

/* The item text, a statement or declaration, as cwrap adds it to a body, with
 * the type Real known. */
static struct ast *parse_stencil_item (const char *text)
{
  struct ast *fn = parse_snippet (arena_printf ("typedef int Real; "
                                                "void f (void) { %s }", text))->next;
  struct ast *item = fn->child->next->next->child->next->child;

  item->parent = item->next = NULL;
  return item;
}

/* A scalar whose block is found once per leaf. */
struct stencil_block {
  const char *key;          /* the scalar as written: s, u.x */
  struct ast *scalar;       /* that expression, from its first access */
  int writes;
  int level;                /* 0 the point's, 1 the coarse, 2 a child's */
};

#define MAX_STENCIL_BLOCKS 256

/* The names that the expansions of the two loops give what the body needs:
 * see rename_leaf (). */
#define STENCIL_NAME(loop, name, expansion) \
  arena_printf ("__cwrap_%s_%s_%d", loop, name, expansion)

/* Has def a name for its block called name? */
static int loop_gives (const struct loop_def *def, const char *name)
{
  return injected (def, name);
}

static int stencil_serial;

static void add_text (struct ast *leaf, void *out);

/* What a stencil loop expands to, marked as one, for check_stencil_calls (). */
static struct ast *stencil_mark (struct ast *n)
{
  return ast_node (AST_STENCIL_LOOP, &n, 1);
}

/* A copy of args, a declaration's arguments to one of its loops, with each
 * of its parameters ps [0..n) named names [k], or NULL. */
static void rename_params (struct ast *leaf, void *arg)
{
  struct ast **pass = arg;
  struct ast **ps = (struct ast **)pass[0];
  const char **names = (const char **)pass[1];
  int n = (int)(long)pass[2];
  struct token *t = &((struct ast_leaf *)leaf)->tok;

  if (t->kind != T_IDENTIFIER || !names_ordinary (leaf))
    return;
  for (int k = 0; k < n; k++)
    if (strcmp (t->text, ast_declared_name (ps[k])) == 0) {
      t->text = names[k];
      t->binding = NULL;
      return;
    }
}

/* In the arguments of a stencil loop's cell loop, the names its block loop
 * gives, as child, are its hidden ones, __cwrap_<loop>_<name>_<expansion>. */
struct given_names {
  const struct loop_def *def;
  int expansion;
};

static void rename_given (struct ast *leaf, void *arg)
{
  const struct given_names *g = arg;
  struct token *t = &((struct ast_leaf *)leaf)->tok;

  if (t->kind == T_IDENTIFIER && names_ordinary (leaf) && injected (g->def, t->text))
    t->text = arena_printf ("__cwrap_%s_%s_%d", g->def->name, t->text,
                            g->expansion);
}

static struct ast *stencil_arguments (const struct ast *args, struct ast **ps,
                                      const char **names, int n)
{
  struct ast *copy;
  void *pass[3] = { ps, names, (void *)(long)n };

  if (!args)
    return NULL;
  copy = ast_copy (args);
  each_leaf (copy, rename_params, pass);
  return copy;
}

/* The tokens under n, one space apart, to compare two definitions. */
static const char *source_text (const struct ast *n)
{
  MEMSTREAM (t);

  if (!t)
    abort ();
  each_leaf ((struct ast *)n, add_text, t);
  return MEMSTREAM_TEXT (t);
}

struct ast *stencil_define (struct stencil_def **list, struct ast *source,
                            struct ast *name, struct ast *params,
                            const struct loop_def *leaf,
                            struct ast *leaf_name, struct ast *leaf_args,
                            const struct loop_def *cell, struct ast *cell_name,
                            struct ast *cell_args, int *rejected)
{
  struct stencil_def *def = allocate (sizeof (*def)), *old;
  struct ast *ps[MAX_LOOP_PARAMS];
  struct loop_def holder = { 0 };
  int n;

  def->name = ast_token (name)->text;
  def->at = ast_token (name);
  def->source = source;
  def->params = params;
  def->leaf = leaf;
  def->cell = cell;
  def->leaf_args = leaf_args;
  def->cell_args = cell_args;
  if (!leaf)
    stencil_error (rejected, leaf_name, "@foreach %s: '%s' is not a @loop",
                   def->name, ast_token (leaf_name)->text);
  if (!cell)
    stencil_error (rejected, cell_name, "@foreach %s: '%s' is not a @loop",
                   def->name, ast_token (cell_name)->text);
  holder.name = def->name;
  holder.at = def->at;
  holder.params = params;
  *rejected += loop_check_params (&holder);
  n = loop_params (&holder, ps, MAX_LOOP_PARAMS);
  def->params = n ? params : NULL;
  for (old = *list; old && strcmp (old->name, def->name) != 0; old = old->next)
    ;
  if (old) {
    if (strcmp (source_text (old->source), source_text (source)) != 0) {
      stencil_error (rejected, name, "@foreach %s is defined again, "
                     "differently", def->name);
      diag (old->at->file, old->at->line, DIAG_NOTE, "'%s' is defined here",
            old->name);
    }
  } else {
    def->next = *list;
    *list = def;
  }
  return loop_stub (source);
}

/* In foreach (), x, y, z and Delta are the cell's: the coordinates of its
 * centre and its size, as in Basilisk. */
static const char *const stencil_reserved[] = { "x", "y", "z", "Delta" };
#define NRESERVED 4

/* Counts in uses[] the names under n that mean the cell's x, y, z or Delta:
 * spelled so, where a variable could be (not u.x), and not declared in the
 * body itself. With clear, their tokens forget the declaration they saw, as
 * they now name cwrap's. Returns how many. */
static int reserved_uses (struct ast *n, const struct ast *body, int uses[],
                          int clear)
{
  int found = 0;

  if (n->kind != AST_TOKEN) {
    for (n = n->child; n; n = n->next)
      found += reserved_uses (n, body, uses, clear);
    return found;
  }
  struct token *t = &((struct ast_leaf *)n)->tok;

  if (t->kind != T_IDENTIFIER || !names_ordinary (n)
      || (t->binding && t->binding->decl && is_under (t->binding->decl, body)))
    return 0;
  for (int r = 0; r < NRESERVED; r++)
    if (strcmp (t->text, stencil_reserved[r]) == 0) {
      uses[r]++;
      if (clear)
        t->binding = NULL;
      return 1;
    }
  return 0;
}

/* An error for each name under n that would be the cell's x, y, z or Delta
 * but was declared outside the body as a scalar, a vector or a tensor (or
 * an array of them, or a pointer): the coordinate would hide it. */
static void reserved_hides_scalar (struct ast *n, const struct ast *body,
                                   int *rejected)
{
  int depth, ind, none[NRESERVED] = { 0 };

  if (n->kind != AST_TOKEN) {
    for (n = n->child; n; n = n->next)
      reserved_hides_scalar (n, body, rejected);
    return;
  }
  const struct token *t = ast_token (n);

  if (!reserved_uses (n, body, none, 0) || !scalar_shape (n, &depth, &ind))
    return;
  stencil_error (rejected, n, "'%s' in foreach () is the cell's %s, which "
                 "hides %s '%s'; rename the %s", t->text,
                 strcmp (t->text, "Delta") == 0 ? "size" : "coordinate",
                 depth == 0 ? "scalar" : depth == 1 ? "vector" : "tensor",
                 t->text, depth == 0 ? "scalar" : depth == 1 ? "vector" : "tensor");
  note_declared (t->binding->decl);
}

/* Does the access reach a cell other than its own: an offset that is not
 * 0 as written? s[i] counts, whatever i is. */
static int stencil_has_offset (struct ast *access)
{
  struct ast *os[4], *off = access->child->next->next;
  int n;

  if (off == access->last)
    return 0;
  n = stencil_offsets (off, os, 4);
  for (int k = 0; k < n && k < 4; k++)
    if (os[k]->kind != AST_TOKEN || ast_token (os[k])->kind != T_CONSTANT
        || strcmp (ast_token (os[k])->text, "0") != 0)
      return 1;
  return 0;
}

/* A scalar whose ghosts a foreach () makes stale, or needs up to date: the
 * scalar itself, if it is the same at every cell; else every element of
 * the array it is one of, A[k] of an array A of a known size; else any
 * scalar at all. */
struct ghost_scalar {
  const char *key;          /* s, A[*] or *, to know two of them the same */
  struct ast *scalar;       /* a copy of the scalar, */
  struct ast *array;        /* or of the array, */
  int all;                  /* or any scalar */
  struct ast *at;           /* the first access, for messages */
};

static struct ghost_scalar ghost_scalar_of (struct ast *access,
                                            const struct ast *body)
{
  struct ast *e = access->child;
  struct ghost_scalar g = { NULL, NULL, NULL, 0, access };

  if (stencil_invariant (e, body)) {
    g.key = expression_text (e);
    g.scalar = ast_copy (e);
    return g;
  }
  if (e->kind == AST_SUBSCRIPT && e->child->kind == AST_TOKEN
      && stencil_invariant (e->child, body)) {
    const struct ast *d = ast_token (e->child)->binding->decl;
    const struct ast *specs = declarator_specs (d);

    if (d->kind == AST_SUFFIXED_DECLARATOR && d->child->kind == AST_DECL_NAME
        && d->last->kind == AST_ARRAY_SUFFIX
        && d->last->child->next != d->last->last      /* [N], not [] */
        && specs && specs->parent->kind == AST_DECLARATION) {
      g.key = arena_printf ("%s[*]", expression_text (e->child));
      g.array = ast_copy (e->child);
      return g;
    }
  }
  g.key = "*";
  g.all = 1;
  return g;
}

/* Adds g to set [0..*n) unless one with its key is there; returns its place. */
static int add_ghost_scalar (struct ghost_scalar *set, int *n,
                             struct ghost_scalar g)
{
  for (int k = 0; k < *n; k++)
    if (strcmp (set[k].key, g.key) == 0)
      return k;
  set[*n] = g;
  return (*n)++;
}

/* Code that cwrap adds as items, text with placeholders __cwrap_e0, ...
 * that take the expressions es [0..), each a copy, all on the line of at. */
static void add_stencil_items (struct ast *items, const char *text,
                               struct ast **es, int nes,
                               const struct token *at)
{
  struct ast *fn = parse_snippet (arena_printf ("typedef int Real; "
                                                "void f (void) { %s }", text))->next;
  struct ast *block = fn->child->next->next, *list = block->child->next;

  for (int k = 0; k < nes; k++) {
    struct ast *e = ast_copy (es[k]);

    each_leaf (e, flatten_leaf, (void *)at);
    fill_placeholder (list, arena_printf ("__cwrap_e%d", k), e);
  }
  for (struct ast *item = list->child, *next; item; item = next) {
    next = item->next;
    item->parent = item->next = NULL;
    each_leaf (item, flatten_leaf, (void *)at);
    ast_add (items, item);
  }
}

#define MAX_GHOST_SCALARS 256

/* The foreach_neighbor () whose statement n is in, its place in loop's, or
 * -1: n is at the loop's own point. */
static int mover_of (const struct stencil_loop *loop, const struct ast *n,
                     const struct ast *body)
{
  for (const struct ast *p = n->parent, *prev = n; p; prev = p, p = p->parent) {
    if (p->kind == AST_STENCIL_MOVER && prev == p->last)
      for (int k = 0; k < loop->nmover; k++)
        if (loop->mover[k] == p)
          return k;
    if (p == body)
      break;
  }
  return -1;
}

/* reserved_uses () at the loop's own point: not in a foreach_neighbor ()'s
 * statement, which has x, y, z and Delta of its own. */
static void reserved_uses_outside (struct ast *n, const struct ast *body,
                                   int uses[])
{
  if (n->kind == AST_TOKEN) {
    reserved_uses (n, body, uses, 0);
    return;
  }
  for (struct ast *k = n->child; k; k = k->next)
    if (!(n->kind == AST_STENCIL_MOVER && k == n->last))
      reserved_uses_outside (k, body, uses);
}

/* s as a C string literal */
static const char *c_string (const char *s)
{
  const char *q = "\"";

  for (; *s; s++)
    q = arena_printf (*s == '"' || *s == '\\' ? "%s\\%c" : "%s%c", q, *s);
  return arena_printf ("%s\"", q);
}

/* child.x, .y, .z under n, child not declared in the body: sgn [0], [1],
 * [2], the half of the coarse cell the point is in. */
static void replace_child_sgn (struct ast *n, const struct ast *body,
                               const char *sgn)
{
  if (n->kind == AST_TOKEN)
    return;
  if (n->kind == AST_MEMBER && n->child->kind == AST_TOKEN
      && strcmp (ast_token (n->child->next)->text, ".") == 0) {
    const struct token *t = ast_token (n->child);
    const char *m = ast_token (n->last)->text;

    if (t->kind == T_IDENTIFIER && strcmp (t->text, "child") == 0
        && !(t->binding && t->binding->decl && is_under (t->binding->decl, body))
        && m[1] == '\0' && (m[0] == 'x' || m[0] == 'y' || m[0] == 'z')) {
      struct ast *e = parse_stencil_item (arena_printf ("(%s[%d]);", sgn,
                                                        m[0] - 'x'))->child;

      e->parent = e->next = NULL;
      each_leaf (e, flatten_leaf, (void *)t);
      first_token (e)->trivia = t->trivia;
      ast_replace (n, e);
      return;
    }
  }
  for (struct ast *k = n->child, *next; k; k = next) {
    next = k->next;
    replace_child_sgn (k, body, sgn);
  }
}

/* The loops' names that the cell's x, y, z and Delta are computed from. */
static int loops_give_position (const struct loop_def *leaf,
                                const struct loop_def *cell)
{
  static const char *const leaf_names[] = { "level", "i", "j", "k" };
  static const char *const cell_names[] = { "ci", "cj", "ck" };

  for (int a = 0; a <= dimension; a++)
    if (!loop_gives (leaf, leaf_names[a]))
      return 0;
  for (int a = 0; a < dimension; a++)
    if (!loop_gives (cell, cell_names[a]))
      return 0;
  return 1;
}

struct ast *stencil_loop_end (struct stencil_loop *loop,
                              const struct stencil_def *def, struct ast *args,
                              int has_ghosts, int has_reach,
                              int has_coarse_reach, struct ast *kw,
                              struct ast *body, int *rejected)
{
  int multi = 0, has_coarse = 0, has_child = 0, *level, *is_child;
  const char *coarse_sblocks, *child_sblocks, *cc, *fc, *sgn;
  int cell_uses[NRESERVED] = { 0 }, (*mover_uses)[NRESERVED], *moved;
  const struct loop_def *leaf = def ? def->leaf : NULL;
  const struct loop_def *cell = def ? def->cell : NULL;
  struct ast *ps[MAX_LOOP_PARAMS], *vs[MAX_LOOP_PARAMS];
  int nparams = 0, nargs = 0;
  struct ghost_scalar written[MAX_GHOST_SCALARS], needed[MAX_GHOST_SCALARS];
  int nwritten = 0, nneeded = 0, *wrote;
  static const char *const strides[] = { "", " * %s.sy", " * %s.sz" };
  const struct token *at = ast_token (kw);
  struct stencil_block blocks[MAX_STENCIL_BLOCKS];
  struct ast *items, *open, *close, *kids[3], *cells;
  const char *sblocks, *c, *shape;
  int nblocks = 0, serial = ++stencil_serial, inner, outer;
  const char *params[MAX_LOOP_PARAMS];
  int uses[NRESERVED] = { 0 }, position = 0;
  int errors = loop->errors;              /* this loop's: those after it */

  if (!def) {
    stencil_error (rejected, kw, "%s () is not declared: include "
                   "<block-loops.h>, which declares it with @foreach",
                   ast_token (kw)->text);
    return body;
  }
  if (!leaf || !cell)
    return body;                          /* reported where it is defined */
  if (!loop_gives (leaf, "sblocks") || !loop_gives (cell, "c")
      || !loop_gives (cell, "shape")) {
    stencil_error (rejected, kw, "%s () needs its loop over the blocks, %s, "
                   "to give sblocks, and its loop over the cells, %s, c and "
                   "shape", def->name, leaf->name, cell->name);
    return body;
  }
  {
    struct loop_def holder = { 0 };

    holder.params = def->params;
    nparams = loop_params (&holder, ps, MAX_LOOP_PARAMS);
  }

  /* Loops over two levels, by what their loops give: see translate.h */
  multi = loop_gives (cell, "cc") || loop_gives (cell, "fc");
  has_coarse = loop_gives (leaf, "coarse_sblocks") && loop_gives (cell, "cc")
               && loop_gives (cell, "sgn");
  has_child = loop_gives (leaf, "child_sblocks") && loop_gives (cell, "fc");
  is_child = allocate ((size_t)loop->nmover * sizeof (*is_child) + 1);
  for (int k = 0; k < loop->nmover; k++) {
    struct ast *mv = loop->mover[k];
    const char *name = ast_token (mv->child)->text;

    is_child[k] = strcmp (name, "foreach_child") == 0;
    if (is_child[k] && !has_child)
      stencil_error (rejected, mv, "foreach_child () in %s (), whose points "
                     "have no children: only in a loop over parents, as "
                     "foreach_restriction ()", def->name);
    if (!is_child[k] && multi)
      stencil_error (rejected, mv, "%s () in %s (), a loop over two levels: "
                     "it reads only the cells of its point", name, def->name);
  }
  level = allocate ((size_t)loop->naccess * sizeof (*level) + 1);
  for (int a = 0; a < loop->naccess; a++) {
    struct ast *access = loop->access[a];
    int m = mover_of (loop, access, body);

    level[a] = access->kind == AST_STENCIL_COARSE ? 1
               : m >= 0 && is_child[m] ? 2 : 0;
    if (level[a] == 1 && !has_coarse)
      stencil_error (rejected, access, "'%s' in %s (), whose points have no "
                     "coarse cell: only in a loop as foreach_prolongation ()",
                     access_text (access), def->name);
    else if (level[a] == 1 && stencil_writes (access))
      stencil_error (rejected, access, "'%s' writes the coarse level: "
                     "coarse () only reads", access_text (access));
    else if (level[a] != 1 && multi && stencil_has_offset (access))
      stencil_error (rejected, access, "'%s' in %s (), a loop over two "
                     "levels: no offset on the point's own level, whose cells "
                     "are being filled", access_text (access), def->name);
  }
  if (*rejected > errors)
    return body;
  for (struct ast *a = args ? args->child : NULL; a; a = a->next)
    if (!(a->kind == AST_TOKEN && strcmp (ast_token (a)->text, ",") == 0)
        && nargs < MAX_LOOP_PARAMS)
      vs[nargs++] = a;
  if (nargs != nparams) {
    stencil_error (rejected, kw, "%s () takes %d argument%s, not %d",
                   def->name, nparams, nparams == 1 ? "" : "s", nargs);
    diag (def->at->file, def->at->line, DIAG_NOTE, "'%s' is defined here",
          def->name);
    return body;
  }
  /* x, y, z and Delta are the cell's, and cannot name a scalar here */
  for (int a = 0; a < loop->naccess; a++) {
    int none[NRESERVED] = { 0 };
    struct ast *scalar = loop->access[a]->child;

    if (reserved_uses (scalar, body, none, 0))
      stencil_error (rejected, scalar, "'%s[...]': in foreach (), x, y, z and "
                     "Delta are the cell's, not a scalar",
                     expression_text (scalar));
  }
  if (*rejected > errors)
    return body;
  reserved_hides_scalar (body, body, rejected);
  if (*rejected > errors)
    return body;
  position = reserved_uses (body, body, uses, 1) > 0;
  reserved_uses_outside (body, body, cell_uses);
  mover_uses = allocate ((size_t)loop->nmover * sizeof (*mover_uses) + 1);
  for (int k = 0; k < loop->nmover; k++) {
    for (int r = 0; r < NRESERVED; r++)
      mover_uses[k][r] = 0;
    reserved_uses (loop->mover[k]->last, body, mover_uses[k], 0);
  }
  {
    int reach = 0, creach = 0;

    for (int k = 0; k < loop->nmover; k++)
      reach |= !is_child[k];
    for (int a = 0; a < loop->naccess; a++)
      if (loop->access[a]->child->next->next != loop->access[a]->last)
        *(level[a] == 1 ? &creach : &reach) = 1;
    if (creach && !has_coarse_reach) {
      stencil_error (rejected, kw, "offsets of coarse () in %s () are checked "
                     "against BAMR_COARSE_REACH, the coarse ghost layers that "
                     "are valid: include <block-loops.h>", def->name);
      return body;
    }
    if (reach && !has_reach) {
      stencil_error (rejected, kw, "offsets and foreach_neighbor () in %s () "
                     "are checked against BAMR_GHOSTS, the ghost layers: "
                     "include <block-loops.h>", def->name);
      return body;
    }
  }
  if (position && !loops_give_position (leaf, cell)) {
    stencil_error (rejected, kw, "x, y, z and Delta in foreach () need "
                   "foreach_leaf to give level, i, j, k and foreach_cell ci, "
                   "cj, ck");
    return body;
  }
  if (*rejected > errors)
    return body;

  /* Ghosts: the scalars written, whose ghosts go stale, and those read at an
   * offset, whose ghosts must be up to date first. A cell writes only its
   * own cell; a neighbour of a scalar the loop writes is read old or new by
   * the order of the cells, and its ghosts stale. */
  wrote = allocate ((size_t)loop->naccess * sizeof (*wrote) + 1);
  moved = allocate ((size_t)loop->naccess * sizeof (*moved) + 1);
  for (int a = 0; a < loop->naccess; a++) {
    struct ast *access = loop->access[a];
    int offset;

    moved[a] = mover_of (loop, access, body);
    offset = moved[a] >= 0 || stencil_has_offset (access);
    wrote[a] = stencil_writes (access);
    if (wrote[a] && moved[a] >= 0) {
      stencil_error (rejected, access, "'%s' in %s () writes a %s: a cell "
                     "writes only its own", expression_text (access),
                     ast_token (loop->mover[moved[a]]->child)->text,
                     is_child[moved[a]] ? "child" : "neighbour");
      continue;
    }
    if (multi)
      continue;                           /* in the ghost update: none tracked */
    if (wrote[a] && offset) {
      stencil_error (rejected, access, "'%s' writes a cell other than its "
                     "own: in %s (), a cell writes only s[]",
                     expression_text (access), def->name);
      continue;
    }
    if (nwritten == MAX_GHOST_SCALARS || nneeded == MAX_GHOST_SCALARS)
      continue;
    if (wrote[a])
      add_ghost_scalar (written, &nwritten, ghost_scalar_of (access, body));
    else if (offset)
      add_ghost_scalar (needed, &nneeded, ghost_scalar_of (access, body));
  }
  if (*rejected > errors)
    return body;
  for (int r = 0; r < nneeded; r++)
    for (int w = 0; w < nwritten; w++)
      if (!needed[r].all && strcmp (needed[r].key, written[w].key) == 0) {
        const struct token *t = ast_first_token (needed[r].at);
        const struct token *u = ast_first_token (written[w].at);

        diag (t->file, t->line, DIAG_WARNING, "'%s' reads a neighbour of '%s', "
              "which this foreach () writes: its ghosts are stale, and its "
              "neighbours old or new by the order of the cells; write it in a "
              "foreach () of its own first", expression_text (needed[r].at),
              needed[r].key);
        diag (u->file, u->line, DIAG_NOTE, "'%s' is written here",
              expression_text (written[w].at));
      }
  if ((nwritten || nneeded) && !has_ghosts) {
    stencil_error (rejected, kw, "foreach () keeps the ghosts of the scalars "
                   "it writes up to date: it needs ghosts_mark_dirty (), "
                   "ghosts_mark_all_dirty () and ghosts_update_dirty (): "
                   "include <block-ghosts.h>");
    return body;
  }

  /* The names the expansions below will have: foreach_cell's first, as it
   * expands first. */
  inner = loop_expansions + 1;
  outer = loop_expansions + 2;
  sblocks = STENCIL_NAME (leaf->name, "sblocks", outer);
  c = STENCIL_NAME (cell->name, "c", inner);
  shape = STENCIL_NAME (cell->name, "shape", inner);
  coarse_sblocks = STENCIL_NAME (leaf->name, "coarse_sblocks", outer);
  child_sblocks = STENCIL_NAME (leaf->name, "child_sblocks", outer);
  cc = STENCIL_NAME (cell->name, "cc", inner);
  fc = STENCIL_NAME (cell->name, "fc", inner);
  sgn = STENCIL_NAME (cell->name, "sgn", inner);

  /* Which scalars get a block, and whether it is written: before any access
   * is rewritten, which would hide the assignments around them. */
  int *which = allocate ((size_t)loop->naccess * sizeof (*which) + 1);

  for (int a = 0; a < loop->naccess; a++) {
    struct ast *scalar = loop->access[a]->child;
    const char *key = arena_printf ("%d %s", level[a], expression_text (scalar));
    int b;

    which[a] = -1;
    if (!stencil_invariant (scalar, body))
      continue;
    for (b = 0; b < nblocks && strcmp (blocks[b].key, key) != 0; b++)
      ;
    if (b == nblocks) {
      if (nblocks == MAX_STENCIL_BLOCKS)
        continue;
      blocks[nblocks].key = key;
      blocks[nblocks].scalar = ast_copy (scalar);
      blocks[nblocks].level = level[a];
      blocks[nblocks++].writes = 0;
    }
    blocks[b].writes |= stencil_writes (loop->access[a]);
    which[a] = b;
  }

  /* Each s[i, j] becomes block[c + (i) + (j) * shape.sy], where block is the
   * scalar's, or block_data (sblocks, s) for one that is not invariant. */
  for (int a = 0; a < loop->naccess; a++) {
    struct ast *access = loop->access[a], *os[4], *e;
    struct ast *scalar = access->child, *rbracket = access->last;
    int n = stencil_offsets (access->child->next->next == rbracket ? NULL
                             : access->child->next->next, os, 4);
    const char *index = level[a] == 1 ? cc : level[a] == 2 ? fc : c, *target;
    const char *levels[3] = { sblocks, coarse_sblocks, child_sblocks };
    struct token *first;

    for (int k = 0; k < (moved[a] >= 0 ? dimension : n) && k < 3; k++)
      index = arena_printf ("%s + (%s%s)%s", index,
                            k < n ? arena_printf ("__cwrap_o%d", k) : "0",
                            moved[a] >= 0
                            ? arena_printf (" + __cwrap_d%c_%d_%d", "ijk"[k],
                                            serial, moved[a]) : "",
                            arena_printf (strides[k], shape));
    /* Each offset is a constant within the ghost layers, with the
     * neighbour's in foreach_neighbor (): gcc checks it, at the access,
     * with a static assertion in the index, which costs nothing. */
    if (n > 0 && level[a] == 1) {
      const char *cond = "1";

      for (int k = 0; k < n && k < 3; k++)
        cond = arena_printf ("%s && (__cwrap_a%d) <= BAMR_COARSE_REACH "
                             "&& (__cwrap_b%d) >= -BAMR_COARSE_REACH", cond, k, k);
      index = arena_printf ("%s + (int) (0 * sizeof (struct { _Static_assert "
                            "(%s, %s); char __cwrap_c; }))", index, cond,
                            c_string (arena_printf ("%s: an offset is not a "
                                                    "constant within "
                                                    "BAMR_COARSE_REACH, the "
                                                    "coarse ghost layers that "
                                                    "are valid",
                                                    access_text (access))));
    }
    else if (n > 0) {
      const char *cond = "1";

      for (int k = 0; k < n && k < 3; k++)
        cond = moved[a] >= 0
          ? arena_printf ("%s && (__cwrap_a%d) + (__cwrap_ra%d) <= BAMR_GHOSTS "
                          "&& (__cwrap_b%d) - (__cwrap_rb%d) >= -BAMR_GHOSTS",
                          cond, k, k, k, k)
          : arena_printf ("%s && (__cwrap_a%d) <= BAMR_GHOSTS "
                          "&& (__cwrap_b%d) >= -BAMR_GHOSTS", cond, k, k);
      index = arena_printf ("%s + (int) (0 * sizeof (struct { _Static_assert "
                            "(%s, %s); char __cwrap_c; }))", index, cond,
                            c_string (arena_printf ("%s: an offset is not a "
                                                    "constant within "
                                                    "BAMR_GHOSTS, the ghost "
                                                    "layers%s",
                                                    expression_text (access),
                                                    moved[a] >= 0
                                                    ? " (adding the radius of foreach_neighbor)"
                                                    : "")));
    }
    target = which[a] >= 0 ? arena_printf ("__cwrap_block_%d_%d", serial, which[a])
                           : arena_printf ("block_data (%s, (__cwrap_s))",
                                           levels[level[a]]);
    /* in a mover, or a level up, a value, not an lvalue: gcc refuses a
     * write there too, whatever cwrap's check sees */
    e = parse_stencil_item (arena_printf (moved[a] >= 0 || level[a] == 1
                                          ? "((Real) %s[%s]);" : "%s[%s];",
                                          target, index))->child;
    e->parent = e->next = NULL;
    each_leaf (e, flatten_leaf, (void *)ast_first_token (access));
    first = first_token (e);
    first->trivia = ast_first_token (access)->trivia;
    last_token (e)->trivia = ast_token (rbracket)->trivia;
    if (which[a] < 0) {
      first_token (scalar)->trivia = "";   /* the line break went to e */
      fill_placeholder (e, "__cwrap_s", scalar);
    }
    /* copies for the assertion first: the offsets themselves move */
    for (int k = 0; k < n && k < 3; k++) {
      const char *names[4] = { "__cwrap_a%d", "__cwrap_b%d",
                               "__cwrap_ra%d", "__cwrap_rb%d" };

      for (int q = 0; q < (moved[a] >= 0 ? 4 : 2); q++) {
        struct ast *copy = ast_copy (q < 2 ? os[k]
                                     : loop->mover[moved[a]]->child->next->next);

        each_leaf (copy, flatten_leaf, (void *)ast_first_token (access));
        fill_placeholder (e, arena_printf (names[q], k), copy);
      }
    }
    for (int k = 0; k < n && k < 3; k++)
      fill_placeholder (e, arena_printf ("__cwrap_o%d", k), os[k]);
    ast_replace (access, e);
  }

  /* Each foreach_neighbor (r) statement becomes
   *   { _Static_assert (r a constant from 0 to BAMR_GHOSTS); const int R = (r);
   *     for (int dk = -R; dk <= R; dk++) ... for (int di ...)
   *       { const double x = ...; statement } }
   * its accesses already offset by di, dj, dk. */
  for (int k = 0; k < loop->nmover; k++) {
    static const char *const cijk[] = { "ci", "cj", "ck" };
    static const char *const first[] = { "ox", "g", "g" };
    struct ast *mv = loop->mover[k], *radius = mv->child->next->next;

    if (is_child[k]) {
      /* foreach_child () { s[] }: { for (int dk = 0; dk <= 1; dk++) ...
       *   { statement } }, in the order of the children, i fastest */
      struct ast *statement = mv->last, *items = ast_new (AST_BLOCK_ITEMS), *hole;
      const struct token *kt = ast_token (mv->child);
      const char *text = "{";

      for (int a = dimension - 1; a >= 0; a--)
        text = arena_printf ("%s for (int __cwrap_d%c_%d_%d = 0; "
                             "__cwrap_d%c_%d_%d <= 1; __cwrap_d%c_%d_%d++)", text,
                             "ijk"[a], serial, k, "ijk"[a], serial, k,
                             "ijk"[a], serial, k);
      add_stencil_items (items, arena_printf ("%s { __cwrap_body; } }", text),
                         NULL, 0, kt);
      hole = find_leaf (items, "__cwrap_body")->parent;
      statement->parent = statement->next = NULL;
      mark_token (first_token (statement));
      ast_replace (hole, statement);
      hole = items->child;
      hole->parent = hole->next = NULL;
      first_token (hole)->trivia = kt->trivia;
      ast_replace (mv, hole);
      continue;
    }
    struct ast *statement = mv->last, *items = ast_new (AST_BLOCK_ITEMS), *hole;
    const struct token *kt = ast_token (mv->child);
    const char *R = arena_printf ("__cwrap_r_%d_%d", serial, k);
    struct ast *rs[3] = { radius, radius, radius };
    const char *text = arena_printf ("{ _Static_assert ((__cwrap_e1) >= 0 && "
                                     "(__cwrap_e2) <= BAMR_GHOSTS, %s); "
                                     "const int %s = (__cwrap_e0);",
                                     c_string ("foreach_neighbor (): its "
                                               "radius is not a constant from "
                                               "0 to BAMR_GHOSTS, the ghost "
                                               "layers"), R);
    const char *decls = "";

    for (int a = dimension - 1; a >= 0; a--)
      text = arena_printf ("%s for (int __cwrap_d%c_%d_%d = -%s; "
                           "__cwrap_d%c_%d_%d <= %s; __cwrap_d%c_%d_%d++)", text,
                           "ijk"[a], serial, k, R, "ijk"[a], serial, k, R,
                           "ijk"[a], serial, k);
    for (int r = 0; r < NRESERVED; r++) {
      if (!mover_uses[k][r])
        continue;
      if (r == 3)
        decls = arena_printf ("%s const double Delta = __cwrap_Delta_%d;",
                              decls, serial);
      else if (r >= dimension)
        decls = arena_printf ("%s const double %s = 0;", decls,
                              stencil_reserved[r]);
      else
        decls = arena_printf ("%s const double %s = __cwrap_%s_%d + (%s - %s.%s"
                              " + __cwrap_d%c_%d_%d) * __cwrap_Delta_%d;", decls,
                              stencil_reserved[r], stencil_reserved[r], serial,
                              STENCIL_NAME (cell->name, cijk[r], inner), shape,
                              first[r], "ijk"[r], serial, k, serial);
    }
    add_stencil_items (items, arena_printf ("%s {%s __cwrap_body; } }", text,
                                            decls), rs, 3, kt);
    hole = find_leaf (items, "__cwrap_body")->parent;   /* the statement */
    statement->parent = statement->next = NULL;
    mark_token (first_token (statement));
    ast_replace (hole, statement);
    hole = items->child;
    hole->parent = hole->next = NULL;
    first_token (hole)->trivia = kt->trivia;
    ast_replace (mv, hole);
  }

  /* { const double x = ..., Delta = ...; body }: those the body uses, at
   * each cell, from what each leaf gives once */
  if (position) {
    static const char *const cijk[] = { "ci", "cj", "ck" };
    static const char *const first[] = { "ox", "g", "g" };
    struct ast *cell_items = ast_new (AST_BLOCK_ITEMS), *wrap[3];

    for (int r = 0; r < NRESERVED; r++) {
      const char *text;
      struct ast *decl;

      if (!cell_uses[r])
        continue;
      if (r == 3)
        text = arena_printf ("const double Delta = __cwrap_Delta_%d;", serial);
      else if (r >= dimension)
        text = arena_printf ("const double %s = 0;", stencil_reserved[r]);
      else
        text = arena_printf ("const double %s = __cwrap_%s_%d + (%s - %s.%s) "
                             "* __cwrap_Delta_%d;", stencil_reserved[r],
                             stencil_reserved[r], serial,
                             STENCIL_NAME (cell->name, cijk[r], inner), shape,
                             first[r], serial);
      decl = parse_stencil_item (text);
      each_leaf (decl, flatten_leaf, (void *)at);
      ast_add (cell_items, decl);
    }
    mark_token (first_token (body));
    ast_add (cell_items, body);
    wrap[0] = new_token ('{', "{", at);
    wrap[1] = cell_items;
    wrap[2] = new_token ('}', "}", last_token (body));
    body = ast_node (AST_COMPOUND_STATEMENT, wrap, 3);
  }

  /* child.x, .y, .z: the half of the coarse cell the point is in */
  if (has_coarse)
    replace_child_sgn (body, body, sgn);

  /* the cell loop, with the arguments the declaration gives it: its
   * parameters as hidden locals, bound below; the names the block loop
   * gives, its own */
  for (int k = 0; k < nparams; k++)
    params[k] = arena_printf ("__cwrap_%s_%s_%d", def->name,
                              ast_declared_name (ps[k]), serial);
  {
    struct ast *cell_args = stencil_arguments (def->cell_args, ps, params, nparams);
    struct given_names g = { leaf, outer };

    if (cell_args)
      each_leaf (cell_args, rename_given, &g);
    cells = expand_loop (cell, new_token (T_IDENTIFIER, "foreach", at),
                         cell_args, body, NULL, 1, rejected);
  }

  /* { the blocks; foreach_cell ... }, the block of foreach_leaf () */
  items = ast_new (AST_BLOCK_ITEMS);
  if (position) {
    /* Delta, and the centre of the leaf's first interior cell */
    static const char *const ijk[] = { "i", "j", "k" };
    const char *level = STENCIL_NAME (leaf->name, "level", outer);
    const char *text = arena_printf (
      "const double __cwrap_Delta_%d = tree_get_size () / ((double) "
      "tree_get_block_m () * (double) (1L << %s));", serial, level);

    for (int a = 0; a < dimension && a < 3; a++)
      if (uses[a])
        text = arena_printf ("%s const double __cwrap_%s_%d = tree_get_origin "
                             "(%d) + ((double) %s * tree_get_block_m () + 0.5) "
                             "* __cwrap_Delta_%d;", text, stencil_reserved[a],
                             serial, a, STENCIL_NAME (leaf->name, ijk[a], outer),
                             serial);
    text = arena_printf ("%s (void) __cwrap_Delta_%d;", text, serial);
    for (struct ast *fn = parse_snippet (arena_printf (
           "typedef int Real; void f (void) { %s }", text))->next,
         *item = fn->child->next->next->child->next->child, *next;
         item; item = next) {
      next = item->next;
      item->parent = item->next = NULL;
      each_leaf (item, flatten_leaf, (void *)at);
      ast_add (items, item);
    }
  }
  for (int b = 0; b < nblocks; b++) {
    struct ast *decl = parse_stencil_item (arena_printf (
      "%sReal *const __cwrap_block_%d_%d = block_data (%s, (__cwrap_s));",
      blocks[b].writes ? "" : "const ", serial, b,
      blocks[b].level == 1 ? coarse_sblocks
      : blocks[b].level == 2 ? child_sblocks : sblocks));

    each_leaf (blocks[b].scalar, flatten_leaf, (void *)at);
    fill_placeholder (decl, "__cwrap_s", blocks[b].scalar);
    each_leaf (decl, flatten_leaf, (void *)at);
    ast_add (items, decl);
  }
  ast_add (items, cells);
  open = new_token ('{', "{", at);
  close = new_token ('}', "}", last_token (cells));   /* where the body ends */
  kids[0] = open;
  kids[1] = items;
  kids[2] = close;
  struct ast *leaf_args = stencil_arguments (def->leaf_args, ps, params, nparams);

  if (!nwritten && !nneeded && !nparams)
    return stencil_mark (expand_loop (leaf, kw, leaf_args,
                                      ast_node (AST_COMPOUND_STATEMENT, kids, 3),
                                      NULL, 1, rejected));

  /* { the parameters, bound; the stale ghosts it reads, updated; the loop;
   *   those it writes, marked } */
  struct ast *leaves = expand_loop (leaf, new_token (T_IDENTIFIER, "foreach", at),
                                    leaf_args,
                                    ast_node (AST_COMPOUND_STATEMENT, kids, 3),
                                    NULL, 1, rejected);
  struct ast *around = ast_new (AST_BLOCK_ITEMS), *es[2 * MAX_GHOST_SCALARS];
  const struct token *end = last_token (leaves);
  int all = 0, nes = 0;
  const char *text = "";

  {
    struct loop_def none = { 0 };
    struct rename r = { &none, 0, NULL, 0 };

    none.serial = -1;                     /* renames nothing */
    for (int k = 0; k < nparams; k++) {
      vs[k]->parent = vs[k]->next = NULL;
      each_leaf (vs[k], flatten_leaf, (void *)at);
      ast_add (around, bind_parameter (ps[k], params[k], vs[k], &r, at));
    }
  }

  for (int r = 0; r < nneeded; r++)
    all |= needed[r].all;
  if (all)
    add_stencil_items (around, "ghosts_update_dirty (-1, 0);", NULL, 0, at);
  else if (nneeded) {
    const char *size = "0", *fill = "";

    for (int r = 0; r < nneeded; r++)
      if (needed[r].scalar) {
        es[nes] = needed[r].scalar;
        size = arena_printf ("%s + 1", size);
        fill = arena_printf ("%s __cwrap_stale_%d[__cwrap_nstale_%d++] = "
                             "(__cwrap_e%d);", fill, serial, serial, nes++);
      } else {
        es[nes] = es[nes + 1] = es[nes + 2] = es[nes + 3] = needed[r].array;
        size = arena_printf ("%s + sizeof (__cwrap_e%d) / sizeof (__cwrap_e%d)[0]",
                             size, nes, nes + 1);
        fill = arena_printf ("%s for (int __cwrap_i_%d = 0; __cwrap_i_%d < (int) "
                             "(sizeof (__cwrap_e%d) / sizeof (__cwrap_e%d)[0]); "
                             "__cwrap_i_%d++) __cwrap_stale_%d[__cwrap_nstale_%d++] "
                             "= (__cwrap_e%d)[__cwrap_i_%d];", fill, serial, serial,
                             nes + 2, nes + 3, serial, serial, serial, nes + 4,
                             serial);
        es[nes + 4] = needed[r].array;
        nes += 5;
      }
    text = arena_printf ("{ int __cwrap_stale_%d[%s], __cwrap_nstale_%d = 0;%s "
                         "ghosts_update_dirty (__cwrap_nstale_%d, "
                         "__cwrap_stale_%d); }", serial, size, serial, fill,
                         serial, serial);
    add_stencil_items (around, text, es, nes, at);
  }
  ast_add (around, leaves);

  all = 0;
  for (int w = 0; w < nwritten; w++)
    all |= written[w].all;
  nes = 0;
  text = "";
  if (all)
    text = "ghosts_mark_all_dirty ();";
  else
    for (int w = 0; w < nwritten; w++)
      if (written[w].scalar) {
        es[nes] = written[w].scalar;
        text = arena_printf ("%s ghosts_mark_dirty (__cwrap_e%d);", text, nes++);
      } else {
        es[nes] = es[nes + 1] = es[nes + 2] = written[w].array;
        text = arena_printf ("%s for (int __cwrap_i_%d = 0; __cwrap_i_%d < (int) "
                             "(sizeof (__cwrap_e%d) / sizeof (__cwrap_e%d)[0]); "
                             "__cwrap_i_%d++) ghosts_mark_dirty ((__cwrap_e%d)"
                             "[__cwrap_i_%d]);", text, serial, serial, nes,
                             nes + 1, serial, nes + 2, serial);
        nes += 3;
      }
  if (*text)
    add_stencil_items (around, text, es, nes, end);

  kids[0] = new_token ('{', "{", at);
  first_token (kids[0])->trivia = at->trivia;
  kids[1] = around;
  kids[2] = new_token ('}', "}", end);
  return stencil_mark (ast_node (AST_COMPOUND_STATEMENT, kids, 3));
}

/* Stencil loops do not nest, even through a call: a function of this file
 * that runs one, itself or through others of this file, must not be called
 * inside one. Returns the number of errors reported. */

struct stencil_fn {
  const char *name;
  struct ast *def;
  const struct token *loop;     /* its first stencil loop, or NULL */
  int runs;                     /* runs one, itself or through a call */
};

static const struct token *first_stencil_loop (struct ast *n)
{
  if (n->kind == AST_TOKEN)
    return NULL;
  if (n->kind == AST_STENCIL_LOOP)
    return ast_first_token (n);
  for (n = n->child; n; n = n->next) {
    const struct token *t = first_stencil_loop (n);

    if (t)
      return t;
  }
  return NULL;
}

static struct stencil_fn *find_fn (struct stencil_fn *fns, int n,
                                   const char *name)
{
  for (int k = 0; k < n; k++)
    if (strcmp (fns[k].name, name) == 0)
      return &fns[k];
  return NULL;
}

/* The function a call names, if it is a plain name. */
static const char *callee (const struct ast *call)
{
  const struct ast *f = call->child;

  while (f->kind == AST_PAREN_EXPRESSION)
    f = f->child->next;
  return f->kind == AST_TOKEN && ast_token (f)->kind == T_IDENTIFIER
         ? ast_token (f)->text : NULL;
}

/* Does anything under n call one that runs a stencil loop? */
static int calls_runner (struct ast *n, struct stencil_fn *fns, int nfns)
{
  struct stencil_fn *f;
  const char *name;

  if (n->kind == AST_TOKEN)
    return 0;
  if (n->kind == AST_CALL && (name = callee (n))
      && (f = find_fn (fns, nfns, name)) && f->runs)
    return 1;
  for (n = n->child; n; n = n->next)
    if (calls_runner (n, fns, nfns))
      return 1;
  return 0;
}

/* An error for each call under n, a stencil loop, of one that runs one. */
static int report_nested_calls (struct ast *n, struct stencil_fn *fns,
                                int nfns, const struct token *outer)
{
  struct stencil_fn *f;
  const char *name;
  int errors = 0;

  if (n->kind == AST_TOKEN)
    return 0;
  if (n->kind == AST_CALL && (name = callee (n))
      && (f = find_fn (fns, nfns, name)) && f->runs) {
    const struct token *t = ast_first_token (n);

    diag (t->file, t->line, DIAG_ERROR, "'%s ()' runs a stencil loop, and is "
          "called inside one: stencil loops do not nest", name);
    if (f->loop)
      diag (f->loop->file, f->loop->line, DIAG_NOTE, "'%s' runs one here",
            name);
    else {
      t = ast_first_token (f->def);
      diag (t->file, t->line, DIAG_NOTE, "'%s' calls one that runs one",
            name);
    }
    diag (outer->file, outer->line, DIAG_NOTE, "the outer loop is here");
    errors++;
  }
  for (n = n->child; n; n = n->next)
    errors += report_nested_calls (n, fns, nfns, outer);
  return errors;
}

static int find_stencil_loops (struct ast *n, struct stencil_fn *fns, int nfns)
{
  int errors = 0;

  if (n->kind == AST_TOKEN)
    return 0;
  if (n->kind == AST_STENCIL_LOOP)
    return report_nested_calls (n, fns, nfns, ast_first_token (n));
  for (n = n->child; n; n = n->next)
    errors += find_stencil_loops (n, fns, nfns);
  return errors;
}

static int check_stencil_calls (struct ast *root)
{
  struct stencil_fn *fns;
  int n = 0, cap = 0, changed = 1, errors = 0;

  for (struct ast *d = root->child; d; d = d->next)
    cap++;
  fns = allocate ((size_t)cap * sizeof (*fns) + 1);
  for (struct ast *d = root->child; d; d = d->next) {
    struct ast *fn = d->kind == AST_EXTENSION_DEFINITION ? d->child->next : d;
    const char *name;

    if (fn->kind != AST_FUNCTION_DEFINITION
        || !(name = ast_declared_name (fn->child->next)))
      continue;
    fns[n].name = name;
    fns[n].def = fn;
    fns[n].loop = first_stencil_loop (fn);
    fns[n].runs = fns[n].loop != NULL;
    n++;
  }
  while (changed) {
    changed = 0;
    for (int k = 0; k < n; k++)
      if (!fns[k].runs && calls_runner (fns[k].def, fns, n))
        fns[k].runs = changed = 1;
  }
  for (int k = 0; k < n; k++)
    errors += find_stencil_loops (fns[k].def, fns, n);
  return errors;
}

/* Ghost operators: @prolongation, @restriction, @boundary */

struct ghost_op {
  const char *kind;         /* prolongation, restriction, boundary */
  const char *name;
  int id;
  int nuses;                /* the scalars it depends on: its parameters
                               after the first */
  const struct token *at;
  struct ghost_op *next;
};
static struct ghost_op *ghost_ops;     /* this file's, newest first */

/* @prolongation -> prolongation */
static const char *operator_kind (const struct ast *at)
{
  return ast_token (at)->text + 1;
}

static int operator_count (const char *kind)
{
  int n = 0;

  for (const struct ghost_op *o = ghost_ops; o; o = o->next)
    n += strcmp (o->kind, kind) == 0;
  return n;
}

/* The trivia of every token under n, to keep the lines of what is dropped. */
static const char *all_trivia (const struct ast *n)
{
  MEMSTREAM (t);

  if (!t)
    abort ();
  subtree_trivia (n, t);
  return MEMSTREAM_TEXT (t);
}

struct ast *operator_define (struct ast *at, struct ast *name, struct ast *open,
                             struct ast *params, struct ast *close,
                             struct ast *deps, struct ast *body, int *rejected)
{
  const char *kind = operator_kind (at), *nm = ast_token (name)->text;
  const struct token *t = ast_token (at);
  struct ghost_op *o;
  struct ast *fn, *decl, *suffix, *kids[2];

  for (o = ghost_ops; o; o = o->next)
    if (strcmp (o->kind, kind) == 0 && strcmp (o->name, nm) == 0) {
      stencil_error (rejected, name, "@%s %s is defined twice", kind, nm);
      diag (o->at->file, o->at->line, DIAG_NOTE, "'%s' is defined here", nm);
      return loop_stub (at);
    }
  o = allocate (sizeof (*o));
  o->kind = kind;
  o->name = nm;
  o->id = operator_count (kind) + 1;
  /* the scalars it depends on, (scalar zb, ...), after its parameters */
  o->nuses = 0;
  if (deps) {
    struct ast *list = deps->child->next;

    for (struct ast *p = list->child; p; p = p->next) {
      const struct ast *d;

      if (p->kind != AST_PARAMETER)
        continue;
      d = p->child->next;
      o->nuses++;
      if (!scalar_type (p->child) || strcmp (scalar_type (p->child), "scalar") != 0
          || !d || d->kind != AST_DECL_NAME)
        stencil_error (rejected, p, "@%s %s: what it depends on is scalars, "
                       "each declared 'scalar name'", kind, nm);
    }
  }
  o->at = ast_token (name);
  o->next = ghost_ops;
  ghost_ops = o;

  /* static inline void __cwrap_KIND_NAME (params) body  const int NAME = id; */
  fn = parse_snippet (arena_printf ("static inline void __cwrap_%s_%s (void) {} "
                                    "const int %s = %d;", kind, nm, nm, o->id));
  decl = fn->next;
  fn->next = decl->next = NULL;
  each_leaf (fn, flatten_leaf, (void *)t);
  each_leaf (decl, flatten_leaf, (void *)t);
  first_token (fn)->trivia = arena_printf ("%s%s", t->trivia,
                                           all_trivia (name));
  /* the parameters, with the trivia of the parentheses, and the body */
  suffix = find_kind (fn->child->next, AST_FUNCTION_SUFFIX);
  first_token (suffix)->trivia = ast_token (open)->trivia;
  last_token (suffix)->trivia = ast_token (close)->trivia;
  if (params)
    ast_replace (suffix->child->next, params);
  if (deps) {
    /* (params) (deps) becomes (params, deps), the trivia of the
     * parentheses dropped kept on the comma */
    struct ast *list = deps->child->next, *comma;
    const char *trivia = arena_printf ("%s%s", ast_token (close)->trivia,
                                       ast_token (deps->child)->trivia);

    if (params) {
      comma = new_token (',', ",", t);
      first_token (comma)->trivia = trivia;
      ast_add (params, comma);
      for (struct ast *p = list->child, *next; p; p = next) {
        next = p->next;
        p->parent = p->next = NULL;
        ast_add (params, p);
      }
    } else {
      list->parent = list->next = NULL;
      first_token (list)->trivia = arena_printf ("%s%s", trivia,
                                                 first_token (list)->trivia);
      ast_replace (suffix->child->next, list);
    }
    last_token (suffix)->trivia = ast_token (deps->last)->trivia;
  }
  ast_replace (fn->last, body);
  kids[0] = fn;
  kids[1] = decl;
  return ast_node (AST_OPERATOR_DEFINITION, kids, 2);
}

struct ast *operator_declare (struct ast *at, struct ast *name, struct ast *semi)
{
  const struct token *t = ast_token (at);
  struct ast *decl = parse_snippet (arena_printf ("extern const int %s;",
                                                  ast_token (name)->text));

  decl->next = NULL;
  each_leaf (decl, flatten_leaf, (void *)t);
  first_token (decl)->trivia = arena_printf ("%s%s", t->trivia,
                                             all_trivia (name));
  last_token (decl)->trivia = ast_token (semi)->trivia;
  return decl;
}

struct ast *operator_dispatch (struct ast *at, struct ast *open, struct ast *id,
                               struct ast *close, struct ast *open2,
                               struct ast *args, struct ast *close2,
                               struct ast *semi, int *rejected)
{
  const char *kind = operator_kind (at), *text;
  const struct token *t = ast_token (at);
  struct ast *vs[MAX_LOOP_PARAMS], **es, *items = ast_new (AST_BLOCK_ITEMS), *sw;
  int nv = 0, n = operator_count (kind), nes = 1;
  struct ghost_op **by_id;

  struct ast *uses = NULL;
  int need_uses = 0;

  (void) open, (void) close, (void) open2, (void) close2;
  if (id->kind == AST_COMMA && id->child->kind != AST_COMMA) {
    uses = id->last;
    id = id->child;
    uses->parent = uses->next = id->parent = id->next = NULL;
  }
  for (struct ghost_op *o = ghost_ops; o; o = o->next)
    need_uses |= strcmp (o->kind, kind) == 0 && o->nuses > 0;
  if (need_uses && !uses) {
    stencil_error (rejected, at, "@%s (id) (...): some @%s depends on scalars: "
                   "pass them, @%s (id, uses) (...)", kind, kind, kind);
    return ast_new (AST_EXPRESSION_STATEMENT);
  }
  if (!n) {
    stencil_error (rejected, at, "@%s (...) (...): no @%s is defined before "
                   "it in this file", kind, kind);
    return ast_new (AST_EXPRESSION_STATEMENT);
  }
  for (struct ast *a = args ? args->child : NULL; a; a = a->next)
    if (!(a->kind == AST_TOKEN && strcmp (ast_token (a)->text, ",") == 0)
        && nv < MAX_LOOP_PARAMS)
      vs[nv++] = a;
  by_id = allocate ((size_t)(n + 1) * sizeof (*by_id));
  for (struct ghost_op *o = ghost_ops; o; o = o->next)
    if (strcmp (o->kind, kind) == 0)
      by_id[o->id] = o;
  {
    int most = 0;                         /* the most scalars one depends on */

    for (int k = 1; k <= n; k++)
      most = by_id[k]->nuses > most ? by_id[k]->nuses : most;
    es = allocate ((size_t)(1 + n * (nv + most)) * sizeof (*es) + 1);
  }
  es[0] = id;
  text = "switch (__cwrap_e0) {";
  for (int k = 1; k <= n; k++) {
    const char *call = "";

    for (int j = 0; j < nv; j++) {
      es[nes] = vs[j];
      call = arena_printf ("%s%s(__cwrap_e%d)", call, j ? ", " : "", nes++);
    }
    /* uses [0], uses [1], ...: a copy of uses for each, as each placeholder
     * is filled once */
    for (int u = 0; u < by_id[k]->nuses; u++) {
      es[nes] = uses;
      call = arena_printf ("%s%s(__cwrap_e%d)[%d]", call, nv || u ? ", " : "",
                           nes++, u);
    }
    text = arena_printf ("%s case %d: __cwrap_%s_%s (%s); break;", text, k,
                         kind, by_id[k]->name, call);
  }
  text = arena_printf ("%s default: __builtin_trap (); }", text);
  add_stencil_items (items, text, es, nes, t);
  sw = items->child;
  sw->parent = sw->next = NULL;
  first_token (sw)->trivia = t->trivia;
  last_token (sw)->trivia = arena_printf ("%s%s", all_trivia (id),
                                          ast_token (semi)->trivia);
  if (args)
    last_token (sw)->trivia = arena_printf ("%s%s", last_token (sw)->trivia,
                                            all_trivia (args));
  return sw;
}

struct ast *operator_uses (struct ast *at, struct ast *word, struct ast *open,
                           struct ast *id, struct ast *close, int *rejected)
{
  const char *kind = operator_kind (at), *text = "";
  const struct token *t = ast_token (at);
  struct ast *items = ast_new (AST_BLOCK_ITEMS), *e;

  (void) open, (void) close;
  if (strcmp (ast_token (word)->text, "uses") != 0) {
    stencil_error (rejected, word, "@%s %s (...) in an expression: only @%s "
                   "uses (id)", kind, ast_token (word)->text, kind);
    return ast_new (AST_EXPRESSION_STATEMENT);
  }
  for (struct ghost_op *o = ghost_ops; o; o = o->next)
    if (strcmp (o->kind, kind) == 0)
      text = arena_printf ("%s__cwrap_u == %d ? %d : ", text, o->id, o->nuses);
  add_stencil_items (items, arena_printf ("({ int __cwrap_u = (__cwrap_e0); "
                                          "(void) __cwrap_u; %s-1; });", text),
                     &id, 1, t);
  e = items->child->child;
  e->parent = e->next = NULL;
  first_token (e)->trivia = t->trivia;
  return e;
}

/* ghost_ops (s, prolongation (op, zb, ...), restriction (...), boundary (...)) */
struct ast *ghost_ops_call (struct ast *fn, struct ast *args, int *rejected)
{
  static const char *const kinds[3] = { "prolongation", "restriction", "boundary" };
  static const char *const consts[3] = { "GHOST_PROLONGATION", "GHOST_RESTRICTION",
                                         "GHOST_BOUNDARY" };
  const struct token *t;
  struct ast *target, *es[4 * MAX_LOOP_PARAMS], *items = ast_new (AST_BLOCK_ITEMS), *e;
  const char *text = "(";
  int nes = 0, errors = *rejected;

  if (!dimension || fn->kind != AST_TOKEN)
    return NULL;
  t = ast_token (fn);
  if (t->kind != T_IDENTIFIER || strcmp (t->text, "ghost_ops") != 0 || t->binding)
    return NULL;
  target = args->child;
  if (!is_scalar_expression (target))
    stencil_error (rejected, target, "ghost_ops (%s, ...): '%s' is not a scalar",
                   expression_text (target), expression_text (target));
  for (struct ast *item = target->next; item; item = item->next) {
    struct ast *vs[MAX_LOOP_PARAMS];
    const char *deps = "";
    int kind = -1, nv = 0;

    if (item->kind == AST_TOKEN)              /* , */
      continue;
    if (item->kind == AST_CALL && item->child->kind == AST_TOKEN)
      for (int k = 0; k < 3; k++)
        if (strcmp (ast_token (item->child)->text, kinds[k]) == 0)
          kind = k;
    if (kind < 0 || item->child->next->next->kind != AST_ARGUMENT_LIST) {
      stencil_error (rejected, item, "ghost_ops (%s, ...): '%s' is not "
                     "prolongation (op, ...), restriction (op, ...) or "
                     "boundary (op, ...)", expression_text (target),
                     expression_text (item));
      continue;
    }
    for (struct ast *a = item->child->next->next->child; a; a = a->next)
      if (!(a->kind == AST_TOKEN && strcmp (ast_token (a)->text, ",") == 0)
          && nv < MAX_LOOP_PARAMS)
        vs[nv++] = a;
    /* the scalars it depends on: global, never freed */
    for (int k = 1; k < nv; k++) {
      const struct ast *specs;
      int member;

      if (!is_scalar_expression (vs[k]) || !global_scalar_decl (vs[k], &specs, &member))
        stencil_error (rejected, vs[k], "%s (%s, ...): '%s', a scalar it "
                       "depends on, must be global (at file scope, static or "
                       "extern), as a global scalar is never freed",
                       kinds[kind], expression_text (vs[0]),
                       expression_text (vs[k]));
      es[nes] = vs[k];
      deps = arena_printf ("%s%s(__cwrap_e%d)", deps, k > 1 ? ", " : "", nes++);
    }
    es[nes] = target;
    es[nes + 1] = vs[0];
    text = arena_printf ("%sghost_ops_set ((__cwrap_e%d), %s, (__cwrap_e%d), %s, %d), ",
                         text, nes, consts[kind], nes + 1,
                         nv > 1 ? arena_printf ("(const scalar []) { %s }", deps)
                                : "(const scalar *) 0", nv - 1);
    nes += 2;
  }
  if (*rejected > errors)
    return NULL;
  add_stencil_items (items, arena_printf ("%s(void) 0);", text), es, nes, t);
  e = items->child->child;
  e->parent = e->next = NULL;
  first_token (e)->trivia = t->trivia;
  return e;
}

/* The translation itself: rewrite the tree here. */
static int transform (struct ast *root)
{
  int errors = check_stencil_calls (root);

  register_scalars (root);
  define_caches (root);
  return errors;
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
  ghost_ops = NULL;
  root = parse_file (in, &rejected);
  if (rejected) {
    fclose (in);
    free_arena ();
    return TRANSLATE_FAILED;
  }
  parsed = root != NULL;
  if (root && transform (root)) {
    fclose (in);
    free_arena ();
    return TRANSLATE_FAILED;
  }
  if (root) {
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
