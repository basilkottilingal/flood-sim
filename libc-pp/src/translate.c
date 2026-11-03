#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "allocator.h"
#include "ast.h"
#include "parser.tab.h"
#include "scope.h"
#include "translate.h"

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
  fprintf (stderr, "%s:%d: warning: '%s' changes scalar '%s', whose value "
           "cwrap gives it\n", at->file, at->line, at->text, name->text);
  at = ast_first_token (v->decl);
  fprintf (stderr, "%s:%d: note: '%s' is declared here\n", at->file, at->line,
           name->text);
}

static void warn_scalar_initializer (const struct ast *init)
{
  const struct token *t = ast_first_token (init);
  const char *name = ast_declared_name (init->child);

  fprintf (stderr, "%s:%d: warning: scalar '%s' has a value of its own, so "
           "cwrap does not number it\n", t->file, t->line, name);
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
    fprintf (stderr, "cwrap: fmemopen: %s\n", strerror (errno));
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

/* Each `scalar a;` at file scope, headers included, gets a slot in the section
 * scalars, on the same line so that no line moves:
 *
 *   scalar a; static scalar *const __cwrap_scalar_slot_0 ... = &a;
 *
 * The linker joins the slots of every file into one array, which the function
 * from add_scalar_constructor () numbers at startup. extern and typedef
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
      char text[512];

      if (d->kind == AST_INIT_DECLARATOR && plain_name (d->child))
        warn_scalar_initializer (d);
      if (!plain_name (d))
        continue;
      snprintf (text, sizeof (text), "static scalar *const __cwrap_scalar_slot_%d "
                "__attribute__ ((section (\"scalars\"), used)) = &%s;",
                slot++, ast_declared_name (d));
      reg = parse_snippet (text);
      first_token (reg)->trivia = " ";
      insert_all_before (next, reg);
    }
  }
}

/* The definition of main at file scope, or NULL. */
static struct ast *find_main (struct ast *root)
{
  for (struct ast *n = root->child; n; n = n->next) {
    struct ast *def = n->kind == AST_EXTENSION_DEFINITION ? n->child->next : n;
    const char *name;

    if (def->kind != AST_FUNCTION_DEFINITION)
      continue;
    name = ast_declared_name (def->child->next);
    if (name && strcmp (name, "main") == 0)
      return n;
  }
  return NULL;
}

/* Is n `typedef ... scalar;`? */
static int declares_scalar_type (const struct ast *n)
{
  const struct ast *list;

  if (n->kind == AST_EXTENSION_DECLARATION)
    n = n->child->next;
  if (n->kind != AST_DECLARATION || !has_token (n->child, T_TYPEDEF, NULL))
    return 0;
  list = n->child->next;
  if (list->kind != AST_INIT_DECLARATOR_LIST)
    return 0;
  for (const struct ast *d = list->child; d; d = d->next) {
    const char *name = ast_declared_name (d);

    if (name && strcmp (name, "scalar") == 0)
      return 1;
  }
  return 0;
}

/* Numbers the slots of register_scalars (), from every file of the program, in
 * the order the linker puts them. __start_scalars and __stop_scalars come from
 * the linker, and are weak for a program with no slot at all. n_global_scalars
 * is where scalars numbered at run time may start. */
static const char scalar_constructor[] =
  "extern scalar *const __start_scalars[] __attribute__ ((weak)); "
  "extern scalar *const __stop_scalars[] __attribute__ ((weak)); "
  "int n_global_scalars; "
  "__attribute__ ((constructor (101))) static void number_global_scalars (void) "
  "{ for (scalar *const *p = __start_scalars; p < __stop_scalars; p++) "
  "**p = p - __start_scalars; "
  "n_global_scalars = __stop_scalars - __start_scalars; }";

/* In the file that defines main, so in one file of the program, adds
 * scalar_constructor on a line of its own just before main. It needs the type
 * scalar, so it is added only if a typedef of scalar comes before main. main
 * gets a linemarker for its own line, so that gcc's line numbers do not move. */
static void add_scalar_constructor (struct ast *root)
{
  struct ast *main_def = find_main (root), *n, *decls;
  struct token *ours, *theirs;
  char *marker;
  int len;

  if (!main_def)
    return;
  for (n = root->child; n != main_def && !declares_scalar_type (n); n = n->next)
    ;
  if (n == main_def)
    return;
  decls = parse_snippet (scalar_constructor);
  ours = first_token (decls);
  theirs = first_token (main_def);
  ours->trivia = theirs->trivia;
  len = snprintf (NULL, 0, "\n# %d \"%s\"\n", theirs->line, theirs->file);
  marker = allocate (len + 1);
  snprintf (marker, len + 1, "\n# %d \"%s\"\n", theirs->line, theirs->file);
  theirs->trivia = marker;
  insert_all_before (main_def, decls);
}

/* The translation itself: rewrite the tree here. */
static void transform (struct ast *root)
{
  register_scalars (root);
  add_scalar_constructor (root);
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
    fprintf (stderr, "cwrap: %s: %s\n", in_path, strerror (errno));
    return TRANSLATE_FAILED;
  }
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
    fprintf (stderr, "cwrap: %s: %s\n", out_path, strerror (errno));
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
    fprintf (stderr, "cwrap: %s: write failed\n", out_path);
    return TRANSLATE_FAILED;
  }
  return parsed ? TRANSLATE_OK : TRANSLATE_SKIPPED;
}
