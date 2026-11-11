#include <stdlib.h>

#include "allocator.h"
#include "ast.h"

static const char *const kind_names[] = {
#define X(name) #name,
  AST_KINDS (X)
#undef X
};

/* allocate () [in allocator.h] returns zeroed memory, so only non-zero
 * fields are set below.
 */

struct ast *ast_new (enum ast_kind kind)
{
  struct ast *n = allocate (sizeof (*n));

  n->kind = kind;
  return n;
}

struct ast *ast_leaf (const struct token *t)
{
  struct ast_leaf *l = allocate (sizeof (*l));

  l->a.kind = AST_TOKEN;
  l->tok = *t;
  return &l->a;
}

struct ast *ast_node (enum ast_kind kind, struct ast **kids, size_t n)
{
  struct ast *node = ast_new (kind);

  for (size_t i = 0; i < n; i++)
    ast_add (node, kids[i]);
  return node;
}

/* create a copy of AST subtree 'n'. utilises recursive calls */
struct ast *ast_copy (const struct ast *n)
{
  struct ast *c;

  if (n->kind == AST_TOKEN)
    return ast_leaf (ast_token (n));
  c = ast_new (n->kind);
  for (n = n->child; n; n = n->next)
    ast_add (c, ast_copy (n));
  return c;
}

/* append a child at the last*/
struct ast *ast_add (struct ast *parent, struct ast *kid)
{
  if (!kid)
    return parent;
  kid->parent = parent;
  if (parent->last)
    parent->last->next = kid;
  else
    parent->child = kid;
  parent->last = kid;
  return parent;
}

struct ast *ast_replace (struct ast *old, struct ast *n)
{
  struct ast *parent = old->parent, **link = &parent->child;

  while (*link != old)
    link = &(*link)->next;
  *link = n;
  n->parent = parent;
  n->next = old->next;
  if (parent->last == old)
    parent->last = n;
  old->parent = old->next = NULL;
  return old;
}

struct ast *ast_insert_before (struct ast *pos, struct ast *n)
{
  struct ast *parent = pos->parent, **link = &parent->child;

  while (*link != pos)
    link = &(*link)->next;
  *link = n;
  n->parent = parent;
  n->next = pos;
  return n;
}

const struct token *ast_token (const struct ast *n)
{
  return &((const struct ast_leaf *)n)->tok;
}

const struct token *ast_first_token (const struct ast *n)
{
  while (n->kind != AST_TOKEN)
    n = n->child;
  return ast_token (n);
}

/* The name always comes before a declarator's suffixes, so the first DECL_NAME
 * in tree order is this declarator's own, not a parameter's. */
const char *ast_declared_name (const struct ast *n)
{
  if (n->kind == AST_DECL_NAME)
    return ast_token (n->child)->text;
  for (n = n->child; n; n = n->next) {
    const char *name = ast_declared_name (n);

    if (name)
      return name;
  }
  return NULL;
}

/* emit the subtree 'n' recursively */
void ast_emit (const struct ast *n, FILE *out)
{
  if (n->kind == AST_TOKEN) {
    const struct token *t = ast_token (n);

    fputs (t->trivia, out);
    fputs (t->text, out);
    return;
  }
  for (n = n->child; n; n = n->next)
    ast_emit (n, out);
}

static void dump (const struct ast *n, int depth, FILE *out)
{
  fprintf (out, "%*s", depth * 2, "");
  if (n->kind == AST_TOKEN) {
    const struct token *t = ast_token (n);

    fprintf (out, "'%s' %s:%d\n", t->text, t->file, t->line);
    return;
  }
  fprintf (out, "%s\n", kind_names[n->kind]);
  for (n = n->child; n; n = n->next)
    dump (n, depth + 1, out);
}

void ast_dump (const struct ast *n, FILE *out)
{
  dump (n, 0, out);
}
