#include "ast.h"

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

typedef struct Allocator
{
  char * mem;
  struct Allocator * next;
  size_t left;
} Allocator;

static Allocator dummy = (Allocator) {NULL, NULL, 0};
static Allocator * allocator = & dummy;

void ast_free_all ()
{
  while (allocator->next)
  {
    void * m = allocator;
    allocator = allocator->next;
    free (m);
  }
  assert (allocator == &dummy);
}

static void * allocate (size_t s)
{
  size_t sptr = sizeof (void *), mask = sptr - 1;
  s = (s + mask) & ~mask;
  if (allocator->left < s)
  {
    char * m = malloc (1<<16);
    assert ("ast.h : malloc ()" && m);
    Allocator * old = allocator;
    allocator = (Allocator *) m;
    size_t salloc = (sizeof (Allocator) + mask) & ~mask;
    *allocator = (Allocator) 
      {
        .mem  = m + salloc,
        .next = old,
        .left = (size_t) (1 << 16) - salloc
      };
  }
  char * m = allocator->mem;
  allocator->left -= s;
  allocator->mem  += s;
  return m;
}

Ast * ast1 (Ast * child)
{
  Ast * p = (Ast *) allocate (sizeof (Ast));
  Ast ** children = (Ast **) allocate (2 * sizeof (Ast *));
  p->children = children;
  p->parent = NULL;
  children [0] = child, children [1] = NULL;
  child->parent = p;
  return p;
}

Ast * ast2 (Ast * c0, Ast * c1)
{
  Ast * p = (Ast *) allocate (sizeof (Ast));
  Ast ** children = (Ast **) allocate (3 * sizeof (Ast *));
  p->children = children;
  p->parent = NULL;
  children [0] = c0, children [1] = c1, children [2] = NULL;
  c0->parent = c1->parent = p;
  return p;
}

Ast * ast3 (Ast * c0, Ast * c1, Ast * c2)
{
  Ast * p = (Ast *) allocate (sizeof (Ast));
  Ast ** children = (Ast **) allocate (4 * sizeof (Ast *));
  p->children = children;
  p->parent = NULL;
  children [0] = c0, children [1] = c1, children [2] = c2, children [3] = NULL;
  c0->parent = c1->parent = c2->parent = p;
  return p;
}

Ast * ast_leaf (char * yytext)
{
  AstLeaf * l = (AstLeaf *) allocate (sizeof (AstLeaf));
  l->start = yytext, l->end = yytext + strlen (yytext);
  Ast * p = (Ast *) l;
  p->children = NULL;
  return p;
}

