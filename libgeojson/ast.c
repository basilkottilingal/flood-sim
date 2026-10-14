#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "ast.h"

/*
.. Abstract Syntax Tree for (Geo)JSON
*/

/*fixme : after testing, remove [almost all of the] asserts */


/*
.. allocation, reuse and free
*/

typedef struct Allocator
{
  char * mem;
  struct Allocator * next;
  size_t left;
} Allocator;
struct Free
{
  struct Free * next;
};

static struct Free * ast_fhead = NULL, * astl_fhead = NULL;
static Allocator dummy = (Allocator) {0};
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
  ast_fhead = astl_fhead = NULL;
}

static void * allocate (size_t s)
{
  size_t sptr = sizeof (void *), mask = sptr - 1;
  s = (s + mask) & ~mask;
  if (allocator->left < s)
  {
    #ifndef _AST_PAGE_SIZE_
    #define _AST_PAGE_SIZE_ (1<<16)
    #endif
    char * m = malloc (_AST_PAGE_SIZE_);
    assert ("ast.h : malloc ()" && m);
    memset (m, 0, _AST_PAGE_SIZE_);
    Allocator * old = allocator;
    allocator = (Allocator *) m;
    
    size_t salloc = (sizeof (Allocator) + mask) & ~mask;
    *allocator = (Allocator) 
      {
        .mem  = m + salloc,
        .next = old,
        .left = (size_t) _AST_PAGE_SIZE_ - salloc
      };
  }
  char * m = allocator->mem;
  allocator->left -= s;
  allocator->mem  += s;
  return m;
}

static Ast * ast ()
{
  if (ast_fhead)
  {
    Ast * a = (Ast *) ast_fhead;
    ast_fhead = ast_fhead->next;
    *a = (Ast) {0};
    return a;
  }
  return (Ast *) allocate (sizeof (Ast));
}

static void ast_free (Ast * a)
{
  struct Free * f = (struct Free *) a;
  f->next = ast_fhead;
  ast_fhead = f;
}

static AstLeaf * astl ()
{
  if (astl_fhead)
  {
    AstLeaf * l = (AstLeaf *) astl_fhead;
    astl_fhead = astl_fhead->next;
    *l = (AstLeaf) {0};
    return l;
  }
  return (AstLeaf *) allocate (sizeof (Ast));
}

static void astl_free (AstLeaf * l)
{
  struct Free * f = (struct Free *) l;
  f->next = astl_fhead;
  astl_fhead = f;
}

Ast * ast_json (Ast * val)
{
  Ast * json = ast ();
  json->child = val, val->parent = json;
  json->type = ast_token (jjson) ;
  return json;
}

Ast * ast_val (Ast * a)
{
  Ast * val = ast ();
  assert (a);
  a->parent = val;
  val->child = a;
  val->type = ast_token (jvalue);
  return val;
}

Ast * ast_obj (Ast * members)
{
  Ast * parent = members ? members->parent : ast ();
  assert (parent);
  parent->type = ast_token (jobject);
  return parent;
}

Ast * ast_append (Ast * prev, Ast * a)
{
  Ast * parent =  prev ? prev->parent : ast ();
  assert (parent && a);
  a->parent = parent;
  if (prev)
  {
    assert (!prev->next);
    prev->next = a;
  }
  else
  {
    assert (!parent->child);
    parent->child = a;
  }
  return a;
}

Ast * ast_pair (Ast * str, Ast * val)
{
  Ast * pair = ast ();
  pair->type = ast_token (jpair);
  pair->child = str, str->next = val;
  assert (str->type == ast_token (jstring) && val->type == ast_token (jvalue));
  /* maybe you can reuse 'str' as 'pair' */
  return pair; 
}

Ast * ast_array (Ast * elements)
{
  Ast * parent = elements ? elements->parent : ast ();
  assert (parent);
  parent->type = ast_token (jarray);
  return parent;
}

Ast * ast_leaf (char * yytext, int yyleng, int tkn)
{
  AstLeaf * l = astl ();
  Ast * a = (Ast *) l;
  l->start = yytext, l->end = yytext + yyleng;
  assert (tkn);
  a->type = tkn;
  return a;
}
