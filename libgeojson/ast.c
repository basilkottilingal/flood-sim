#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "ast.h"
#include "allocator.h"

/*
.. Abstract Syntax Tree for JSON
*/

/*fixme : after testing, remove [almost all of the] asserts */

#define MAGIC 13947582016482937511ULL

Allocator * ast_allocator (Ast * root)
{
  AstRoot * r = (AstRoot *) root;
  assert (r->magic == MAGIC && r->alloc);
  return r->alloc;
}

void ast_free (Ast * root)
{
  if (root)
    allocator_free (ast_allocator (root));
}

static char * alloc_str (Allocator * alloc, const char * yytext, int yyleng, int is_string)
{
  if (is_string)
  {
    /* remove the quotes */
    yyleng -= 2;
    yytext += 1; 
  }
  char * s = allocate (alloc, yyleng + 1);
  memcpy (s, yytext, yyleng);
  s [yyleng] = '\0';
  return s;
}

static inline Ast * ast (Allocator * alloc)
{
  return allocate (alloc, sizeof (Ast));
}

static inline Ast * astr (Allocator * alloc)
{
  AstRoot * r = allocate (alloc, sizeof (AstRoot));
  r->alloc = alloc;
  r->magic = MAGIC;
  return (Ast *) r;
}

static inline AstLeaf * astl (Allocator * alloc)
{
  return allocate (alloc, sizeof (AstLeaf));
}

Ast * ast_json (Allocator * alloc, Ast * val)
{
  Ast * json = astr (alloc);
  json->child = val, val->parent = json;
  json->type = ast_token (jjson) ;
  return json;
}

Ast * ast_val (Allocator * alloc, Ast * a)
{
  Ast * val = ast (alloc);
  assert (a);
  a->parent = val;
  val->child = a;
  val->type = ast_token (jvalue);
  return val;
}

Ast * ast_obj (Allocator * alloc, Ast * members)
{
  Ast * parent = members ? members->parent : ast (alloc);
  assert (parent);
  parent->type = ast_token (jobject);
  return parent;
}

Ast * ast_append (Allocator * alloc, Ast * prev, Ast * a)
{
  Ast * parent =  prev ? prev->parent : ast (alloc);
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

Ast * ast_pair (Allocator * alloc, Ast * str, Ast * val)
{
  Ast * pair = ast (alloc);
  pair->type = ast_token (jpair);
  pair->child = str, str->next = val;
  assert (str->type == ast_token (jstring) && val->type == ast_token (jvalue));
  str->parent = val->parent = pair;
  /* maybe you can reuse 'str' as 'pair' */
  return pair; 
}

Ast * ast_array (Allocator * alloc, Ast * elements)
{
  Ast * parent = elements ? elements->parent : ast (alloc);
  assert (parent);
  parent->type = ast_token (jarray);
  return parent;
}

Ast * ast_leaf (Allocator * alloc, char * yytext, int yyleng, int tkn)
{
  AstLeaf * l = astl (alloc);
  Ast * a = (Ast *) l;
  l->token = alloc_str (alloc, yytext, yyleng, tkn == ast_token (jstring));
  assert (ast_is_leaf (tkn));
  a->type = tkn;
  return a;
}

void ast_print (Ast * a)
{
  if (!a || ((AstRoot *)a)->magic != MAGIC)
  {
    fprintf (stderr, "\nast_print () : missing AST root");
    return;
  }

  /* warning : not formatted properly */
  int depth = 0, depthmax = 15;
  do {
    while (a->child)
    {
      if (a->type == ast_token (jpair))
      {
        printf ("\n"); for (int i=0; i<depth/3; ++i) printf ("  ");
        printf ("\"%s\" :", ((AstLeaf *)a->child)->token);
        depth++, a = a->child->next;
        continue;
      }
      //for (int i=0; i<depth; ++i) printf ("  ");
      printf ("%s", a->type == ast_token (jobject) ? " {" :
                    a->type == ast_token (jarray)  ? " [" : "");

      if (depth == depthmax)
      {
        printf (" ... ");
        break;
      }

      depth++, a = a->child;

      if (ast_is_leaf (a->type))
      {
        if (a->parent->type != ast_token (jvalue))
        {
          printf ("\n");
          for (int i=0; i<depth/3; ++i)
            printf ("  ");
        }
        if (a->type == ast_token (jstring))
            printf ("\"%s\"", ((AstLeaf *)a)->token);
        else
          printf ("%s", ((AstLeaf *)a)->token);
      }
      else if (!a->child)
      {
        printf ("%s", a->type == ast_token (jobject) ? " {}" :
          a->type == ast_token (jarray)  ? " []" : "");
      }
    }
    /* pop from stack */
    while (depth)
    {
      if (a->next) {
        printf (",");
        a = a->next; break;
      }
      depth--, a = a->parent;
      if (a->type == ast_token (jobject) || a->type == ast_token (jarray))
      {
        printf ("\n"); for (int i=0; i<depth/3; ++i) printf ("  ");
        printf ("%s", a->type == ast_token (jobject) ? "}" :
                    a->type == ast_token (jarray)  ? "]" : "");
      }
    }
  } while (depth);
  
}
