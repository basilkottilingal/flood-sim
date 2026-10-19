#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <assert.h>

#include "ast.h"
#include "json.tab.h"
#include "allocator.h"
#include "geojson.h"

static uint64_t hash (const char *s)
{
  /* FNV-1a, 64-bit */
  uint64_t h = 1469598103934665603ULL;
  for (const unsigned char *p = (const unsigned char *)s; *p; p++)
    h = (h ^ *p) * 1099511628211ULL;
  return h;
}

static size_t bin_size (size_t n)
{
  size_t p = 16;
  while (p < n)
    p <<= 1;
  return p;
}

int index_init (Allocator * alloc, GJHashTable *ix, size_t feature_count)
{
  /* hash table with a load factor <= 0.5 */
  ix->nb = bin_size (feature_count * 2);
  ix->n  = 0;
  ix->b  = allocate (alloc, (ix->nb) * (sizeof (GJHashEntry *)));
  ix->mask = ix->nb - 1;
  return 0;
}

int index_add (Allocator * alloc, GJHashTable *ix, Ast * g, const char * key)
{
  GJHashEntry * e = allocate (alloc, sizeof (GJHashEntry));
  e->key = key; /* assumes json AST survives till the end */
  e->geometry = g;
  size_t i = hash (key) & ix->mask;
  e->next = ix->b [i];
  ix->b [i] = e;
  ix->n++;
  return 0;
}

# ifndef AST_MAX_DEPTH
  # define AST_MAX_DEPTH 256
# endif
#define TREE_START                            \
  do {                                        \
    int depth = 0;                            \
    do {                                      \
      while (a->child)                        \
      {                                       \

#define TREE_END                              \
        /* push to stack */                   \
        depth++, a = a->child;                \
      }                                       \
      /* pop from stack */                    \
      while (depth)                           \
      {                                       \
        if (a->next) { a = a->next; break; }  \
        depth--, a = a->parent;               \
      }                                       \
    } while (depth);                          \
  } while (0)

static char * find_tag (Ast * prop, const char * tag)
{
  /* error : properties not attached with geometry */
  if (!prop)
  {
    fprintf (stderr, "\nGeoJSON \"properties\" not found alongside \"geometry\"");
    return NULL;
  }

  Ast * a = prop;
  TREE_START
    if (a->type == ast_token (jpair) && 
       !strcmp (((AstLeaf *)a->child)->token, tag))
    {
      Ast * id = a->child->next->child;
      /* expects only a string ID */
      if (!id || id->type != ast_token (jstring))
        return NULL;
      char * token = ((AstLeaf *) id)->token;
      return token;
    }
  TREE_END;

  /* cannot find that tag */
  fprintf (stderr, "\n[unique] tag not found in GeoJSON \"properties\"");
  return NULL;
}

int geojson (const char * json_file, GeoJSON * gj, const char * uniqueTag)
{
  /* parse & create AST */
  if (yyparse_file (&gj->tree, json_file))
    return -1;
  Ast * root = gj->tree;

  /* reusing the ast allocator */
  Allocator * alloc = ast_allocator (root);

  /* count the number of geometry features */ 
  int nf = 0;
  Ast * a = gj->tree;
  assert (a->type == ast_token (jjson));

  TREE_START
    if (a->type == ast_token (jpair) &&
        !strcmp (((AstLeaf *)a->child)->token, "geometry"))
    {
      nf++; break;
    }
    /* push to stack */
    assert (depth < AST_MAX_DEPTH && "ast stack limit");
  TREE_END;

  /* create a hashtable */
  GJHashTable * idx = &gj->hashtable; 
  index_init (alloc, idx, nf);

  /* locate tag in the "properties" and add the tag Id to hashtable */
  a = gj->tree;
  TREE_START
    if (a->type == ast_token (jpair) &&
        !strcmp (((AstLeaf *)a->child)->token, "geometry"))
    {
      assert (depth);
      Ast * parent = a->parent, * prop = NULL, * child = parent->child;

      /* note that order of "type", "geometry", "properties" is not enforced in GJSON*/
      while (child)
      {
        if (child->type == ast_token (jpair) &&
            !strcmp (((AstLeaf *)child->child)->token, "properties"))
        { prop = child; break; }
        child = child->next;
      }

      /* error : cannot find the unique id or the paired [string] ID */
      char * tag_id = find_tag (prop, uniqueTag);
      if (!tag_id)
      {
        ast_free (gj->tree);
        *gj = (GeoJSON) {0};
        return -1;
      }
      index_add (alloc, idx, a, tag_id);
      a = a->parent, depth--;
      break;
    }
  TREE_END;

  printf ("\n nfeatures %d : %zd", nf, idx->n);
  return 0;
}

/* First match; pass the previous result as `after` to walk duplicates. */
const Ast * geojson_locate (GeoJSON * gj, const char * key)
{
  assert (gj && key);
  GJHashTable * ix = & gj->hashtable;
  const GJHashEntry * e = ix->b [hash(key) & ix->mask];
  for (; e; e = e->next)
    if (!strcmp(e->key, key))
      return e->geometry;
  return NULL;
}
