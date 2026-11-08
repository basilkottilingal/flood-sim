/*
.. voxel-memindex: a sparse array of pointers as a tree of 64-way nodes.
.. See voxel-memindex.h.
*/

#include <stdlib.h>
#include <assert.h>

#include "voxel-memindex.h"

VMemindex vmem_new (int len)
{
  (void) len;
  return calloc (1, sizeof (struct _VMemindex));
}

/* every node goes with the pool; the struct itself stays, as in memindex */
void vmem_destroy (VMemindex m, int len)
{
  (void) len;
  pool_destroy (m->pool);
  *m = (struct _VMemindex) {0};
}

static void start (VMemindex m, int len)
{
  assert (len > 0 && !(len & (len - 1)));
  int level = __builtin_ctz ((unsigned) len);
  int depth = level ? (level + VMEM_BITS - 1) / VMEM_BITS : 1;
  m->len = len;
  m->top = (depth - 1) * VMEM_BITS;
  m->pool = pool_new (sizeof (VNode));
  m->root = pool_alloc (m->pool);
}

void vmem_assign (VMemindex m, VMEM_ARGS, int len, void * b)
{
  assert (b != NULL);
  if (!m->root)
    start (m, len);
  assert (VMEM_IN (m));
  VNode * n = m->root;
  for (int s = m->top; s > 0; s -= VMEM_BITS) {
    void ** c = &n->c [VMEM_DIGIT (s)];
    if (!*c)
      *c = pool_alloc (m->pool);      /* zeroed */
    n = *c;
  }
  n->c [VMEM_DIGIT (0)] = b;
}

static int empty (const VNode * n)
{
  for (int d = 0; d < VMEM_NCHILDREN; d++)
    if (n->c [d])
      return 0;
  return 1;
}

void vmem_free (VMemindex m, VMEM_ARGS, int len, void * b)
{
  (void) len;
  (void) b;
  if (!VMEM_IN (m))
    return;
  /* the path from the root: path [h] and the digit taken there */
  VNode * path [32];
  int digit [32], h = 0;
  VNode * n = m->root;
  for (int s = m->top; ; s -= VMEM_BITS) {
    path [h] = n;
    digit [h++] = VMEM_DIGIT (s);
    if (s == 0)
      break;
    if (!(n = n->c [digit [h - 1]]))
      return;
  }
  if (!path [h - 1]->c [digit [h - 1]])
    return;
  path [h - 1]->c [digit [h - 1]] = NULL;
  /* empty nodes go back to the pool, all but the root */
  for (int l = h - 1; l > 0 && empty (path [l]); l--) {
    pool_free (m->pool, path [l]);
    path [l - 1]->c [digit [l - 1]] = NULL;
  }
}
