#ifndef _VOXEL_MEMINDEX_H_
#define _VOXEL_MEMINDEX_H_

  /*
  .. A sparse (1 << level)^dimension array of pointers, as a tree of
  .. 64-way nodes: each node takes VMEM_BITS bits of every index (6 in 1D,
  .. 3 in 2D, 2 in 3D), so a level-L array is ceil (L / VMEM_BITS) nodes
  .. deep (at least one). The same interface as memindex.h, prefixed vmem_.
  ..
  .. Memory follows what is stored: one 512-byte node per occupied group
  .. of 64, plus the nodes above it, all from a pool. Nothing depends on
  .. the page size, transparent huge pages or overcommit, unlike memindex,
  .. which pays a page (and a page-table page) per row touched.
  ..
  .. vmem_get ()      the pointer at an index, or NULL (out of range too):
  ..                  one walk, where memindex needs mem_allocated () and
  ..                  then mem_data ().
  .. vmem_get_cursor () the same, but remembers the last bottom node: an
  ..                  index in the same group of 64 costs one load. For
  ..                  sweeps over neighbours. A cursor belongs to one array,
  ..                  starts as {0}, and must be reset ({0}) after any
  ..                  vmem_free () or vmem_destroy () of that array.
  ..
  .. The root is kept until vmem_destroy (); empty nodes below it go back
  .. to the pool as soon as their last entry is freed.
  */

  #include <stdint.h>
  #include "pool.h"

  #ifndef dimension
  #  error "dimension must be defined as 1, 2 or 3"
  #elif (dimension < 1) || (dimension > 3)
  #  error "dimension must be 1, 2 or 3"
  #endif

  #if dimension == 1
  #  define VMEM_BITS  6
  #elif dimension == 2
  #  define VMEM_BITS  3
  #else
  #  define VMEM_BITS  2
  #endif
  #define VMEM_MASK      ((1 << VMEM_BITS) - 1)
  #define VMEM_NCHILDREN 64

  typedef struct _VNode {
    void * c [VMEM_NCHILDREN];  /* children, or the data at the bottom */
  } VNode;

  struct _VMemindex {
    int len;          /* 1 << level; 0 while nothing was ever assigned */
    int top;          /* shift of the root's digit, a multiple of VMEM_BITS */
    VNode * root;
    Pool * pool;
  };
  #define VMemindex struct _VMemindex *

  /* the indices as arguments, the indices as values, and the range test */
  #if dimension == 1
  #  define VMEM_ARGS     int i
  #  define VMEM_IDX      i
  #  define VMEM_IN(m)    ((unsigned) i < (unsigned) (m)->len)
  #  define VMEM_DIGIT(s) (((i) >> (s)) & VMEM_MASK)
  #elif dimension == 2
  #  define VMEM_ARGS     int i, int j
  #  define VMEM_IDX      i, j
  #  define VMEM_IN(m)    ((unsigned) i < (unsigned) (m)->len &&          \
                           (unsigned) j < (unsigned) (m)->len)
  #  define VMEM_DIGIT(s) ((((i) >> (s)) & VMEM_MASK) |                   \
                           (((j) >> (s)) & VMEM_MASK) << VMEM_BITS)
  #else
  #  define VMEM_ARGS     int i, int j, int k
  #  define VMEM_IDX      i, j, k
  #  define VMEM_IN(m)    ((unsigned) i < (unsigned) (m)->len &&          \
                           (unsigned) j < (unsigned) (m)->len &&          \
                           (unsigned) k < (unsigned) (m)->len)
  #  define VMEM_DIGIT(s) ((((i) >> (s)) & VMEM_MASK) |                   \
                           (((j) >> (s)) & VMEM_MASK) << VMEM_BITS |      \
                           (((k) >> (s)) & VMEM_MASK) << (2*VMEM_BITS))
  #endif

  typedef struct {
    VNode * leaf;           /* the bottom node last reached, or NULL */
    int key [dimension];    /* its indices >> VMEM_BITS */
  } VMemCursor;

  VMemindex vmem_new     (int len);
  void      vmem_destroy (VMemindex m, int len);
  void      vmem_assign  (VMemindex m, VMEM_ARGS, int len, void * b);
  void      vmem_free    (VMemindex m, VMEM_ARGS, int len, void * b);

  static inline void * vmem_get (VMemindex m, VMEM_ARGS)
  {
    if (!VMEM_IN (m))
      return NULL;
    VNode * n = m->root;
    for (int s = m->top; s > 0; s -= VMEM_BITS)
      if (!(n = n->c [VMEM_DIGIT (s)]))
        return NULL;
    return n->c [VMEM_DIGIT (0)];
  }

  static inline int vmem_allocated (VMemindex m, VMEM_ARGS)
  {
    return vmem_get (m, VMEM_IDX) != NULL;
  }
  #define vmem_data(m, ...)  vmem_get (m, __VA_ARGS__)

  static inline void * vmem_get_cursor (VMemindex m, VMemCursor * c,
                                        VMEM_ARGS)
  {
    if (!VMEM_IN (m))
      return NULL;
  #if dimension == 1
    int same = c->key [0] == i >> VMEM_BITS;
  #elif dimension == 2
    int same = c->key [0] == i >> VMEM_BITS && c->key [1] == j >> VMEM_BITS;
  #else
    int same = c->key [0] == i >> VMEM_BITS && c->key [1] == j >> VMEM_BITS &&
      c->key [2] == k >> VMEM_BITS;
  #endif
    if (c->leaf && same)
      return c->leaf->c [VMEM_DIGIT (0)];
    VNode * n = m->root;
    for (int s = m->top; s > 0; s -= VMEM_BITS)
      if (!(n = n->c [VMEM_DIGIT (s)]))
        return NULL;
    c->leaf = n;
    c->key [0] = i >> VMEM_BITS;
  #if dimension > 1
    c->key [1] = j >> VMEM_BITS;
  #endif
  #if dimension > 2
    c->key [2] = k >> VMEM_BITS;
  #endif
    return n->c [VMEM_DIGIT (0)];
  }

#endif
