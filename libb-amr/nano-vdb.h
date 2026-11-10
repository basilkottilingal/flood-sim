#ifndef _NANO_VDB_H_
#define _NANO_VDB_H_

  /*
  .. A sparse (1 << level)^dimension array of pointers, organised as a VDB
  .. tree (OpenVDB, NanoVDB): fixed depth, very wide nodes.
  ..
  ..   root    a hash table of upper nodes, keyed by the indices >> 12; one
  ..           upper node, without a table, while len <= 4096
  ..   upper   2^(5 d) slots: bits 11..7 of every index
  ..   lower   2^(4 d) slots: bits  6..3
  ..   leaf    2^(3 d) slots: bits  2..0, holding the data pointers
  ..
  .. In 3D that is VDB's 32^3 / 16^3 / 8^3; in 2D 32^2 / 16^2 / 8^2. The
  .. first index is the most significant in a slot, as in VDB. Every node
  .. keeps a mask of its occupied slots and their count.
  ..
  .. Two layouts of a node's children, chosen when compiling:
  ..   dense (default)   a table of every slot, NULL where empty, as
  ..                     OpenVDB's internal nodes: one load per level
  ..   -DNVDB_COMPACT    only the occupied children, in slot order, found
  ..                     by rank: the bits set before the slot, a running
  ..                     count per 64-bit word of the mask plus a popcount
  ..                     in the word (NanoVDB's index grids). Build with
  ..                     -mpopcnt or -march=native, or each popcount is a
  ..                     call.
  ..
  .. The same interface as memindex.h and voxel-memindex.h, prefixed nvdb_:
  ..   nvdb_get ()          the pointer at an index, or NULL (out of range
  ..                        too), in one walk
  ..   nvdb_get_cursor ()   the same through an accessor, which keeps the
  ..                        leaf, lower and upper node last reached: an
  ..                        index in the same leaf costs one load, in the
  ..                        same lower node two. Reset it ({0}) after any
  ..                        nvdb_free () or nvdb_destroy () of that array.
  ..   nvdb_assign (), nvdb_free (), nvdb_destroy (), nvdb_allocated ()
  .. A node that becomes empty goes back to its pool at once.
  */

  #include <stdint.h>
  #include "pool.h"

  #ifndef dimension
  #  error "dimension must be defined as 1, 2 or 3"
  #elif (dimension < 1) || (dimension > 3)
  #  error "dimension must be 1, 2 or 3"
  #endif

  /* level 0 leaf, 1 lower, 2 upper: 3 + l bits of each index */
  #define NVDB_LOG2(l)   (3 + (l))
  #define NVDB_SHIFT(l)  ((l) == 0 ? 0 : (l) == 1 ? 3 : 7)  /* bits below  */
  #define NVDB_SPAN(l)   (NVDB_SHIFT (l) + NVDB_LOG2 (l))   /* 3, 7, 12    */
  #define NVDB_SLOTS(l)  (1 << (dimension * NVDB_LOG2 (l)))
  #define NVDB_WORDS(l)  (NVDB_SLOTS (l) >= 64 ? NVDB_SLOTS (l) / 64 : 1)
  #define NVDB_ROOT_SPAN 12                                 /* upper: 4096 */

  typedef struct {
    uint32_t n;          /* occupied slots                                  */
    uint32_t cap;        /* compact: room in c                              */
    void ** c;           /* dense: the table after the mask; compact: the
                            occupied children, in slot order (malloc)      */
    uint64_t mask [];    /* NVDB_WORDS (l) words, then dense: the table of
                            NVDB_SLOTS (l) pointers; compact: NVDB_WORDS (l)
                            uint32_t, the set bits of the words before     */
  } NvNode;

  typedef struct {
    uint64_t key;        /* the indices >> 12, 21 bits each                 */
    NvNode * node;       /* NULL: an empty entry                            */
  } NvTile;

  struct _NanoVDB {
    int len;             /* 1 << level; 0 while nothing was ever assigned   */
    NvNode * single;     /* the upper node, while len <= 4096               */
    NvTile * root;       /* len > 4096: open addressing, linear probing    */
    uint32_t root_cap, root_n;
    Pool * pool [3];     /* leaf, lower and upper nodes                     */
  };
  #define NanoVDB struct _NanoVDB *

  /* the indices as arguments, as values, and the range test */
  #if dimension == 1
  #  define NVDB_ARGS     int i
  #  define NVDB_IDX      i
  #  define NVDB_IN(m)    ((unsigned) i < (unsigned) (m)->len)
  #elif dimension == 2
  #  define NVDB_ARGS     int i, int j
  #  define NVDB_IDX      i, j
  #  define NVDB_IN(m)    ((unsigned) i < (unsigned) (m)->len &&          \
                           (unsigned) j < (unsigned) (m)->len)
  #else
  #  define NVDB_ARGS     int i, int j, int k
  #  define NVDB_IDX      i, j, k
  #  define NVDB_IN(m)    ((unsigned) i < (unsigned) (m)->len &&          \
                           (unsigned) j < (unsigned) (m)->len &&          \
                           (unsigned) k < (unsigned) (m)->len)
  #endif

  typedef struct {
    NvNode * node [3];          /* the leaf, lower and upper node last reached */
    int key [3][dimension];     /* their indices >> 3, >> 7, >> 12             */
  } NvdbCursor;

  NanoVDB nvdb_new     (int len);
  void    nvdb_destroy (NanoVDB m, int len);
  void    nvdb_assign  (NanoVDB m, NVDB_ARGS, int len, void * b);
  void    nvdb_free    (NanoVDB m, NVDB_ARGS, int len, void * b);

  /* the slot of indices x in a node of level l: x [0] most significant */
  static inline int nvdb_slot (const int * x, int l)
  {
    int b = NVDB_LOG2 (l), s = NVDB_SHIFT (l), n = 0;
    for (int d = 0; d < dimension; d++)
      n = (n << b) | ((x [d] >> s) & ((1 << b) - 1));
    return n;
  }

  /* the child (or data) in slot s of node p of level l, or NULL */
  static inline void * nvdb_child (const NvNode * p, int l, int s)
  {
  #ifdef NVDB_COMPACT
    int w = s >> 6;
    uint64_t bit = 1ull << (s & 63), word = p->mask [w];
    if (!(word & bit))
      return NULL;                                  /* absent: one bit */
    const uint32_t * prefix = (const uint32_t *) (p->mask + NVDB_WORDS (l));
    return p->c [prefix [w] + __builtin_popcountll (word & (bit - 1))];
  #else
    (void) l;
    return p->c [s];
  #endif
  }

  static inline uint64_t nvdb_key (const int * x)
  {
    uint64_t key = 0;
    for (int d = 0; d < dimension; d++)
      key = key << 21 | (uint64_t) (x [d] >> NVDB_ROOT_SPAN);
    return key;
  }

  static inline uint32_t nvdb_hash (uint64_t key, uint32_t cap)
  {
    return (uint32_t) ((key * 0x9E3779B97F4A7C15ull) >> 32) & (cap - 1);
  }

  static inline NvNode * nvdb_root_find (const struct _NanoVDB * m,
                                         uint64_t key)
  {
    if (!m->root_cap)
      return NULL;
    for (uint32_t h = nvdb_hash (key, m->root_cap); ;
         h = (h + 1) & (m->root_cap - 1)) {
      if (!m->root [h].node)
        return NULL;
      if (m->root [h].key == key)
        return m->root [h].node;
    }
  }

  /* the upper node holding indices x, or NULL */
  static inline NvNode * nvdb_upper (const struct _NanoVDB * m, const int * x)
  {
    if (m->len <= 1 << NVDB_ROOT_SPAN)
      return m->single;
    return nvdb_root_find (m, nvdb_key (x));
  }

  static inline void * nvdb_get (NanoVDB m, NVDB_ARGS)
  {
    if (!NVDB_IN (m))
      return NULL;
    int x [dimension] = { NVDB_IDX };
    NvNode * u = nvdb_upper (m, x);
    if (!u)
      return NULL;
    NvNode * l = nvdb_child (u, 2, nvdb_slot (x, 2));
    if (!l)
      return NULL;
    NvNode * f = nvdb_child (l, 1, nvdb_slot (x, 1));
    if (!f)
      return NULL;
    return nvdb_child (f, 0, nvdb_slot (x, 0));
  }

  static inline int nvdb_allocated (NanoVDB m, NVDB_ARGS)
  {
    return nvdb_get (m, NVDB_IDX) != NULL;
  }
  #define nvdb_data(m, ...)  nvdb_get (m, __VA_ARGS__)

  static inline int nvdb_same (const int * key, const int * x, int span)
  {
    for (int d = 0; d < dimension; d++)
      if (key [d] != x [d] >> span)
        return 0;
    return 1;
  }

  static inline void nvdb_keep (NvdbCursor * c, int l, NvNode * n,
                                const int * x)
  {
    c->node [l] = n;
    for (int d = 0; d < dimension; d++)
      c->key [l][d] = x [d] >> NVDB_SPAN (l);
  }

  /* VDB's ValueAccessor: from the deepest cached node that holds x */
  static inline void * nvdb_get_cursor (NanoVDB m, NvdbCursor * c, NVDB_ARGS)
  {
    if (!NVDB_IN (m))
      return NULL;
    int x [dimension] = { NVDB_IDX };
    if (c->node [0] && nvdb_same (c->key [0], x, NVDB_SPAN (0)))
      return nvdb_child (c->node [0], 0, nvdb_slot (x, 0));
    NvNode * l;
    if (c->node [1] && nvdb_same (c->key [1], x, NVDB_SPAN (1)))
      l = c->node [1];
    else {
      NvNode * u;
      if (c->node [2] && nvdb_same (c->key [2], x, NVDB_SPAN (2)))
        u = c->node [2];
      else {
        if (!(u = nvdb_upper (m, x)))
          return NULL;
        nvdb_keep (c, 2, u, x);
      }
      if (!(l = nvdb_child (u, 2, nvdb_slot (x, 2))))
        return NULL;
      nvdb_keep (c, 1, l, x);
    }
    NvNode * f = nvdb_child (l, 1, nvdb_slot (x, 1));
    if (!f)
      return NULL;
    nvdb_keep (c, 0, f, x);
    return nvdb_child (f, 0, nvdb_slot (x, 0));
  }

#endif
