#ifndef _BAMR_H_
#define _BAMR_H_

  #include "memindex.h"
  #include "real.h"

  /*
  .. 2:1 balanced, [block based] Adaptive Mesh Refinement (AMR), for
  .. dimension 1, 2 or 3 (a binary tree, a quadtree or an octree).
  .. Note that the AMR grid is a global grid, and you can't
  .. create multiple tree in a single run
  */

  #define BLOCK_AMR        1
  #define BALANCE_2_1      1
  #ifndef BAMR_MAX_DEPTH
  #define BAMR_MAX_DEPTH   20
  #endif
  #ifndef GHOSTS
    # define GHOSTS 2
  #endif

  /* children of a node */
  #define TREE_NCHILDREN   (1 << dimension)

  #if dimension == 1
    # define IJK           i
    # define IJKargs       i
    # define TREE_ROOT     0, 0
  #elif dimension == 2
    # define IJK           i, j
    # define IJKargs       i, int j
    # define TREE_ROOT     0, 0, 0
  #else
    # define IJK           i, j, k
    # define IJKargs       i, int j, int k
    # define TREE_ROOT     0, 0, 0, 0
  #endif
  #define IJKL             IJK, level
  #define IJKLargs         IJKargs, int level

  extern struct _Memindex grid_info_layers [BAMR_MAX_DEPTH + 1];

  typedef struct
  {
    int IJKL;
  } Index;

  typedef struct
  {
    Index * indices;
    int n;
  } IndexStack;

  typedef enum
  {
    BAMR_NODE_LEAF   = 1,
    BAMR_NODE_IS_PRUNED = 2,
    BAMR_NODE_HAS_CHILDREN = 4,
  } BAMRNodeFlags;

  typedef struct
  {
    int flags;
    int pid;
    void *** s;     /* s [id/w][id%w]: the block of scalar id, or NULL */
  } BAMRNode;

  #if dimension == 1
    #define tree_node(i, level)                                   \
      ((BAMRNode *) mem_data (&grid_info_layers [level], i))
    #define tree_is_allocated(i, level)                           \
      ((level) <= BAMR_MAX_DEPTH &&                               \
        mem_allocated (&grid_info_layers [level], i))
  #elif dimension == 2
    #define tree_node(i, j, level)                                \
      ((BAMRNode *) mem_data (&grid_info_layers [level], i, j))
    #define tree_is_allocated(i, j, level)                        \
      ((level) <= BAMR_MAX_DEPTH &&                               \
        mem_allocated (&grid_info_layers [level], i, j))
  #else
    #define tree_node(i, j, k, level)                             \
      ((BAMRNode *) mem_data (&grid_info_layers [level], i, j, k))
    #define tree_is_allocated(i, j, k, level)                     \
      ((level) <= BAMR_MAX_DEPTH &&                               \
        mem_allocated (&grid_info_layers [level], i, j, k))
  #endif

  #define tree_has_children(n)                                    \
    (((n)->flags & BAMR_NODE_HAS_CHILDREN) &&                     \
    !((n)->flags & BAMR_NODE_IS_PRUNED))

  /* f (a, b, ...) after the arguments are expanded, so IJKL can be passed */
  #define TREE_CALL(f, ...)  f (__VA_ARGS__)

  /*
  .. In the scope of IJKL (int i, j, k as the dimension has them, and level):
  .. node, the node there or NULL, and sblocks, its scalar blocks or NULL.
  */
  #define BAMR_NODE_ENV()                                         \
      BAMRNode * node = TREE_CALL (tree_is_allocated, IJKL) ?     \
        TREE_CALL (tree_node, IJKL) : NULL;                       \
      void *** sblocks = node ? node->s : NULL;                   \
      (void) sblocks;

  /*
  .. Index helpers, the same in every dimension.
  ..   index_child (p, c)    child c (0 to TREE_NCHILDREN - 1) of p; bit 0 of
  ..                         c steps in i, bit 1 in j, bit 2 in k
  ..   index_is_allocated (p), index_node (p)
  ..   INDEX_UNPACK (p);     declares IJKL from p, for BAMR_NODE_ENV ()
  */
  #if dimension == 1
    #define INDEX_IJK(p)     (p).i
    #define INDEX_UNPACK(p)                                       \
      int i = (p).i, level = (p).level;                           \
      (void) i; (void) level
  #elif dimension == 2
    #define INDEX_IJK(p)     (p).i, (p).j
    #define INDEX_UNPACK(p)                                       \
      int i = (p).i, j = (p).j, level = (p).level;                \
      (void) i; (void) j; (void) level
  #else
    #define INDEX_IJK(p)     (p).i, (p).j, (p).k
    #define INDEX_UNPACK(p)                                       \
      int i = (p).i, j = (p).j, k = (p).k, level = (p).level;     \
      (void) i; (void) j; (void) k; (void) level
  #endif
  #define index_is_allocated(p)                                   \
    TREE_CALL (tree_is_allocated, INDEX_IJK (p), (p).level)
  #define index_node(p)                                           \
    TREE_CALL (tree_node, INDEX_IJK (p), (p).level)

  static inline Index index_child (Index p, int c)
  {
    Index q = p;
    q.level = p.level + 1;
    q.i = 2 * p.i + (c & 1);
    #if dimension > 1
    q.j = 2 * p.j + ((c >> 1) & 1);
    #endif
    #if dimension > 2
    q.k = 2 * p.k + ((c >> 2) & 1);
    #endif
    return q;
  }

  /*
  .. Loops over the indices around a node, cwrap @loops. Each gives its block
  .. IJKL, node and sblocks of the index visited, as BAMR_NODE_ENV () does:
  .. node is NULL where no node is allocated, as outside the domain. Each is
  .. one loop, so 'break' leaves it and 'continue' goes to the next index, in
  .. any dimension.
  ..
  ..   children_of (p)              the TREE_NCHILDREN children of p, in Morton
  ..                                order: index_child (p, c), c = 0, 1, ...
  ..   neighborhood_of (p, radius)  every index of p's level within radius of p
  ..                                in each direction, p included: the box of
  ..                                (2 radius + 1)^dimension, i fastest
  ..   neighbors_of (p, radius)     the same box, without p
  */
  @loop children_of (Index parent) : IJKL, node, sblocks
  {
    for (int c = 0; c < TREE_NCHILDREN; c++)
    {
      Index q = index_child (parent, c);
      INDEX_UNPACK (q);
      BAMR_NODE_ENV ();
      {...}
    }
  }

  /* the m-th index of the box of width w around p, m from 0 to w^dimension */
  static inline Index index_in_box (Index p, int radius, int m)
  {
    int w = 2 * radius + 1;
    p.i += m % w - radius;
    #if dimension > 1
    p.j += m / w % w - radius;
    #endif
    #if dimension > 2
    p.k += m / (w * w) - radius;
    #endif
    return p;
  }

  static inline int box_size (int radius)
  {
    int w = 2 * radius + 1, n = w;
    #if dimension > 1
    n *= w;
    #endif
    #if dimension > 2
    n *= w;
    #endif
    return n;
  }

  @loop neighborhood_of (Index center, int radius) : IJKL, node, sblocks
  {
    for (int m = 0, size = box_size (radius); m < size; m++)
    {
      Index q = index_in_box (center, radius, m);
      INDEX_UNPACK (q);
      BAMR_NODE_ENV ();
      {...}
    }
  }

  @loop neighbors_of (Index center, int radius) : IJKL, node, sblocks
  {
    for (int m = 0, size = box_size (radius); m < size; m++)
    {
      if (m == size / 2)          /* the center: the middle of the box */
        continue;
      Index q = index_in_box (center, radius, m);
      INDEX_UNPACK (q);
      BAMR_NODE_ENV ();
      {...}
    }
  }

  /*
  .. API (s)
  */
  int  tree (int depthmax, int blockdepth);
  void tree_shift_resize (double x0, double y0, double z0, double h0);
  void tree_destroy (void);
  int  tree_depth (void);       /* deepest level with a node */
  int  tree_nblocks (void);     /* nodes, leaves or not      */
  int  tree_nleaves (void);     /* leaf nodes                */
  int  tree_block_m (void);     /* cells per block side: 2^blockdepth; a
                                   block stores m^dimension, no ghosts */

#endif
