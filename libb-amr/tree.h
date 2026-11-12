#ifndef _BAMR_H_
#define _BAMR_H_

  #include "level-index.h"
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
  /*
  .. Every block stores GHOSTS ghost layers on each side, kept up to date by
  .. ghosts_update () (block-ghosts.h). At least 2: adapt_wavelet () reads
  .. two layers. At most m / 2.
  ..
  .. BAMR_ALIGN: the vector width, in Reals, that every interior row of a
  .. stored block starts on (block_shape_aligned () in block-loops.h). 1, the
  .. default, is the compact layout, no padding.
  */
  #ifndef GHOSTS
    # define GHOSTS 2
  #endif
  #ifndef BAMR_ALIGN
    # define BAMR_ALIGN 1
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

  extern LevelIndex grid_info_layers [BAMR_MAX_DEPTH + 1];

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
    BAMR_NODE_REFINE = 8,         /* needs refining: adapt.h          */
    BAMR_NODE_COARSEN = 16,       /* may be coarsened, with its family */
  } BAMRNodeFlags;

  typedef struct BAMRNode_
  {
    int flags;
    int pid;
    void *** s;     /* s [id/w][id%w]: the block of scalar id, or NULL */
  } BAMRNode;

  /*
  .. The tree's nodes in lists, made again only after the tree changes
  .. (tree (), tree_refine (), tree_coarsen (), tree_destroy (), and so
  .. adapt_wavelet ()): what the loops over the tree read, instead of
  .. walking it each time. Each list is in Morton order (morton-curve.h: a
  .. leaf as it is reached, a parent after its children), so that a loop over
  .. one visits the nodes as the walk did.
  ..
  ..   tree_caches ()   the lists, made again first if the tree changed
  ..   tree_generation  changes with the tree: a loop over a list checks it,
  ..                    and stops the run (tree_changed_under ()) if the
  ..                    tree changed during the loop: refine or coarsen
  ..                    after the loop, not in it.
  */
  typedef struct
  {
    Index p;
    struct BAMRNode_ * node;
  } TreeItem;

  typedef struct
  {
    const TreeItem * all;          /* every node                       */
    int nall;
    const TreeItem * leaves;       /* the leaves                       */
    int nleaves;
    const TreeItem * by_level;     /* every node of level l: by_level
                                      [level_first [l] .. level_first
                                      [l + 1])                         */
    int level_first [BAMR_MAX_DEPTH + 2];
    const TreeItem * leaves_by_level;  /* the same, leaves only        */
    int leaf_first [BAMR_MAX_DEPTH + 2];
    int depth;                     /* the deepest level with a node    */
  } TreeCaches;

  const TreeCaches * tree_caches (void);
  extern unsigned tree_generation;
  void tree_changed_under (const char * loop);

  /*
  .. tree_node (IJKL): the node there, or NULL. tree_is_allocated (IJKL):
  .. whether there is one. Both through the level's index (level-index.h).
  */
  #if dimension == 1
    #define tree_node(i, level)                                   \
      ((BAMRNode *) lidx_get (&grid_info_layers [level], i))
    #define tree_is_allocated(i, level)                           \
      ((level) <= BAMR_MAX_DEPTH &&                               \
        lidx_allocated (&grid_info_layers [level], i))
  #elif dimension == 2
    #define tree_node(i, j, level)                                \
      ((BAMRNode *) lidx_get (&grid_info_layers [level], i, j))
    #define tree_is_allocated(i, j, level)                        \
      ((level) <= BAMR_MAX_DEPTH &&                               \
        lidx_allocated (&grid_info_layers [level], i, j))
  #else
    #define tree_node(i, j, k, level)                             \
      ((BAMRNode *) lidx_get (&grid_info_layers [level], i, j, k))
    #define tree_is_allocated(i, j, k, level)                     \
      ((level) <= BAMR_MAX_DEPTH &&                               \
        lidx_allocated (&grid_info_layers [level], i, j, k))
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
      BAMRNode * node = level <= BAMR_MAX_DEPTH ?                 \
        TREE_CALL (tree_node, IJKL) : NULL;     /* one lookup */  \
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
  double tree_get_origin (int axis); /* x0, y0, z0: the domain's corner  */
  double tree_get_size (void);       /* h0: its side                     */
  void tree_destroy (void);
  int  tree_get_depth (void);       /* the finest level of a leaf: the
                                       deepest with a node             */
  int  tree_get_depth_max (void);   /* the depth tree () was given     */
  int  tree_get_depth_min (void);   /* the coarsest level of a leaf    */
                                    /* both from a bit set of the levels
                                       with leaves, O(1): BAMR_MAX_DEPTH
                                       must be below 32 */
  int  tree_get_nblocks (void);     /* nodes, leaves or not            */
  int  tree_get_nleaves (void);     /* leaf nodes                      */
  int  tree_get_block_m (void);     /* interior cells per block side,
                                       2^blockdepth; a block stores them
                                       and GHOSTS layers around:
                                       tree_block_shape () */

  /*
  .. Refining and coarsening, keeping the tree 2:1 balanced: two leaves that
  .. touch, by a face, an edge or a corner, differ by at most one level.
  ..
  ..   tree_refine (p)    divides the leaf p into its TREE_NCHILDREN
  ..                      children, refining first any coarser leaf beside
  ..                      it, as far as needed. Returns how many leaves were
  ..                      divided: 0 if p is not a leaf, or at depth max.
  ..   tree_coarsen (p)   merges the children of p, all leaves, back into p,
  ..                      unless a neighbour of p has grandchildren, which
  ..                      would then touch p two levels finer. Returns 1 if
  ..                      it did, 0 if not.
  ..
  .. The children of a node are allocated together, nodes and blocks: one
  .. object of TREE_NCHILDREN nodes, and for each scalar one object of
  .. TREE_NCHILDREN blocks. A refinement fills the children's interiors by
  .. each scalar's prolongation from p, its interior and ghosts, then the
  .. children's ghosts: so their blocks are not zeroed first. A coarsening
  .. gives p each scalar's restriction of its children, then fills the
  .. ghosts of p and of the leaves beside it that read p's children.
  ..
  .. So both keep every leaf's ghosts valid, if they were before: call
  .. ghosts_update () after changing values and before refining. See
  .. block-ghosts.h.
  */
  int  tree_refine (Index p);
  int  tree_coarsen (Index p);

#endif
