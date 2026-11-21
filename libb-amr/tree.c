#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <stdint.h>
#include <assert.h>

#include "bamr-error.h"
#include "tree.h"
#include "grid-scalar.h"
#include "block-ghosts.h"
#include "pool.h"
#include "libutils/bits.h"
#ifdef _MPI
# include <mpi.h>
# include "tree-mpi.h"
#endif

#define G(x) grid_info_##x

/* the levels that hold terminal cells are the bits of a uint32_t */
#if BAMR_MAX_DEPTH >= 32
# error "BAMR_MAX_DEPTH must be below 32: tree.c keeps levels as bits of a uint32_t"
#endif
#if dimension * BAMR_MAX_DEPTH > 64
# error "dimension * BAMR_MAX_DEPTH above 64: tree_morton_key () is a uint64_t"
#endif

/*
.. Leaves per level, and as bits: bit l of level_bits is set while level l
.. has any. The coarsest level of a leaf is then the lowest bit, and the
.. depth the highest. Pruned nodes are not leaves.
*/
static double G(x0) = 0.;     /* left corner of domain            */
static double G(y0) = 0.;     /* bottom corner of domain          */
static double G(z0) = 0.;
static double G(h0) = 1.;     /* domain [x0, x0+h] x [y0, y0+h]   */
static int G(nb) = 0;         /* number of blocks                 */
static int G(nlb) = 0;        /* number of leaf blocks            */
static int G(depthmax) = 1;   /* max depth forced by user         */
static int G(blockdepth) = 1; /* 2^blockdepth cells per block side */
static int G(nleaves) [BAMR_MAX_DEPTH + 1];
static uint32_t G(level_bits) = 0;
static int G(npruned) = 0;    /* pruned nodes: tree_prune ()      */
static int G(depth0) = 0;     /* tree ()'s level, 0: its default  */

/* the processes (tree.h): MPI_COMM_WORLD's, read by tree (); or constants */
#ifdef _MPI
static int G(npe) = 1, G(rank) = 0;
#else
static const int G(npe) = 1, G(rank) = 0;
#endif

int tree_npe (void)
{
  return G(npe);
}

int tree_rank (void)
{
  return G(rank);
}

/*
.. *v, of type (TREE_INT ... TREE_DOUBLE), replaced by op (TREE_SUM,
.. TREE_PROD, TREE_MIN, TREE_MAX) over every process: one MPI_Allreduce.
.. A no-op on one process. cwrap's reduction clause calls it after a
.. foreach (); the value is copied first, as MPI forbids the same buffer
.. for input and output.
*/
void tree_reduce (void * v, int type, int op)
{
  assert (type >= TREE_INT && type <= TREE_DOUBLE);
  assert (op >= TREE_SUM && op <= TREE_MAX);
  #ifdef _MPI
  if (G(npe) > 1)
  {
    static const MPI_Datatype types [] = {
      MPI_INT, MPI_LONG, MPI_LONG_LONG, MPI_UNSIGNED, MPI_UNSIGNED_LONG,
      MPI_UNSIGNED_LONG_LONG, MPI_FLOAT, MPI_DOUBLE,
    };
    static const size_t size [] = {
      sizeof (int), sizeof (long), sizeof (long long), sizeof (unsigned),
      sizeof (unsigned long), sizeof (unsigned long long), sizeof (float),
      sizeof (double),
    };
    MPI_Op ops [] = { MPI_SUM, MPI_PROD, MPI_MIN, MPI_MAX };
    char in [sizeof (long double)];
    memcpy (in, v, size [type]);
    MPI_Allreduce (in, v, 1, types [type], ops [op], MPI_COMM_WORLD);
  }
  #else
  (void) v;
  #endif
}
LevelIndex G(layers) [BAMR_MAX_DEPTH + 1] = {{0}};
static Pool * pool_node = NULL;       /* the root, which has no siblings */

/*
.. The lists of the tree's nodes (tree.h), made again by one walk when read
.. after a change. Their arrays only grow: no allocation for each update.
*/
static int caches_are_dirty = 1;
unsigned tree_generation = 0;
static TreeItem * c_all, * c_leaves, * c_by_level, * c_leaves_by_level,
  * c_parents_by_level, * c_owned;
static int c_cap = 0;
static TreeCaches caches;

/* the tree changes: at each node made or freed, and with tree (), tree_destroy () */
static void caches_changed (void)
{
  caches_are_dirty = 1;
  tree_generation++;
}
static Pool * pool_family = NULL;     /* TREE_NCHILDREN nodes at once    */

/*
.. Interior cells per block side, m = 2^blockdepth. A block stores them and
.. GHOSTS layers around: tree_block_shape ().cells Reals (block-loops.h).
*/
#define bm           (1 << G (blockdepth))

#define is_leaf()    (node->flags & BAMR_NODE_LEAF)
#define is_pruned()  (node->flags & BAMR_NODE_IS_PRUNED)
#define is_parent()  (node->flags & BAMR_NODE_HAS_CHILDREN)

/*
.. n leaves more (or fewer, n < 0) at level l: every change of the counts
.. goes through here, so that level_bits follows them.
*/
static void terminal_add (int l, int n)
{
  assert (l >= 0 && l <= BAMR_MAX_DEPTH);
  G(nleaves) [l] += n;
  assert (G(nleaves) [l] >= 0);
  if (G(nleaves) [l])
    G(level_bits) |= (uint32_t) 1 << l;
  else
    G(level_bits) &= ~((uint32_t) 1 << l);
}

/*
.. The tree gone, back to the state before tree (): the MPI plan first (it
.. points into the tree), the caches, the level index, then every node and
.. block with their pools; the domain and the counts reset. MPI itself
.. stays up: a new tree () may follow.
*/
void tree_destroy (void)
{
  #ifdef _MPI
  tree_mpi_end ();                    /* its plan points into the tree */
  #endif
  caches_changed ();
  free (c_all), free (c_leaves), free (c_by_level), free (c_leaves_by_level);
  free (c_parents_by_level), free (c_owned);
  c_all = c_leaves = c_by_level = c_leaves_by_level = c_parents_by_level = NULL;
  c_owned = NULL;
  c_cap = 0;
  for (int l = 0; l <= BAMR_MAX_DEPTH; ++l)
  {
    lidx_destroy (& G(layers) [l]);
    G(nleaves) [l] = 0;
  }
  /* every node and block goes with its pool */
  grid_scalar_destroy ();
  pool_destroy (pool_family);
  pool_destroy (pool_node);
  pool_node = pool_family = NULL;
  G(x0) = 0.;
  G(y0) = 0.;
  G(z0) = 0.;
  G(h0) = 1.;
  G(nb) = 0;
  G(nlb) = 0;
  G(level_bits) = 0;
  G(npruned) = 0;
  G(depthmax) = 1;
  G(blockdepth) = 1;
}

/* the root: a node and blocks of its own */
static BAMRNode * allocate_root (void)
{
  BAMRNode * node = pool_alloc (pool_node);   /* zeroed: flags 0 */
  node->pid = G(rank);
  node->s = grid_scalar_table_new ();
  grid_scalar_root_blocks (node->s);
  return node;
}

/*
.. The family of children of the leaf at IJKL made, below depth max: one
.. object for the 2^D nodes, one for each scalar's blocks (not zeroed:
.. the caller fills them), the pid the parent's. The leaf counts and the
.. LEAF flags are the caller's (divide_leaf_ghosts ()), and so is depth
.. max. 0 if the split fails: no memory, a child already there or not in
.. the level index.
*/
static int add_children (int IJKLargs)
{
  assert (level >= tree_get_depth_min () && level < G(depthmax));
  BAMR_NODE_ENV ();
  assert (node && !is_parent ());
  caches_changed ();
  node->flags |= BAMR_NODE_HAS_CHILDREN;
  LevelIndex * m = & G (layers) [level + 1];
  int len = 1 << (level + 1);
  /*
  .. The family at once: its nodes in one object, child c the c-th, and for
  .. each scalar its blocks in one object. The blocks are not zeroed:
  .. divide_leaf () fills every cell.
  */
  BAMRNode * family = pool_alloc (pool_family);     /* zeroed: flags 0 */
  if (!family)
    return 0;
  void *** t [TREE_NCHILDREN];
  for (int c = 0; c < TREE_NCHILDREN; c++)
  {
    t [c] = family [c].s = grid_scalar_table_new ();
    if (!t [c])
      return 0;
    family [c].pid = node->pid;       /* the owner's, until partitioned */
  }
  grid_scalar_family_new (t, 1);
  int c = 0;
  children_of ((Index) {IJKL})
  {
    if (node)
      return 0;
    lidx_assign (m, IJK, len, &family [c++]);
  }
  Index p = {IJKL};                   /* not in the loop: IJKL holds k */
  for (int q = 0; q < TREE_NCHILDREN; q++)
    if (!index_is_allocated (index_child (p, q)))
      return 0;
  G (nb) += TREE_NCHILDREN;
  return 1;
}

/* the leaf at IJKL divided, its children's interiors prolonged from it;
.. their ghosts made too if ghosts (else the caller's, as tree_divide ()).
.. 0 at depth max; below it the split must succeed, or the run stops */
static int divide_leaf_ghosts (int IJKLargs, int ghosts)
{
  if (level >= G(depthmax))
    return 0;
  if (!add_children (IJKL))
    BAMR_ABORT ("divide_leaf_ghosts (): the split failed (no memory, or a "
                "child already there or not in the level index)");
  BAMR_NODE_ENV ();
  assert (is_leaf ());
  node->flags &= ~BAMR_NODE_LEAF;
  children_of ((Index) {IJKL})
    node->flags |= BAMR_NODE_LEAF;
  terminal_add (level, -1);
  terminal_add (level + 1, TREE_NCHILDREN);
  G (nlb) += TREE_NCHILDREN - 1;
  /*
  .. The children's interiors: each scalar's prolongation from this node,
  .. its interior and ghosts. Then their ghosts: from their siblings, the
  .. nodes of their level around, and the leaves of this level, whose
  .. ghosts are valid. The leaves around that read this node's cells, by a
  .. copy or a prolongation, read the same values from the children now.
  */
  foreach_scalar_ordered (s)
    block_prolong ((Index) {IJKL}, s);
  if (ghosts)
    children_of ((Index) {IJKL})
      foreach_scalar_ordered (s)
        block_ghosts ((Index) {IJKL}, s);
  return 1;
}

static int divide_leaf (int IJKLargs)
{
  return divide_leaf_ghosts (IJKL, 1);
}

/* the children of the leaf-parent at IJKL, gone: it is a leaf again */
static void remove_children (int IJKLargs)
{
  BAMR_NODE_ENV ();
  assert (node && is_parent ());
  caches_changed ();
  BAMRNode * parent = node, * family = NULL;
  void *** t [TREE_NCHILDREN];
  int c = 0;
  children_of ((Index) {IJKL})
  {
    assert (node && (node->flags & BAMR_NODE_LEAF));
    if (!c)
      family = node;                  /* child 0 starts the family's object */
    t [c++] = node->s;
    lidx_free (& G (layers) [level], IJK, 1 << level, node);
  }
  grid_scalar_family_free (t);
  pool_free (pool_family, family);
  parent->flags = (parent->flags & ~BAMR_NODE_HAS_CHILDREN) | BAMR_NODE_LEAF;
  G (nb) -= TREE_NCHILDREN;
  G (nlb) -= TREE_NCHILDREN - 1;
  terminal_add (level + 1, - TREE_NCHILDREN);
  terminal_add (level, 1);
}

/* is q inside the domain at its level? */
static int inside (Index q)
{
  int n = 1 << q.level;
  int in = q.i >= 0 && q.i < n;
  #if dimension > 1
  in &= q.j >= 0 && q.j < n;
  #endif
  #if dimension > 2
  in &= q.k >= 0 && q.k < n;
  #endif
  return in;
}

/* refining and coarsening read the whole tree around a node */
static void not_pruned (const char * f)
{
  if (G(npruned))
  {
    fprintf (stderr, "%s (): the tree is pruned (tree_prune ()): it cannot "
             "be refined or coarsened yet\n", f);
    abort ();
  }
}

/*
.. Refining p, at level l, puts leaves of level l + 1 beside whatever touches
.. p. Each index of level l around p must therefore hold a node already, or
.. the coarser leaf covering it would touch leaves two levels finer: that
.. leaf is refined first, which may in turn refine others, coarser still.
.. By 2:1 balance (faces, edges and corners) the leaf covering a missing
.. index is its parent, one level up, and refining it makes the index: one
.. step, no search (asserted).
*/
int tree_refine (Index p)
{
  not_pruned ("tree_refine");
  if (!index_is_allocated (p) || p.level >= G(depthmax)
      || !(index_node (p)->flags & BAMR_NODE_LEAF))
    return 0;
  int count = 0;
  neighbors_of (p, 1)
  {
    Index q = {IJKL};
    if (inside (q) && !index_is_allocated (q))
    {
      Index c = index_parent (q);     /* the covering leaf, by 2:1 */
      assert (index_is_allocated (c)  /* else the tree is not 2:1 */
              && (index_node (c)->flags & BAMR_NODE_LEAF));
      int n = tree_refine (c);
      assert (n);                     /* coarser than p: below depth max */
      assert (index_is_allocated (q));
      count += n;
    }
  }
  INDEX_UNPACK (p);
  count += divide_leaf (IJKL);
  return count;
}

/* does the node q, one level finer than p, touch p? */
static int touches (Index p, Index q)
{
  int t = q.i >= 2 * p.i - 1 && q.i <= 2 * p.i + 2;
  #if dimension > 1
  t &= q.j >= 2 * p.j - 1 && q.j <= 2 * p.j + 2;
  #endif
  #if dimension > 2
  t &= q.k >= 2 * p.k - 1 && q.k <= 2 * p.k + 2;
  #endif
  return t;
}

/*
.. Does child c of b, a neighbour of p at p's level, touch p? Along each
.. axis where b is beside p, the child must sit in the half of b towards
.. p; along the others, either half. Exact for faces, edges and corners.
*/
int tree_child_touches (Index p, Index b, int c)
{
  /* along each axis where b is beside p, the half of b towards p */
  int d [3] = { b.i - p.i, 0, 0 }, touch = 1;
  #if dimension > 1
  d [1] = b.j - p.j;
  #endif
  #if dimension > 2
  d [2] = b.k - p.k;
  #endif
  for (int a = 0; a < dimension; a++)
    touch &= !d [a] || ((c >> a) & 1) == (d [a] < 0);
  return touch;
}

/*
.. The family below p merged, on one process (a pruned tree aborts):
.. 1. every child must be a leaf;
.. 2. the exact 2:1 rule: no neighbour of p has a child touching p
..    (tree_child_touches ()) that has children of its own, which would
..    touch p two levels finer;
.. 3. p gets the restriction of its children, which are freed;
.. 4. the ghosts again of p and of the finer leaves touching it, which
..    read p's children before and prolong p now.
.. Returns 1 if merged, 0 if not. No cascade: coarsening never breaks 2:1
.. elsewhere.
*/
int tree_coarsen (Index p)
{
  not_pruned ("tree_coarsen");
  if (!index_is_allocated (p) || !tree_has_children (index_node (p)))
    return 0;
  children_of (p)
    if (tree_has_children (node))
      return 0;                       /* a child is not a leaf */
  /* a neighbour's grandchildren would touch p two levels finer: those of
   * its children that touch p (2:1 balance needs no more) */
  neighbors_of (p, 1)
    if (node && tree_has_children (node))
    {
      Index b = {IJKL};
      for (int c = 0; c < TREE_NCHILDREN; c++)
        if (tree_child_touches (p, b, c)
            && tree_has_children (index_node (index_child (b, c))))
          return 0;
    }
  /* p gets the restriction of its children, then they go */
  foreach_scalar_ordered (s)
    block_restrict (p, s);
  {
    INDEX_UNPACK (p);
    remove_children (IJKL);
  }
  /*
  .. The ghosts of p, a leaf now; then of the leaves one level finer that
  .. touch p, which copied p's children and now prolong p. The leaves of
  .. p's level read p, its restriction as before, and the coarser ones
  .. p's parent, unchanged.
  */
  foreach_scalar_ordered (s)
    block_ghosts (p, s);
  neighbors_of (p, 1)
    if (node && tree_has_children (node))
      children_of ((Index) {IJKL})
        if ((node->flags & BAMR_NODE_LEAF) && touches (p, (Index) {IJKL}))
          foreach_scalar_ordered (s)
            block_ghosts ((Index) {IJKL}, s);
  return 1;
}

#ifdef _MPI
/* at exit, if tree () started MPI and nothing finalised it */
static void mpi_finalize (void)
{
  int done = 0;
  MPI_Finalized (&done);
  if (!done)
    MPI_Finalize ();
}
#endif

/*
.. The tree made, depthmax and blockdepth already clamped (tree_replicated
.. ()): MPI started if the program did not (FUNNELED, finalised at exit),
.. npe and rank read; the pools of nodes, families and blocks; the root,
.. a leaf, divided once, so the tree starts at level 1. Aborts if a tree
.. exists already.
*/
static int grid_init (int depthmax, int blockdepth)
{
  if (pool_node)
  {
    fprintf (stderr, "tree: cannot reinit grid; call tree_destroy () first\n");
    exit (-1);
  }
  G (depthmax)   = depthmax;
  G (level_bits) = 0;
  G (npruned)    = 0;
  #ifdef _MPI
  /* MPI, started here if the program did not: finalised at exit */
  int on = 0;
  MPI_Initialized (&on);
  if (!on)
  {
    int provided;
    MPI_Init_thread (NULL, NULL, MPI_THREAD_FUNNELED, &provided);
    atexit (mpi_finalize);
  }
  MPI_Comm_size (MPI_COMM_WORLD, &G(npe));
  MPI_Comm_rank (MPI_COMM_WORLD, &G(rank));
  #endif
  G (blockdepth) = blockdepth;

  assert (GHOSTS >= 1 && 2 * GHOSTS <= bm);
  pool_node = pool_new (sizeof (BAMRNode));
  pool_family = pool_new (TREE_NCHILDREN * sizeof (BAMRNode));
  grid_scalar_init ((size_t) tree_block_shape ().cells);

  Index root = {TREE_ROOT};
  INDEX_UNPACK (root);
  LevelIndex * m = & G (layers) [level];
  BAMRNode * node = allocate_root ();
  node->flags |= BAMR_NODE_LEAF;
  caches_changed ();
  lidx_assign (m, IJK, 1 << level, node);
  G (nb)  = 1;
  G (nlb) = 1;
  terminal_add (level, 1);
  divide_leaf (IJKL);
  return 0;
}

void tree_initial_depth (int level)
{
  G(depth0) = level > 0 ? level : 0;
}

/* every leaf refined, level after level, until level L: a uniform tree */
static void refine_uniform (int L)
{
  for (int l = tree_get_depth_min (); l < L; l++)
  {
    const TreeCaches * tc = tree_caches ();
    int n = tc->nleaves;
    Index * leaves = malloc ((size_t) n * sizeof (Index));
    if (!leaves)
      BAMR_ABORT ("out of memory");
    for (int q = 0; q < n; q++)
      leaves [q] = tc->leaves [q].p;
    for (int q = 0; q < n; q++)
      tree_refine (leaves [q]);
    free (leaves);
  }
}

/*
.. The tree, the one entry point: uniform at level L, the level of
.. tree_initial_depth () (below depth max), by default 1 on one process
.. and on several the first level with at least 4 npe families of leaves
.. (the cut keeps families whole: so every process gets some). On several
.. processes each makes only its part and the skeleton around it, from
.. the cut of a uniform level, a formula (tree_mpi_start ()): no process
.. holds the whole level, none sends its tree. Every process must call it.
*/
int tree (int depthmax, int blockdepth)
{
  tree_replicated (depthmax, blockdepth);
  int L = G(depth0) ? G(depth0) : 1;
  #ifdef _MPI
  if (G(npe) > 1 && !G(depth0))
    while (L < G(depthmax) && (1L << (dimension * (L - 1))) < 4L * G(npe))
      L++;
  #endif
  L = L < G(depthmax) ? L : G(depthmax);
  #ifdef _MPI
  if (G(npe) > 1)
  {
    tree_mpi_start (L);
    return 1;
  }
  #endif
  refine_uniform (L);
  return 1;
}

/*
.. The tree as on one process, also under MPI: every process holds all of
.. it, nothing partitioned (tests, and tree () before its cut). depthmax
.. clamped to 1 ... BAMR_MAX_DEPTH, blockdepth to 4 ... the largest block
.. the dimension allows.
*/
int tree_replicated (int depthmax, int blockdepth)
{
  depthmax = depthmax > BAMR_MAX_DEPTH ? BAMR_MAX_DEPTH :
    depthmax <= 0 ? 1 : depthmax;
  /*
  .. A block has (2^blockdepth + 2 GHOSTS)^dimension Reals, so the largest
  .. blockdepth shrinks with the dimension: 3D at 6 is already 68^3
  .. doubles, 2.4 MiB a block. fixme : tune better
  */
  #if dimension == 1
  enum { blockdepth_max = 16 };
  #elif dimension == 2
  enum { blockdepth_max = 10 };
  #else
  enum { blockdepth_max = 6 };
  #endif
  blockdepth = blockdepth < 4 ? 4 :
    blockdepth > blockdepth_max ? blockdepth_max : blockdepth;
  grid_init (depthmax, blockdepth);
  return 1;
}

/*
.. Every family below p freed, nodes and blocks, deepest first: the nodes
.. freed. The leaves and pruned nodes among them leave the counts.
*/
static int free_below (Index p)
{
  BAMRNode * family = NULL;
  void *** t [TREE_NCHILDREN];
  int n = TREE_NCHILDREN;
  for (int c = 0; c < TREE_NCHILDREN; c++)
  {
    Index q = index_child (p, c);
    BAMRNode * child = index_node (q);
    assert (child);
    if (child->flags & BAMR_NODE_HAS_CHILDREN)
      n += free_below (q);
    if (child->flags & BAMR_NODE_LEAF)
      terminal_add (q.level, -1), G(nlb)--;
    if (child->flags & BAMR_NODE_IS_PRUNED)
      G(npruned)--;
    if (!c)
      family = child;                 /* child 0 starts the family's object */
    t [c] = child->s;
    INDEX_UNPACK (q);
    lidx_free (& G (layers) [level], IJK, 1 << level, child);
  }
  grid_scalar_family_free (t);
  pool_free (pool_family, family);
  G(nb) -= TREE_NCHILDREN;
  return n;
}

/*
.. Hold less: p's subtree freed here, nodes and blocks; p stays, flagged
.. PRUNED (no children here: none, or held by other processes). The tree
.. itself does not change, only what this process holds of it. Returns
.. the nodes freed, 0 if p has no children here.
*/
int tree_prune (Index p)
{
  if (!index_is_allocated (p))
    return 0;
  BAMRNode * node = index_node (p);
  if (!(node->flags & BAMR_NODE_HAS_CHILDREN))
    return 0;
  caches_changed ();
  int n = free_below (p);
  node->flags = (node->flags & ~(BAMR_NODE_HAS_CHILDREN | BAMR_NODE_LEAF))
                | BAMR_NODE_IS_PRUNED;
  G(npruned)++;
  return n;
}

int tree_pruned (void)
{
  return G(npruned);
}

/*
.. A real refinement of this process's leaf p: children made and their
.. interiors prolonged from p. Unlike tree_refine (), no 2:1 cascade and
.. no ghosts: the caller (tree_mpi_refine ()) sees to both. 0 if p is not
.. a leaf here or is at depth max.
*/
int tree_divide (Index p)
{
  if (!index_is_allocated (p) || p.level >= G(depthmax)
      || !(index_node (p)->flags & BAMR_NODE_LEAF))
    return 0;
  INDEX_UNPACK (p);
  return divide_leaf_ghosts (IJKL, 0);
}

/*
.. A real coarsening of this process's family below p: p gets the
.. restriction of its children, which are freed. Unlike tree_coarsen (),
.. no 2:1 check and no ghosts: the caller (tree_mpi_coarsen ()) sees to
.. both. 0 if a child is not a leaf.
*/
int tree_merge (Index p)
{
  if (!index_is_allocated (p)
      || !(index_node (p)->flags & BAMR_NODE_HAS_CHILDREN))
    return 0;
  children_of (p)
    if (!node || !(node->flags & BAMR_NODE_LEAF))
      return 0;
  foreach_scalar_ordered (s)
    block_restrict (p, s);
  INDEX_UNPACK (p);
  remove_children (IJKL);
  return 1;
}

/*
.. Know less: the leaf p, of another process and out of this one's
.. skeleton, is no longer claimed a leaf: PRUNED, its structure unknown,
.. since no one will tell this process when it changes. Its blocks stay
.. (the family's). Nothing if p is not a leaf.
*/
void tree_forget (Index p)
{
  BAMRNode * node = index_node (p);
  if (!node || !(node->flags & BAMR_NODE_LEAF))
    return;
  caches_changed ();
  node->flags = (node->flags & ~BAMR_NODE_LEAF) | BAMR_NODE_IS_PRUNED;
  G(nlb)--;
  G(npruned)++;
  terminal_add (p.level, -1);
}

/*
.. Catch up with another process's division: the leaf p, held here but
.. owned elsewhere, gets children, structure only (inactive, no blocks,
.. the owner's pid). The owner divided it (tree_divide ()) and told this
.. process (tree_mpi_refine ()).
*/
void tree_graft (Index p)
{
  BAMRNode * node = index_node (p);
  assert (node && (node->flags & BAMR_NODE_LEAF) && p.level < BAMR_MAX_DEPTH);
  caches_changed ();
  node->flags = (node->flags & ~BAMR_NODE_LEAF) | BAMR_NODE_HAS_CHILDREN;
  terminal_add (p.level, -1);
  BAMRNode * family = pool_alloc (pool_family);     /* zeroed */
  for (int c = 0; c < TREE_NCHILDREN; c++)
  {
    family [c].s = grid_scalar_table_new ();        /* no blocks: inactive */
    family [c].flags = BAMR_NODE_LEAF | BAMR_NODE_IS_INACTIVE;
    family [c].pid = node->pid;
    Index q = index_child (p, c);
    INDEX_UNPACK (q);
    assert (!index_is_allocated (q));
    lidx_assign (& G (layers) [level], IJK, 1 << level, &family [c]);
  }
  terminal_add (p.level + 1, TREE_NCHILDREN);
  G(nb) += TREE_NCHILDREN;
  G(nlb) += TREE_NCHILDREN - 1;
}

/*
.. Catch up with another process's merge: the family below p, held here
.. but owned elsewhere, freed; p a leaf again. The owner merged it
.. (tree_merge ()) and told this process (tree_mpi_coarsen ()).
*/
void tree_ungraft (Index p)
{
  BAMRNode * node = index_node (p);
  assert (node && (node->flags & BAMR_NODE_HAS_CHILDREN));
  caches_changed ();
  free_below (p);
  node->flags = (node->flags & ~(BAMR_NODE_HAS_CHILDREN | BAMR_NODE_IS_PRUNED))
                | BAMR_NODE_LEAF;
  G(nlb)++;
  terminal_add (p.level, 1);
}

/*
.. Hold more: p made a node of this process's tree, because a message
.. named it (a node moved here, or the skeleton around one). Its missing
.. ancestors are made first, each with its whole family, the new
.. children PRUNED, inactive, without blocks and with pid -1 (the caller
.. sets pids). If leaf and p is pruned, p becomes a leaf. Returns the
.. nodes made.
*/
int tree_insert (Index p, int leaf)
{
  int made = 0;
  if (!index_is_allocated (p))
  {
    /* the parent first, then its family, every child pruned (unknown) */
    assert (p.level > 0);
    Index up = index_parent (p);
    made += tree_insert (up, 0);
    BAMRNode * parent = index_node (up);
    assert (!(parent->flags & (BAMR_NODE_LEAF | BAMR_NODE_HAS_CHILDREN)));
    caches_changed ();
    if (parent->flags & BAMR_NODE_IS_PRUNED)
      G(npruned)--;
    parent->flags = (parent->flags & ~BAMR_NODE_IS_PRUNED)
                    | BAMR_NODE_HAS_CHILDREN;
    BAMRNode * family = pool_alloc (pool_family);     /* zeroed */
    for (int c = 0; c < TREE_NCHILDREN; c++)
    {
      family [c].s = grid_scalar_table_new ();      /* no blocks: inactive */
      family [c].flags = BAMR_NODE_IS_PRUNED | BAMR_NODE_IS_INACTIVE;
      family [c].pid = -1;
    }
    for (int c = 0; c < TREE_NCHILDREN; c++)
    {
      Index q = index_child (up, c);
      INDEX_UNPACK (q);
      assert (!index_is_allocated (q));
      lidx_assign (& G (layers) [level], IJK, 1 << level, &family [c]);
    }
    G(nb) += TREE_NCHILDREN;
    G(npruned) += TREE_NCHILDREN;
    made += TREE_NCHILDREN;
  }
  BAMRNode * node = index_node (p);
  if (leaf && (node->flags & BAMR_NODE_IS_PRUNED))
  {
    assert (!(node->flags & BAMR_NODE_HAS_CHILDREN));
    caches_changed ();
    node->flags = (node->flags & ~BAMR_NODE_IS_PRUNED) | BAMR_NODE_LEAF;
    G(npruned)--;
    G(nlb)++;
    terminal_add (p.level, 1);
  }
  return made;
}

void tree_shift_resize (double x0, double y0, double z0, double h0)
{
  assert (h0 > 0.);
  G(x0) = x0, G(y0) = y0, G(z0) = z0, G (h0) = h0;
}

/* the lists again, by one walk: every node and the leaves in Morton order,
.. then each by level, still in that order */
static void caches_update (void)
{
  int n = G(nb), na = 0, nl = 0, no = 0, depth = 0;
  if (n > c_cap) {
    c_all = realloc (c_all, (size_t) n * sizeof (TreeItem));
    c_leaves = realloc (c_leaves, (size_t) n * sizeof (TreeItem));
    c_by_level = realloc (c_by_level, (size_t) n * sizeof (TreeItem));
    c_leaves_by_level = realloc (c_leaves_by_level, (size_t) n * sizeof (TreeItem));
    c_parents_by_level = realloc (c_parents_by_level, (size_t) n * sizeof (TreeItem));
    c_owned = realloc (c_owned, (size_t) n * sizeof (TreeItem));
    if (!c_all || !c_leaves || !c_by_level || !c_leaves_by_level || !c_parents_by_level
        || !c_owned)
      BAMR_ABORT ("out of memory");
    c_cap = n;
  }
  int count [BAMR_MAX_DEPTH + 1] = {0}, lcount [BAMR_MAX_DEPTH + 1] = {0};
  int pcount [BAMR_MAX_DEPTH + 1] = {0};
  morton (TREE_ROOT_INDEX) {
    TreeItem it = { (Index) {IJKL}, node };
    assert (na < n);
    c_all [na++] = it;
    count [level]++;
    if (node->flags & BAMR_NODE_LEAF) {
      c_leaves [nl++] = it;
      lcount [level]++;
      if (node->pid == G(rank))
        c_owned [no++] = it;
    }
    else if (tree_has_children (node))
      pcount [level]++;
    if (level > depth)
      depth = level;
  }
  caches.level_first [0] = caches.leaf_first [0] = caches.parent_first [0] = 0;
  for (int l = 0; l <= BAMR_MAX_DEPTH; l++) {
    caches.level_first [l + 1] = caches.level_first [l] + count [l];
    caches.leaf_first [l + 1] = caches.leaf_first [l] + lcount [l];
    caches.parent_first [l + 1] = caches.parent_first [l] + pcount [l];
  }
  int next [BAMR_MAX_DEPTH + 1], lnext [BAMR_MAX_DEPTH + 1], pnext [BAMR_MAX_DEPTH + 1];
  for (int l = 0; l <= BAMR_MAX_DEPTH; l++)
    next [l] = caches.level_first [l], lnext [l] = caches.leaf_first [l],
    pnext [l] = caches.parent_first [l];
  for (int q = 0; q < na; q++) {
    int l = c_all [q].p.level;
    c_by_level [next [l]++] = c_all [q];
    if (c_all [q].node->flags & BAMR_NODE_LEAF)
      c_leaves_by_level [lnext [l]++] = c_all [q];
    else if (tree_has_children (c_all [q].node))
      c_parents_by_level [pnext [l]++] = c_all [q];
  }
  caches.all = c_all, caches.nall = na;
  caches.leaves = c_leaves, caches.nleaves = nl;
  caches.owned = c_owned, caches.nowned = no;
  caches.by_level = c_by_level;
  caches.leaves_by_level = c_leaves_by_level;
  caches.parents_by_level = c_parents_by_level;
  caches.depth = depth;
  caches_are_dirty = 0;
}

const TreeCaches * tree_caches (void)
{
  if (caches_are_dirty)
    caches_update ();
  return &caches;
}

void tree_caches_changed (void)
{
  caches_changed ();
}

/*
.. The bits of x spread apart, dimension - 1 zeros between two: bit b to
.. bit dimension * b. Each step moves half of the groups left of the last
.. step, so the bits are placed in log2 (bits) steps, not one at a time
.. ("binary magic numbers"; libhmesh's spread () is the same). 2D: the
.. low 32 bits, 3D: the low 21, the most BAMR_MAX_DEPTH allows (the
.. checks at the top).
*/
#if dimension == 2
static inline uint64_t spread (uint64_t x)
{
  x &= UINT64_C (0xffffffff);
  x = (x | x << 16) & UINT64_C (0x0000ffff0000ffff);
  x = (x | x << 8)  & UINT64_C (0x00ff00ff00ff00ff);
  x = (x | x << 4)  & UINT64_C (0x0f0f0f0f0f0f0f0f);
  x = (x | x << 2)  & UINT64_C (0x3333333333333333);
  x = (x | x << 1)  & UINT64_C (0x5555555555555555);
  return x;
}
#elif dimension == 3
static inline uint64_t spread (uint64_t x)
{
  x &= UINT64_C (0x1fffff);
  x = (x | x << 32) & UINT64_C (0x001f00000000ffff);
  x = (x | x << 16) & UINT64_C (0x001f0000ff0000ff);
  x = (x | x << 8)  & UINT64_C (0x100f00f00f00f00f);
  x = (x | x << 4)  & UINT64_C (0x10c30c30c30c30c3);
  x = (x | x << 2)  & UINT64_C (0x1249249249249249);
  return x;
}
#endif

/*
.. The position of p on the Morton curve: the bits of its lower corner at
.. the finest level, interleaved (x lowest). A node and its first leaf
.. share the key, and a subtree's keys are one interval, so the key of
.. any index says which process's piece holds it (tree_partition_owner
.. ()).
*/
uint64_t tree_morton_key (Index p)
{
  assert (p.level >= 0 && p.level <= BAMR_MAX_DEPTH);
  int s = BAMR_MAX_DEPTH - p.level;
  #if dimension == 1
  return (uint64_t) p.i << s;         /* nothing to interleave */
  #elif dimension == 2
  return spread ((uint64_t) p.i << s) | spread ((uint64_t) p.j << s) << 1;
  #else
  return spread ((uint64_t) p.i << s) | spread ((uint64_t) p.j << s) << 1
         | spread ((uint64_t) p.k << s) << 2;
  #endif
}

void tree_changed_under (const char * loop)
{
  fprintf (stderr, "%s (): the tree changed during the loop: refine or "
           "coarsen after it, not in it\n", loop);
  abort ();
}

double tree_get_origin (int axis)
{
  assert (axis >= 0 && axis < 3);
  return axis == 0 ? G(x0) : axis == 1 ? G(y0) : G(z0);
}

double tree_get_size (void)
{
  return G(h0);
}

int tree_get_depth (void)
{
  int l = bits_highest_u32 (G(level_bits));     /* -1: no tree */
  return l < 0 ? 0 : l;
}

int tree_get_depth_max (void)
{
  return G(depthmax);
}

int tree_get_depth_min (void)
{
  int l = bits_lowest_u32 (G(level_bits));      /* -1: no tree */
  return l < 0 ? 0 : l;
}

int tree_get_nblocks (void)
{
  return G(nb);
}

int tree_get_block_m (void)
{
  return bm;
}

int tree_get_nleaves (void)
{
  return G(nlb);
}
