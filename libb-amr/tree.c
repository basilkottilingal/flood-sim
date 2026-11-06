#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

#include "tree.h"
#include "grid-scalar.h"
#include "block-ghosts.h"
#include "pool.h"

#define G(x) grid_info_##x

static double G(x0) = 0.;     /* left corner of domain            */
static double G(y0) = 0.;     /* bottom corner of domain          */
static double G(z0) = 0.;
static double G(h0) = 1.;     /* domain [x0, x0+h] x [y0, y0+h]   */
static int G(nb) = 0;         /* number of blocks                 */
static int G(nlb) = 0;        /* number of leaf blocks            */
static int G(depth) = 0;      /* depth of the current bAMR        */
static int G(depthmin) = 0;   /* min {level of leaves}            */
static int G(depthmax) = 1;   /* max depth forced by user         */
static int G(blockdepth) = 1; /* 2^blockdepth cells per block side */
static int G(nleaves) [BAMR_MAX_DEPTH + 1];   /* leaves per level */
struct _Memindex G(layers) [BAMR_MAX_DEPTH + 1] = {{0}};
static Pool * pool_node = NULL;       /* the root, which has no siblings */
static Pool * pool_family = NULL;     /* TREE_NCHILDREN nodes at once    */

/*
.. cells a block stores: m^dimension, m = 2^blockdepth, and no ghost cells.
.. Code that needs ghosts fills a scratch block of its own from the
.. neighbours: block_gather () in block-ghosts.h.
*/
#define bm           (1 << G (blockdepth))
#if dimension == 1
# define block_size()     (bm)
#elif dimension == 2
# define block_size()     (bm * bm)
#else
# define block_size()     (bm * bm * bm)
#endif

#define is_leaf()   (node->flags & BAMR_NODE_LEAF)
#define is_pruned() (node->flags & BAMR_NODE_IS_PRUNED)
#define is_parent() (node->flags & BAMR_NODE_HAS_CHILDREN)

static void update_depthmin (void)
{
  int l = 0;
  while (l < G(depth) && !G(nleaves) [l])
    l++;
  G(depthmin) = l;
}

void tree_destroy (void)
{
  for (int l = 0; l <= BAMR_MAX_DEPTH; ++l)
  {
    Memindex m = & G(layers) [l];
    if (m->b)
      mem_destroy (m, m->len);
    *m = (struct _Memindex) {0};
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
  G(depth) = 0;
  G(depthmin) = 0;
  G(depthmax) = 1;
  G(blockdepth) = 1;
}

/* the root: a node and blocks of its own */
static BAMRNode * allocate_root (void)
{
  BAMRNode * node = pool_alloc (pool_node);   /* zeroed: flags 0 */
  node->s = grid_scalar_table_new ();
  grid_scalar_root_blocks (node->s);
  return node;
}

/* recount the deepest and coarsest levels of leaves */
static void update_depth (void)
{
  while (G(depth) > 0 && !G(nleaves) [G(depth)])
    G(depth)--;
  update_depthmin ();
}

static int add_children (int IJKLargs)
{
  assert (level >= G(depthmin));
  if (level >= G(depthmax))
    return 0;
  BAMR_NODE_ENV ();
  assert (node && !is_parent ());
  node->flags |= BAMR_NODE_HAS_CHILDREN;
  if (level == G (depth))
    G (depth)++;
  Memindex m = & G (layers) [level + 1];
  int len = 1 << (level + 1);
  /*
  .. The family at once: its nodes in one object, child c the c-th, and for
  .. each scalar its blocks in one object. The blocks are not zeroed:
  .. divide_leaf () fills every cell.
  */
  BAMRNode * family = pool_alloc (pool_family);     /* zeroed: flags 0 */
  void *** t [TREE_NCHILDREN];
  for (int c = 0; c < TREE_NCHILDREN; c++)
    t [c] = family [c].s = grid_scalar_table_new ();
  grid_scalar_family_new (t, 1);
  int c = 0;
  children_of ((Index) {IJKL})
  {
    assert (!node);
    mem_assign (m, IJK, len, &family [c++]);
  }
  G (nb) += TREE_NCHILDREN;
  return 1;
}

static int divide_leaf (int IJKLargs)
{
  if (!add_children (IJKL))
    return 0;
  BAMR_NODE_ENV ();
  assert (is_leaf ());
  node->flags &= ~BAMR_NODE_LEAF;
  children_of ((Index) {IJKL})
    node->flags |= BAMR_NODE_LEAF;
  G (nleaves) [level]--;
  G (nleaves) [level + 1] += TREE_NCHILDREN;
  G (nlb) += TREE_NCHILDREN - 1;
  update_depthmin ();
  /* the children's values: linear prolongation, which keeps the mean */
  foreach_scalar (s)
    block_prolong ((Index) {IJKL}, s);
  return 1;
}

/* the children of the leaf-parent at IJKL, gone: it is a leaf again */
static void remove_children (int IJKLargs)
{
  BAMR_NODE_ENV ();
  assert (node && is_parent ());
  BAMRNode * parent = node, * family = NULL;
  void *** t [TREE_NCHILDREN];
  int c = 0;
  children_of ((Index) {IJKL})
  {
    assert (node && (node->flags & BAMR_NODE_LEAF));
    if (!c)
      family = node;                  /* child 0 starts the family's object */
    t [c++] = node->s;
    mem_free (& G (layers) [level], IJK, 1 << level, node);
  }
  grid_scalar_family_free (t);
  pool_free (pool_family, family);
  parent->flags = (parent->flags & ~BAMR_NODE_HAS_CHILDREN) | BAMR_NODE_LEAF;
  G (nb) -= TREE_NCHILDREN;
  G (nlb) -= TREE_NCHILDREN - 1;
  G (nleaves) [level + 1] -= TREE_NCHILDREN;
  G (nleaves) [level]++;
  update_depth ();
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

/* the node covering q, q or one of its ancestors: the root covers all */
static Index covering (Index q)
{
  while (!index_is_allocated (q))
  {
    q.level--;
    q.i >>= 1;
  #if dimension > 1
    q.j >>= 1;
  #endif
  #if dimension > 2
    q.k >>= 1;
  #endif
  }
  return q;
}

/*
.. Refining p, at level l, puts leaves of level l + 1 beside whatever touches
.. p. Each index of level l around p must therefore hold a node already, or
.. the coarser leaf covering it would touch leaves two levels finer: that
.. leaf is refined first, which may in turn refine others, coarser still.
*/
int tree_refine (Index p)
{
  if (!index_is_allocated (p) || p.level >= G(depthmax)
      || !(index_node (p)->flags & BAMR_NODE_LEAF))
    return 0;
  int count = 0;
  neighbors_of (p, 1)
  {
    Index q = {IJKL};
    while (inside (q) && !index_is_allocated (q))
    {
      int n = tree_refine (covering (q));
      assert (n);                     /* coarser than p: below depth max */
      count += n;
    }
  }
  INDEX_UNPACK (p);
  count += divide_leaf (IJKL);
  return count;
}

int tree_coarsen (Index p)
{
  if (!index_is_allocated (p) || !tree_has_children (index_node (p)))
    return 0;
  children_of (p)
    if (tree_has_children (node))
      return 0;                       /* a child is not a leaf */
  /* a neighbour's grandchildren would touch p two levels finer */
  neighbors_of (p, 1)
    if (node && tree_has_children (node))
      children_of ((Index) {IJKL})
        if (tree_has_children (node))
          return 0;
  /* p gets the mean of its children, then they go */
  foreach_scalar (s)
    block_restrict (p, s);
  INDEX_UNPACK (p);
  remove_children (IJKL);
  return 1;
}

static int grid_init (int depthmax, int blockdepth)
{
  if (pool_node)
  {
    fprintf (stderr, "tree: cannot reinit grid; call tree_destroy () first\n");
    exit (-1);
  }
  G (depthmax)   = depthmax;
  G (depth)      = 0;
  G (depthmin)   = 0;
  G (blockdepth) = blockdepth;

  pool_node = pool_new (sizeof (BAMRNode));
  pool_family = pool_new (TREE_NCHILDREN * sizeof (BAMRNode));
  grid_scalar_init (block_size ());

  Index root = {TREE_ROOT};
  INDEX_UNPACK (root);
  Memindex m = & G (layers) [level];
  BAMRNode * node = allocate_root ();
  node->flags |= BAMR_NODE_LEAF;
  mem_assign (m, IJK, 1 << level, node);
  G (nb)  = 1;
  G (nlb) = 1;
  G (nleaves) [level] = 1;
  divide_leaf (IJKL);
  return 0;
}

int tree (int depthmax, int blockdepth)
{
  depthmax = depthmax > BAMR_MAX_DEPTH ? BAMR_MAX_DEPTH :
    depthmax <= 0 ? 1 : depthmax;
  /*
  .. A block has (2^blockdepth)^dimension Reals, and a scratch block with
  .. ghosts (2^blockdepth + 2 GHOSTS)^dimension, so the largest blockdepth
  .. shrinks with the dimension: 3D at 6 is already 64^3 doubles, 2 MiB a
  .. block. fixme : tune better
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

void tree_shift_resize (double x0, double y0, double z0, double h0)
{
  assert (h0 > 0.);
  G(x0) = x0, G(y0) = y0, G(z0) = z0, G (h0) = h0;
}

int tree_depth (void)
{
  return G(depth);
}

int tree_depth_max (void)
{
  return G(depthmax);
}

int tree_depth_min (void)
{
  return G(depthmin);
}

int tree_nblocks (void)
{
  return G(nb);
}

int tree_block_m (void)
{
  return bm;
}

int tree_nleaves (void)
{
  return G(nlb);
}
