#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

#include "tree.h"
#include "grid-scalar.h"
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
static Pool * pool_node = NULL;

/* cells in a block, ghosts included: bb^dimension */
#define bb           ((1 << G (blockdepth)) + 2 * GHOSTS)
#if dimension == 1
# define blockwidth()     (bb)
#elif dimension == 2
# define blockwidth()     (bb * bb)
#else
# define blockwidth()     (bb * bb * bb)
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
  pool_destroy (pool_node);
  pool_node = NULL;
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

static BAMRNode * allocate_node (void)
{
  BAMRNode * node = pool_alloc (pool_node);   /* zeroed: flags 0 */
  node->s = grid_scalar_node_new ();
  return node;
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
  children_of ((Index) {IJKL})
  {
    assert (!node);
    node = allocate_node ();
    mem_assign (m, IJK, len, node);
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
  grid_scalar_init (blockwidth ());

  Index root = {TREE_ROOT};
  INDEX_UNPACK (root);
  Memindex m = & G (layers) [level];
  BAMRNode * node = allocate_node ();
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
  .. A block has (2^blockdepth + 2 GHOSTS)^dimension doubles, so the largest
  .. blockdepth shrinks with the dimension: 3D at 6 is already 68^3 doubles,
  .. 2.4 MB a block. fixme : tune better
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

int tree_nblocks (void)
{
  return G(nb);
}

int tree_nleaves (void)
{
  return G(nlb);
}
