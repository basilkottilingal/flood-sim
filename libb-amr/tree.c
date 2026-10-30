#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <assert.h>

#include "tree.h"
#include "liballoc/allocator.h"

#define G(x) grid_info_##x

static double G(x0) = 0.;     /* left corner of domain            */
static double G(y0) = 0.;     /* bottom corner of domain          */
static double G(z0) = 0.;
static double G(h0) = 1.;     /* domain [x0, x0+h] x [y0, y0+h]   */
static int G(nb) = 1;         /* number of blocks                 */
static int G(nlb) = 1;        /* number of leaf blocks            */
static int G(depth) = 1;      /* depth of the current bAMR        */
static int G(depthmin) = 1;   /* min {level of leaves}            */
static int G(depthmax) = 1;   /* max depth forced by user         */
static int G(blockdepth) = 1; /* number of cells in block 2^(2*blockdepth) */
struct _Memindex G(layers) [BAMR_MAX_DEPTH + 1] = {0};
static Allocator * alloc_node  = NULL;
static Allocator * alloc_block = NULL;

BAMRNode * allocate_node ()
{
  BAMRNode * node = allocate (alloc_node, sizeof (BAMRNode));
  // fixme : 1) allocator. 2)append scalars 3) allocate block
  return node;
}

#if dimension == 2
#define is_node()   tree_is_allocated (i,j,level)
#define node_print(msg)                        \
  printf("%s : [%d %d %d]", msg, i,j,level);   \
  fflush (stdout);
#endif

#define is_leaf()   (node->flags & BAMR_NODE_LEAF)
#define is_pruned() (node->flags & BAMR_NODE_IS_PRUNED)
#define is_parent() (node->flags & BAMR_NODE_HAS_CHILDREN)
#define children() do {                        \
    int ib = i, jb = j, lb = level++;          \
    for (int i = 2*ib; i < 2*ib + 2; ++i)      \
      for (int j = 2*jb; j < 2*jb + 2; ++j)    \
      {                                        \
        BAMR_NODE_ENV ();  
#define end()                                  \
      }                                        \
    i = ib, j = jb, level = lb;                \
  } while (0)
#define neighborhood(N) do {                   \
    int ib = i; jb = j, lb = level;            \
    for (int i = ib-N; i <= ib + N; ++i)       \
      for (int j = jb-N; j <= jb + N; ++j)     \
      {                                        \
        BAMR_NODE_ENV ();
#define dummy                                  \
      }                                        \
  }
#undef dummy
#define neighbors(N) do {                      \
    int ib = i, jb = j, lb = level;            \
    level++;                                   \
    for (int i = ib-N; i <= ib + N; ++i)       \
      for (int j = jb-N; j <= jb + N; ++j)     \
        if (i != ib  && j != jb)               \
        {                                      \
          BAMR_NODE_ENV ();
#define dummy                                  \
        }                                      \
  }
#undef dummy

static int add_children (int IJKLargs)
{
  assert (level >= G(depthmin));
  if (level == G(depthmax))
    return 0;
  BAMR_NODE_ENV ();
  assert (node && !is_parent ());
  node->flags |= BAMR_NODE_HAS_CHILDREN;
  if (level == G (depth))
    G (depth)++;
  Memindex m = & G (layers) [level + 1];
  int len = 1 << (level+1);
  children ()
  {
    assert (!node);
    node = allocate_node ();
    mem_assign (m, i, j, len, node);
  } end ();
  //fixme : depthmin is not updated
  return 1;    
}

static int divide_leaf (int IJKLargs)
{
  if (!add_children (IJKL))
    return 0;
  BAMR_NODE_ENV ();
  assert (is_leaf ());
  node->flags &= ~BAMR_NODE_LEAF;
  children ()
  {
    node->flags |= BAMR_NODE_LEAF;
  } end ();
  return 1;
}

static int grid_init (int depthmax, int blockdepth)
{
  G (depthmax)   = depthmax;
  G (depth)      = 0;
  G (depthmin)   = 0;
  G (blockdepth) = blockdepth;
  G (nlb)        = 1;
  G (nb)         = 1;

  int i = 0, j = 0, level = 0;
  BAMR_NODE_ENV ();
  if (node) {
    exit (-1);
    //fixme : rather warning : reinitializing grid;
  }
  Memindex m = & G (layers) [level];
  m->len = 1 << level; // fixme : add ghost layer
  alloc_node = allocator ();
  alloc_block = allocator (); //fixme : use a large allocator
  node = allocate_node ();
  node->flags |= BAMR_NODE_LEAF;
  mem_assign (m, i, j, m->len, node);
  divide_leaf (i, j, level);
  //fixme 
  return 0;
}

int tree (int depthmax, int blockdepth)
{
  depthmax = depthmax > BAMR_MAX_DEPTH ? BAMR_MAX_DEPTH :
    depthmax <= 0 ? 1 : depthmax;
  blockdepth = blockdepth < 4 ? 4 : blockdepth > 10 ? 10 : blockdepth;
  grid_init (depthmax, blockdepth);
  return 1;
}

void tree_shift_resize (double x0, double y0, double z0, double h0)
{
  assert (h0 > 0.f);
  G(x0) = x0, G(y0) = y0, G(z0) = z0, G (h0) = h0;
}
