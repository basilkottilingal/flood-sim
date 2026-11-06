#include <assert.h>
#include <stdlib.h>

#include "block-ghosts.h"

/* the index of cell (x, y, z) in a stored block of m^dimension cells */
#if dimension == 1
# define STORED(m, x, y, z)  (x)
#elif dimension == 2
# define STORED(m, x, y, z)  ((x) + (m) * (y))
#else
# define STORED(m, x, y, z)  ((x) + (m) * ((y) + (m) * (z)))
#endif

/* log2 of m, a power of 2: its trailing zero bits */
static inline int log2i (unsigned m)
{
  assert (m && !(m & (m - 1)));
#if defined (__GNUC__) || defined (__clang__)
  return __builtin_ctz (m);
#else
  /* de Bruijn lookup: m has a single bit set */
  static const unsigned char tbl [32] =
  {
     0,  1, 28,  2, 29, 14, 24,  3, 30, 22, 20, 15, 25, 17,  4,  8,
    31, 27, 13, 23, 21, 19, 16,  7, 26, 12, 18,  6, 11,  5, 10,  9
  };
  return tbl [(unsigned) (m * 0x077CB531u) >> 27];
#endif
}

/* a coordinate outside [0, n), reflected back inside */
static inline int reflect (int x, int n)
{
  return x < 0 ? -1 - x : x >= n ? 2 * n - 1 - x : x;
}

/*
.. The value of scalar s in the cell (x, y, z) of level l, in cells of that
.. level across the whole domain: from the node of level l that holds it, or
.. the nearest coarser one.
*/
static Real sample (int s, int l, int x, int y, int z, int m, int bd)
{
  for (int L = l; L >= 0; L--)
  {
    int sh = l - L, X = x >> sh, Y = y >> sh, Z = z >> sh;
    Index b = { X >> bd,
  #if dimension > 1
                Y >> bd,
  #endif
  #if dimension > 2
                Z >> bd,
  #endif
                L };
    (void) Y, (void) Z;
    if (index_is_allocated (b))
    {
      const Real * blk = block_data (index_node (b)->s, s);
      return blk [STORED (m, X & (m - 1), Y & (m - 1), Z & (m - 1))];
    }
  }
  assert (0);                   /* the root covers the domain */
  return 0;
}

void block_gather (Real * w, BlockShape shape, Index p, int s)
{
  int m = tree_block_m (), bd = log2i (m), l = p.level, n = m << l;
  assert (shape.m == m && shape.g <= m && index_is_allocated (p));
  const Real * own = block_data (index_node (p)->s, s);

  /* the block's own cells */
  block_interior (shape)
    w [c] = own [o];

  /* the ghosts: what is not interior */
  int x0 = p.i * m, y0 = 0, z0 = 0;
  #if dimension > 1
  y0 = p.j * m;
  #endif
  #if dimension > 2
  z0 = p.k * m;
  #endif
  (void) y0, (void) z0;
  block_cells (shape, shape.g)
  {
    int x = ci - shape.ox, y = 0, z = 0;
  #if dimension > 1
    y = cj - shape.g;
  #endif
  #if dimension > 2
    z = ck - shape.g;
  #endif
    if (x >= 0 && x < m && y >= 0 && y < m && z >= 0 && z < m)
      continue;
    w [c] = sample (s, l, reflect (x0 + x, n),
  #if dimension > 1
                    reflect (y0 + y, n),
  #else
                    0,
  #endif
  #if dimension > 2
                    reflect (z0 + z, n),
  #else
                    0,
  #endif
                    m, bd);
  }
}

void block_scatter (Index p, int s, const Real * w, BlockShape shape)
{
  assert (shape.m == tree_block_m () && index_is_allocated (p));
  Real * own = block_data (index_node (p)->s, s);
  block_interior (shape)
    own [o] = w [c];
}

/*
.. morton () visits children before their parent, so each parent is
.. averaged from children already up to date: one pass from the leaves up.
*/
void tree_restrict (int s)
{
  int m = tree_block_m (), h = m / 2;
  morton (TREE_ROOT_INDEX)
  {
    if (!tree_has_children (node))
      continue;
    Real * parent = block_data (sblocks, s);
    Index self = {IJKL};
    const Real * child [TREE_NCHILDREN];
    for (int c = 0; c < TREE_NCHILDREN; c++)
      child [c] = block_data (index_node (index_child (self, c))->s, s);
    for (int x = 0; x < m; x++)
  #if dimension > 1
    for (int y = 0; y < m; y++)
  #endif
  #if dimension > 2
    for (int z = 0; z < m; z++)
  #endif
    {
  #if dimension < 2
      int y = 0;
  #endif
  #if dimension < 3
      int z = 0;
  #endif
      /* the child holding this cell, and the first of its 2^d fine cells */
      int c = (x >= h) | (dimension > 1 && y >= h) << 1
              | (dimension > 2 && z >= h) << 2;
      int fx = 2 * (x & (h - 1)), fy = 2 * (y & (h - 1)), fz = 2 * (z & (h - 1));
      (void) fy, (void) fz;
      RealSum sum = 0;
      for (int a = 0; a < 2; a++)
  #if dimension > 1
      for (int b = 0; b < 2; b++)
  #endif
  #if dimension > 2
      for (int e = 0; e < 2; e++)
  #endif
      {
  #if dimension < 2
        int b = 0;
  #endif
  #if dimension < 3
        int e = 0;
  #endif
        (void) b, (void) e;
        sum += child [c] [STORED (m, fx + a, fy + b, fz + e)];
      }
      parent [STORED (m, x, y, z)] = (Real) (sum / TREE_NCHILDREN);
    }
  }
}
