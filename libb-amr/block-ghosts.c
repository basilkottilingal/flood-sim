#include <assert.h>
#include <stdlib.h>
#include <string.h>

#include "block-ghosts.h"

/* the index of cell (x, y, z) in a stored block of m^dimension cells */
#if dimension == 1
# define STORED(m, x, y, z)  (x)
#elif dimension == 2
# define STORED(m, x, y, z)  ((x) + (m) * (y))
#else
# define STORED(m, x, y, z)  ((x) + (m) * ((y) + (m) * (z)))
#endif

/*
.. Where the cell x of a ghost region comes from, along one axis: x is in
.. the block's own cells (0 to m - 1 inside), d the direction of the region
.. (-1, 0 or 1), mirror whether that side is outside the domain, and off,
.. sh where the source block sits: the cell x' of the block in direction d
.. (or of the block itself, reflected, if mirror) is the cell (off + x') >> sh
.. of the source block, sh levels coarser.
*/
static inline int source_cell (int x, int d, int mirror, int m, int off, int sh)
{
  int xl = mirror ? (d < 0 ? -1 - x : 2 * m - 1 - x) : x - d * m;
  return (off + xl) >> sh;
}

/*
.. The ghosts come by region: one for each of the 3^dimension - 1 directions
.. around the block, its faces, edges and corners. For each:
..
..   - the sides outside the domain are reflected: that component of the
..     direction becomes 0, and the cells are mirrored, so a corner region
..     across a wall reads the neighbour along the wall, not the block itself;
..   - the source is the node in the direction left, at the same level (a
..     leaf, or a parent holding the mean of its children: tree_restrict ())
..     or else the nearest coarser node covering it: found once per region,
..     not once per cell;
..   - each row of the region is then copied, with memcpy where it is a run
..     of the same level not reflected along i.
*/
void block_gather (Real * w, BlockShape shape, Index p, int s)
{
  int m = tree_get_block_m (), l = p.level, g = shape.g, nb = 1 << l;
  assert (shape.m == m && g <= m && index_is_allocated (p));
  const Real * own = block_data (index_node (p)->s, s);

  /* the block's own cells */
  block_interior (shape)
    w [c] = own [o];
  if (!g)
    return;

  int P [3] = { p.i, 0, 0 }, ndir = 1;
  #if dimension > 1
  P [1] = p.j;
  #endif
  #if dimension > 2
  P [2] = p.k;
  #endif
  for (int a = 0; a < dimension; a++)
    ndir *= 3;

  for (int dir = 0; dir < ndir; dir++)
  {
    int d [3] = { 0, 0, 0 }, mirror [3] = { 0, 0, 0 }, q [3] = { 0, 0, 0 };
    int lo [3] = { 0, 0, 0 }, hi [3] = { 1, 1, 1 }, t = dir, centre = 1;
    for (int a = 0; a < dimension; a++)
    {
      d [a] = t % 3 - 1;
      t /= 3;
      centre &= d [a] == 0;
      q [a] = P [a] + d [a];
      if (q [a] < 0 || q [a] >= nb)       /* outside: reflect */
      {
        mirror [a] = 1;
        q [a] = P [a];
      }
      lo [a] = d [a] < 0 ? -g : d [a] == 0 ? 0 : m;
      hi [a] = d [a] < 0 ? 0 : d [a] == 0 ? m : m + g;
    }
    if (centre)
      continue;

    /* the source: the node at q, or the nearest coarser one covering it */
    const Real * src = NULL;
    int sh = 0, off [3] = { 0, 0, 0 };
    for (int L = l; L >= 0 && !src; L--)
    {
      sh = l - L;
      Index b = { q [0] >> sh,
  #if dimension > 1
                  q [1] >> sh,
  #endif
  #if dimension > 2
                  q [2] >> sh,
  #endif
                  L };
      if (index_is_allocated (b))
      {
        src = block_data (index_node (b)->s, s);
        /* where q's cells start in it, in cells of level l */
        for (int a = 0; a < dimension; a++)
          off [a] = q [a] * m - (((q [a] >> sh) * m) << sh);
      }
    }
    assert (src);                       /* the root covers the domain */

    int run = !mirror [0] && sh == 0;   /* rows contiguous in both */
    for (int z = lo [2]; z < hi [2]; z++)
    for (int y = lo [1]; y < hi [1]; y++)
    {
      int sy = source_cell (y, d [1], mirror [1], m, off [1], sh);
      int sz = source_cell (z, d [2], mirror [2], m, off [2], sh);
      (void) sy, (void) sz;
      const Real * srow = src + STORED (m, 0, sy, sz);
      Real * drow = w + shape.ox + (y + g) * shape.sy + (z + g) * shape.sz;
      if (run)
        memcpy (drow + lo [0], srow + lo [0] - d [0] * m,
                (size_t) (hi [0] - lo [0]) * sizeof (Real));
      else
        for (int x = lo [0]; x < hi [0]; x++)
          drow [x] = srow [source_cell (x, d [0], mirror [0], m, off [0], sh)];
    }
  }
}

void block_scatter (Index p, int s, const Real * w, BlockShape shape)
{
  assert (shape.m == tree_get_block_m () && index_is_allocated (p));
  Real * own = block_data (index_node (p)->s, s);
  block_interior (shape)
    own [o] = w [c];
}

/* the scratch index of cell x [] of a block, interior-relative */
static inline long scratch_index (BlockShape b, const int x [3])
{
  return (x [0] + b.ox) + (long) (x [1] + b.g) * b.sy
         + (long) (x [2] + b.g) * b.sz;
}

/* the stride of axis a in a scratch block */
static inline long scratch_stride (BlockShape b, int a)
{
  return a == 0 ? 1 : a == 1 ? b.sy : b.sz;
}

/*
.. Linear prolongation, which keeps the mean: the fine cell at offset
.. sgn = -1 or +1 from the centre of its coarse cell C, along each axis, is
..
..   u (C) + sum over axes of sgn / 8 * (u (C + 1) - u (C - 1))
..
.. (a quarter of a coarse cell times the central gradient). The 2^d fine
.. cells of C have offsets of both signs on every axis, so their mean is
.. u (C). The coarse values come with one ghost layer, from block_gather ().
*/
void block_prolong (Index p, int s)
{
  int m = tree_get_block_m ();
  assert (index_is_allocated (p) && tree_has_children (index_node (p)));
  BlockShape sh = tree_scratch_shape (1);
  Real * u = block_scratch_new (sh);
  block_gather (u, sh, p, s);
  int c = 0;
  children_of (p)
  {
    Real * fine = block_data (sblocks, s);
    int cb [3] = { c & 1, (c >> 1) & 1, (c >> 2) & 1 };
    int hi [3] = { m, dimension > 1 ? m : 1, dimension > 2 ? m : 1 };
    int x [3];
    for (x [2] = 0; x [2] < hi [2]; x [2]++)
    for (x [1] = 0; x [1] < hi [1]; x [1]++)
    for (x [0] = 0; x [0] < hi [0]; x [0]++)
    {
      int C [3] = { 0, 0, 0 }, sgn [3] = { 0, 0, 0 };
      for (int a = 0; a < dimension; a++)
      {
        int f = cb [a] * m + x [a];         /* in the parent's fine cells */
        C [a] = f >> 1;
        sgn [a] = f & 1 ? 1 : -1;
      }
      long k = scratch_index (sh, C);
      RealSum v = u [k];
      for (int a = 0; a < dimension; a++)
      {
        long st = scratch_stride (sh, a);
        v += sgn [a] * (u [k + st] - u [k - st]) / 8.;
      }
      fine [STORED (m, x [0], x [1], x [2])] = (Real) v;
    }
    c++;
  }
  free (u);
}

/* the parent p: the mean of its children's cells */
void block_restrict (Index p, int s)
{
  int m = tree_get_block_m (), h = m / 2;
  assert (index_is_allocated (p) && tree_has_children (index_node (p)));
  Real * parent = block_data (index_node (p)->s, s);
  const Real * child [TREE_NCHILDREN];
  for (int c = 0; c < TREE_NCHILDREN; c++)
    child [c] = block_data (index_node (index_child (p, c))->s, s);
  int hi [3] = { m, dimension > 1 ? m : 1, dimension > 2 ? m : 1 };
  int x [3];
  for (x [2] = 0; x [2] < hi [2]; x [2]++)
  for (x [1] = 0; x [1] < hi [1]; x [1]++)
  for (x [0] = 0; x [0] < hi [0]; x [0]++)
  {
    /* the child holding this cell, and the first of its 2^d fine cells */
    int c = 0, f [3] = { 0, 0, 0 };
    for (int a = 0; a < dimension; a++)
    {
      c |= (x [a] >= h) << a;
      f [a] = 2 * (x [a] & (h - 1));
    }
    RealSum sum = 0;
    for (int k = 0; k < TREE_NCHILDREN; k++)
      sum += child [c] [STORED (m, f [0] + (k & 1), f [1] + ((k >> 1) & 1),
                                f [2] + ((k >> 2) & 1))];
    parent [STORED (m, x [0], x [1], x [2])] = (Real) (sum / TREE_NCHILDREN);
  }
}

/*
.. morton () visits children before their parent, so each parent is
.. averaged from children already up to date: one pass from the leaves up.
*/
void tree_restrict (int s)
{
  morton (TREE_ROOT_INDEX)
    if (tree_has_children (node))
      block_restrict ((Index) {IJKL}, s);
}
