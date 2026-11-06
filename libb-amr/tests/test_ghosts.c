/* Ghost cells for blocks stored without them: restriction, gathering at
 * level jumps and at the boundary, scattering. Includes tree.c to refine. */
#include "tree.c"
#include "block-ghosts.h"
#include <math.h>
static int checks, fails;
#define CHECK(c) do { checks++; if (!(c)) { fails++; if (fails < 10) printf ("FAIL line %d: %s\n", __LINE__, #c); } } while (0)
scalar u, v;

/* f at the centre of cell (x, y, z) of level l, in cells of the domain */
static double f (int l, int x, int y, int z)
{
  double n = (double) (tree_block_m () << l);
  return 1 + 2 * (x + 0.5) / n + 3 * (dimension > 1) * (y + 0.5) / n
           + 5 * (dimension > 2) * (z + 0.5) / n;
}
static int refl (int x, int n) { return x < 0 ? -1 - x : x >= n ? 2 * n - 1 - x : x; }

/* the expected value of a ghost: f at the finest existing cell covering it */
static double expected (int l, int x, int y, int z)
{
  int m = tree_block_m (), n = m << l;
  x = refl (x, n); y = refl (y, n); z = refl (z, n);
  for (int L = l; L >= 0; L--) {
    int X = x >> (l - L), Y = y >> (l - L), Z = z >> (l - L);
    Index b = { X / m,
#if dimension > 1
                Y / m,
#endif
#if dimension > 2
                Z / m,
#endif
                L };
    if (index_is_allocated (b))
      return f (L, X, Y, Z);
  }
  return NAN;
}

int main (void)
{
  int m;
  tree (4, 4);
  m = tree_block_m ();
  CHECK (m == 16 && tree_block_shape ().cells == (long) pow (16, dimension));
  /* refine: corner chain to level 3, and one more level-1 leaf to level 2 */
  Index a1 = index_child (TREE_ROOT_INDEX, 0), a2 = index_child (a1, TREE_NCHILDREN - 1);
  { INDEX_UNPACK (a1); CHECK (divide_leaf (IJKL)); }
  { INDEX_UNPACK (a2); CHECK (divide_leaf (IJKL)); }
  Index b1 = index_child (TREE_ROOT_INDEX, TREE_NCHILDREN - 1);
  { INDEX_UNPACK (b1); CHECK (divide_leaf (IJKL)); }

  /* every leaf: f at its cell centres */
  BlockShape store = tree_block_shape ();
  foreach_leaf () {
    Real * b = block_data (sblocks, u);
    block_cells (store, 0) {
      int x = i * m + ci, y = 0, z = 0;
#if dimension > 1
      y = j * m + cj;
#endif
#if dimension > 2
      z = k * m + ck;
#endif
      b [c] = f (level, x, y, z);
    }
  }
  /* parents: the mean of their children, which is f at their centres */
  tree_restrict (u);
  double err = 0;
  morton (TREE_ROOT_INDEX) {
    Real * b = block_data (sblocks, u);
    block_cells (store, 0) {
      int x = i * m + ci, y = 0, z = 0;
#if dimension > 1
      y = j * m + cj;
#endif
#if dimension > 2
      z = k * m + ck;
#endif
      err = fmax (err, fabs (b [c] - f (level, x, y, z)));
    }
  }
  CHECK (err < 1e-5);

  /* every leaf, gathered with 2 ghosts: interior is its own, ghosts as expected */
  for (int g = 1; g <= 2; g++) {
    BlockShape shape = tree_scratch_shape (g);
    Real * w = malloc (shape.cells * sizeof (Real));
    double gerr = 0;
    long ghosts = 0;
    foreach_leaf () {
      block_gather (w, shape, (Index) {IJKL}, u);
      block_cells (shape, g) {
        int x = ci - shape.ox, y = 0, z = 0;
#if dimension > 1
        y = cj - shape.g;
#endif
#if dimension > 2
        z = ck - shape.g;
#endif
        int gx = i * m + x, gy = 0, gz = 0;
#if dimension > 1
        gy = j * m + y;
#endif
#if dimension > 2
        gz = k * m + z;
#endif
        ghosts += !(x >= 0 && x < m && y >= 0 && y < m && z >= 0 && z < m);
        gerr = fmax (gerr, fabs (w [c] - expected (level, gx, gy, gz)));
      }
    }
    CHECK (gerr < 1e-5);
    CHECK (ghosts == tree_nleaves () * ((long) pow (m + 2 * g, dimension) - (long) pow (m, dimension)));
    free (w);
  }

  /* a stencil through foreach_leaf_gathered, against the same through
   * block_gather by hand; then scatter round trip */
  {
    BlockShape sh = tree_scratch_shape (1);
    Real * w2 = malloc (sh.cells * sizeof (Real));
    double diff = 0;
    foreach_leaf_gathered (u, 1) {
      Real * out = block_data (sblocks, v);
      block_interior (shape)
        out [o] = cell_at (w, 1
#if dimension > 1
                           , 0
#endif
#if dimension > 2
                           , 0
#endif
                           ) - cell_at (w, -1
#if dimension > 1
                           , 0
#endif
#if dimension > 2
                           , 0
#endif
                           );
    }
    foreach_leaf () {
      block_gather (w2, sh, (Index) {IJKL}, u);
      Real * out = block_data (sblocks, v);
      block_interior (sh) {
        Real d = cell_at (w2, 1
#if dimension > 1
                          , 0
#endif
#if dimension > 2
                          , 0
#endif
                          ) - cell_at (w2, -1
#if dimension > 1
                          , 0
#endif
#if dimension > 2
                          , 0
#endif
                          );
        diff = fmax (diff, fabs (out [o] - d));
      }
    }
    CHECK (diff == 0);
    /* scatter: doubling the interior of w2 doubles the block */
    foreach_leaf () {
      block_gather (w2, sh, (Index) {IJKL}, u);
      block_interior (sh) w2 [c] *= 2;
      block_scatter ((Index) {IJKL}, u, w2, sh);
    }
    double serr = 0;
    foreach_leaf () {
      Real * b = block_data (sblocks, u);
      block_cells (store, 0) {
        int x = i * m + ci, y = 0, z = 0;
#if dimension > 1
        y = j * m + cj;
#endif
#if dimension > 2
        z = k * m + ck;
#endif
        serr = fmax (serr, fabs (b [c] - 2 * f (level, x, y, z)));
      }
    }
    CHECK (serr < 1e-5);
    free (w2);
  }

  /* the parallel loop gives what the serial one gives, into another
   * scalar, and scratch blocks are aligned */
  {
    scalar ser, par;
    double d = 0;
    foreach_leaf_gathered (u, 2) {
      CHECK (((size_t) w & 63) == 0);
      Real * out = block_data (sblocks, ser);
      block_interior (shape)
        out [o] = cell_at (w, -2
#if dimension > 1
                           , 1
#endif
#if dimension > 2
                           , -1
#endif
                           ) + 3 * cell_at (w, 0
#if dimension > 1
                           , 0
#endif
#if dimension > 2
                           , 0
#endif
                           );
    }
    foreach_leaf_gathered_parallel (u, 2) {
      Real * out = block_data (sblocks, par);
      block_interior (shape)
        out [o] = cell_at (w, -2
#if dimension > 1
                           , 1
#endif
#if dimension > 2
                           , -1
#endif
                           ) + 3 * cell_at (w, 0
#if dimension > 1
                           , 0
#endif
#if dimension > 2
                           , 0
#endif
                           );
    }
    BlockShape store = tree_block_shape ();
    foreach_leaf () {
      const Real * a = block_data (sblocks, ser), * b = block_data (sblocks, par);
      block_cells (store, 0) d = fmax (d, fabs (a [c] - b [c]));
    }
    CHECK (d == 0);
    /* a stored block needs no padding to be aligned: g = 0, m a multiple of V */
    BlockShape al = block_shape_aligned (tree_block_m (), 0, 8);
    CHECK (al.cells == store.cells && al.ox == 0 && al.nx == tree_block_m ());
  }

  /* blocks freed and reused: the hot list, with the same accounting */
  {
    long before = grid_scalar_blocks_attached ();
    { scalar t; (void) t; }               /* blocks for every node, then freed */
    CHECK (grid_scalar_blocks_attached () == before && grid_scalar_check ());
    { scalar t; (void) t; }
    CHECK (grid_scalar_check ());
    grid_scalar_trim (0);
    CHECK (grid_scalar_check ());
  }
  tree_destroy ();
  printf ("dimension %d: %d checks, %d failed\n", dimension, checks, fails);
  return fails != 0;
}
