/* block-loops.h: every back end against block_cells, in both layouts and
 * every halo, with the padding untouched; and the leaf loops on a tree */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "block-loops.h"

static int checks, fails;
#define CHECK(c) do { checks++; if (!(c)) { fails++; \
  printf ("FAIL line %d: %s\n", __LINE__, #c); } } while (0)

#if dimension == 1
# define LAP(u) (cell_at (u, 1) + cell_at (u, -1) - 2 * cell_at (u, 0))
#elif dimension == 2
# define LAP(u) (cell_at (u, 1, 0) + cell_at (u, -1, 0) + cell_at (u, 0, 1) \
                 + cell_at (u, 0, -1) - 4 * cell_at (u, 0, 0))
#else
# define LAP(u) (cell_at (u, 1, 0, 0) + cell_at (u, -1, 0, 0) + cell_at (u, 0, 1, 0) \
                 + cell_at (u, 0, -1, 0) + cell_at (u, 0, 0, 1) + cell_at (u, 0, 0, -1) \
                 - 6 * cell_at (u, 0, 0, 0))
#endif

static long ipow (long x, int e) { long r = 1; while (e--) r *= x; return r; }

static void layout (BlockShape shape)
{
  long n = shape.cells;
  Real * u = malloc (n * sizeof (Real)), * ref = malloc (n * sizeof (Real));
  Real * v = malloc (n * sizeof (Real));
  const Real sentinel = R (-12345);
#define RESET(x) for (long q = 0; q < n; q++) (x) [q] = sentinel
  RESET (u);
  long visited = 0;
  block_cells (shape, shape.g) {
    u [c] = (Real) sin (0.37 * c);
    visited++;
  }
  CHECK (visited == ipow (shape.m + 2 * shape.g, dimension));
  long pad = 0;
  for (long q = 0; q < n; q++) pad += u [q] == sentinel;
  CHECK (pad == n - visited);

  for (int h = 0; h < shape.g; h++) {
    RESET (ref);
    block_cells (shape, h) ref [c] = LAP (u);
    RESET (v); block_cells_simd (shape, h) v [c] = LAP (u);
    CHECK (!memcmp (v, ref, n * sizeof (Real)));
    RESET (v);
    block_rows (shape, h)
      for (int x = 0; x < n; x++) { int c = row + x; v [c] = LAP (u); }
    CHECK (!memcmp (v, ref, n * sizeof (Real)));
    RESET (v); block_cells_tiled (shape, h, 3) v [c] = LAP (u);
    CHECK (!memcmp (v, ref, n * sizeof (Real)));
    RESET (v); block_cells_threads (shape, h) v [c] = LAP (u);
    CHECK (!memcmp (v, ref, n * sizeof (Real)));
    RESET (v); foreach_cell (shape, h) v [c] = LAP (u);
    CHECK (!memcmp (v, ref, n * sizeof (Real)));
    Real * us [2] = { u, u }, * vs [2] = { malloc (n * sizeof (Real)), malloc (n * sizeof (Real)) };
    RESET (vs [0]); RESET (vs [1]);
    gpu_cells (shape, 2, h) { const Real * uu = us [blk]; vs [blk] [c] = LAP (uu); }
    CHECK (!memcmp (vs [0], ref, n * sizeof (Real)) && !memcmp (vs [1], ref, n * sizeof (Real)));
    free (vs [0]); free (vs [1]);
  }
  /* block_at: from interior coordinates, ghosts included, to c */
  long ok = 0, cells = 0;
  block_cells (shape, shape.g) {
    int x [3] = { ci - shape.ox, 0, 0 };
#if dimension > 1
    x [1] = cj - shape.g;
#endif
#if dimension > 2
    x [2] = ck - shape.g;
#endif
    ok += block_at (shape, x) == c;
    cells++;
  }
  CHECK (ok == cells && cells == ipow (shape.m + 2 * shape.g, dimension));
  long k = 0;
  block_cells (shape, 0) { if (k == 3) break; k++; }
  CHECK (k == 3);
  free (u); free (ref); free (v);
}

scalar s;

int main (void)
{
  for (int m = 4; m <= 32; m *= 2) {
    layout (block_shape (m, 2));
    for (int vl = 2; vl <= 8; vl *= 2)
      layout (block_shape_aligned (m, 2, vl));
  }
  /* the leaf loops, on stored blocks, which hold GHOSTS layers */
  tree (3, 4);
  BlockShape store = tree_block_shape ();
  CHECK (store.g == GHOSTS && store.m == tree_get_block_m ());
  CHECK (store.ox % BAMR_ALIGN == 0 && store.nx % BAMR_ALIGN == 0);
  int leaves = 0, par = 0;
  foreach_leaf () {
    Real * u = block_data (sblocks, s);
    block_cells_simd (store, 0) u [c] = level;
    leaves++;
  }
  CHECK (leaves == tree_get_nleaves ());
  foreach_leaf_parallel () {
    Real * u = block_data (sblocks, s);
    RealSum sum = 0;
    block_cells (store, 0) sum += u [c];
    BLOCK_PRAGMA (omp atomic)
    par += sum == ipow (store.m, dimension) * level;
  }
  CHECK (par == tree_get_nleaves ());
  Real * table [64];
  CHECK (block_table (s, table, NULL) == tree_get_nleaves ());
  tree_destroy ();
  printf ("dimension %d: %d checks, %d failed\n", dimension, checks, fails);
  return fails != 0;
}
