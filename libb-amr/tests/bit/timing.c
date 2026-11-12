/* repeated ghosts_update_levels (): restriction, prolongation at level
 * jumps, boundary conditions; the best of several runs */
#include <stdio.h>
#include <time.h>
#include "adapt.h"
#include "block-ghosts.h"

scalar a, b, h;

static double now (void)
{
  struct timespec t;
  clock_gettime (CLOCK_MONOTONIC, &t);
  return t.tv_sec + 1e-9 * t.tv_nsec;
}

int main (void)
{
  tree (5, 4);
  for (int r = 0; r < 40; r++) {
    Index first = { 0 };
    int got = 0;
    foreach_leaf ()
      if (!got && level < 4 && (r * 7 + i) % 3 == 0)
        first = (Index) {IJKL}, got = 1;
    if (got)
      tree_refine (first);
  }
  scalar_set_ops (h, (ScalarOps) { .prolongation = prolong_minmod });
  BlockShape shape = tree_block_shape ();
  foreach_leaf () {
    Real * pa = block_data (sblocks, a), * pb = block_data (sblocks, b);
    Real * ph = block_data (sblocks, h);
    block_cells (shape, 0)
      pa [c] = (Real) (c % 17), pb [c] = (Real) (c % 5), ph [c] = (Real) (c % 3);
  }
  double best = 1e30;
  for (int rep = 0; rep < 7; rep++) {
    double t0 = now ();
    for (int k = 0; k < 10; k++)
      ghosts_update_levels (a, b, h);
    double t = now () - t0;
    if (t < best)
      best = t;
  }
  printf ("%d blocks, %.3f ms per update\n", tree_get_nblocks (), 1e3 * best / 10);
  tree_destroy ();
  return 0;
}
