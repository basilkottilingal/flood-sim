/* a hash of every cell of every block, ghosts included, after refining,
 * ghosts_update () with given operators, and one adapt_wavelet () */
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "adapt.h"
#include "block-ghosts.h"

scalar a, b, h;

static unsigned long long hash = 1469598103934665603ull;
static void add (const void * p, size_t n)
{
  const unsigned char * c = p;
  for (size_t k = 0; k < n; k++)
    hash = (hash ^ c [k]) * 1099511628211ull;
}

static void set_fields (void)
{
  BlockShape shape = tree_block_shape ();
  foreach_leaf () {
    Real * pa = block_data (sblocks, a), * pb = block_data (sblocks, b);
    Real * ph = block_data (sblocks, h);
    double n = (double) (shape.m << level);
    block_cells (shape, 0) {
      double x = (i * shape.m + ci - shape.ox + 0.5) / n, y = 0.3;
#if dimension > 1
      y = (j * shape.m + cj - shape.g + 0.5) / n;
#endif
      pa [c] = (Real) sin (6 * x + 2 * y);
      pb [c] = (Real) (x * y - 0.2);
      ph [c] = (Real) fmax (0, 0.1 - fabs (x - 0.4));
    }
  }
}

static void hash_all (void)
{
  BlockShape shape = tree_block_shape ();
  morton (TREE_ROOT_INDEX) {
    add (&level, sizeof level);
    add (block_data (sblocks, a), shape.cells * sizeof (Real));
    add (block_data (sblocks, b), shape.cells * sizeof (Real));
    add (block_data (sblocks, h), shape.cells * sizeof (Real));
  }
}

int main (void)
{
  tree (5, 4);
  for (int r = 0; r < 3; r++) {
    Index first = { 0 };
    int got = 0;
    foreach_leaf ()
      if (!got && level == r + 1)
        first = (Index) {IJKL}, got = 1;
    if (!got)
      foreach_leaf ()
        if (!got)
          first = (Index) {IJKL}, got = 1;
    tree_refine (first);
  }
  scalar_set_ops (h, (ScalarOps) { .prolongation = prolong_minmod });
  scalar_set_ops (b, (ScalarOps) { .prolongation = prolong_injection });
  ghost_ops (b, boundary.all (bc_odd));
  set_fields ();
  ghosts_update_levels (a, b, h);
  hash_all ();
  int s [3] = { a, b, h };
  double max [3] = { 1e-2, 1e-2, 1e-3 };
  AdaptStats st = adapt_wavelet (3, s, max, 1, 5);
  add (&st, sizeof st);
  hash_all ();
  printf ("%d %d %016llx\n", tree_get_nleaves (), tree_get_nblocks (), hash);
  tree_destroy ();
  return 0;
}
