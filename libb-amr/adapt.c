#include <assert.h>
#include <math.h>
#include <stdlib.h>

#include "adapt.h"
#include "block-ghosts.h"

double block_wavelet (Index p, int s)
{
  int m = tree_get_block_m (), h = m / 2, nc = h + 2;
  assert (GHOSTS >= 2 && m >= 4);
  BlockShape sh = tree_scratch_shape (2);
  Real * u = block_scratch_new (sh);
  block_gather (u, sh, p, s);

  /*
  .. The coarse cells, -1 to h along each axis: the means of 2 x 2 (x 2)
  .. fine cells, the ghost layer of coarse cells from the two of fine ones.
  */
  long cs [3] = { 1, dimension > 1 ? nc : 0, dimension > 2 ? (long) nc * nc : 0 };
  long ncells = nc;
  for (int a = 1; a < dimension; a++)
    ncells *= nc;
  RealSum * coarse = malloc ((size_t) ncells * sizeof (RealSum));
  int lo [3] = { -1, dimension > 1 ? -1 : 0, dimension > 2 ? -1 : 0 };
  int hi [3] = { h + 1, dimension > 1 ? h + 1 : 1, dimension > 2 ? h + 1 : 1 };
  int C [3];
  for (C [2] = lo [2]; C [2] < hi [2]; C [2]++)
  for (C [1] = lo [1]; C [1] < hi [1]; C [1]++)
  for (C [0] = lo [0]; C [0] < hi [0]; C [0]++)
  {
    RealSum sum = 0;
    for (int k = 0; k < TREE_NCHILDREN; k++)
    {
      int x [3] = { 0, 0, 0 };
      for (int a = 0; a < dimension; a++)
        x [a] = 2 * C [a] + ((k >> a) & 1);
      sum += u [(x [0] + sh.ox) + (long) (x [1] + sh.g) * sh.sy
                + (long) (x [2] + sh.g) * sh.sz];
    }
    long q = 0;
    for (int a = 0; a < dimension; a++)
      q += (C [a] + 1) * cs [a];
    coarse [q] = sum / TREE_NCHILDREN;
  }

  /* each fine cell against its prediction from the coarse cells */
  double err = 0;
  int fhi [3] = { m, dimension > 1 ? m : 1, dimension > 2 ? m : 1 };
  int x [3];
  for (x [2] = 0; x [2] < fhi [2]; x [2]++)
  for (x [1] = 0; x [1] < fhi [1]; x [1]++)
  for (x [0] = 0; x [0] < fhi [0]; x [0]++)
  {
    long q = 0;
    for (int a = 0; a < dimension; a++)
      q += ((x [a] >> 1) + 1) * cs [a];
    RealSum pred = coarse [q];
    for (int a = 0; a < dimension; a++)
      pred += (x [a] & 1 ? 1 : -1) * (coarse [q + cs [a]] - coarse [q - cs [a]]) / 8.;
    Real v = u [(x [0] + sh.ox) + (long) (x [1] + sh.g) * sh.sy
                + (long) (x [2] + sh.g) * sh.sz];
    err = fmax (err, fabs (v - pred));
  }
  free (coarse);
  free (u);
  return err;
}

static int by_level_up (const void * a, const void * b)
{
  return ((const Index *) a)->level - ((const Index *) b)->level;
}

static int by_level_down (const void * a, const void * b)
{
  return by_level_up (b, a);
}

AdaptStats adapt_wavelet (int n, const int * s, const double * max,
                          int minlevel, int maxlevel)
{
  AdaptStats st = { 0, 0 };
  if (maxlevel > tree_get_depth_max ())
    maxlevel = tree_get_depth_max ();
  if (minlevel < 0)
    minlevel = 0;

  /* the parents, which ghosts read where a neighbour is finer */
  for (int i = 0; i < n; i++)
    tree_restrict (s [i]);

  /* 1. flags */
  cache (Index) refine;
  foreach_leaf () {
    node->flags &= ~(BAMR_NODE_REFINE | BAMR_NODE_COARSEN);
    /* not i: that is the leaf's index, which the block sees */
    int big = 0, small = 1;
    for (int q = 0; q < n; q++) {
      double e = block_wavelet ((Index) {IJKL}, s [q]);
      big |= e > max [q];
      small &= e < max [q] / 1.5;
    }
    if (level < maxlevel && (big || level < minlevel)) {
      node->flags |= BAMR_NODE_REFINE;
      cache_push (refine, (Index) {IJKL});
    }
    else if (level > maxlevel || (small && level > minlevel))
      node->flags |= BAMR_NODE_COARSEN;       /* above maxlevel: always */
  }

  /* 2. refine, coarsest first: a cascade for balance then goes down */
  if (refine.n)
    qsort (refine.a, refine.n, sizeof (Index), by_level_up);
  foreach (p in refine)
    st.refined += tree_refine (p);

  /*
  .. 3. coarsen, finest first: families whose children are all leaves
  .. flagged. A child refined in 2, or made for balance, is not one.
  ..
  .. And the parent must not be too coarse itself: its own wavelet at most
  .. max, as for a leaf that would not be refined. Otherwise a leaf refined
  .. for an error e > max has children of error about e / 4, which for e
  .. below 2.7 max would be flagged, coarsened, and the leaf refined again,
  .. and so on: the parent's error is what decided the refinement, so it
  .. decides the coarsening too.
  */
  cache (Index) merge;
  morton (TREE_ROOT_INDEX) {
    if (!tree_has_children (node) || level < minlevel)
      continue;
    int all = 1;
    children_of ((Index) {IJKL})
      all &= (node->flags & BAMR_NODE_LEAF) && (node->flags & BAMR_NODE_COARSEN);
    for (int q = 0; all && q < n && level + 1 <= maxlevel; q++)
      all &= block_wavelet ((Index) {IJKL}, s [q]) <= max [q];
    if (all)
      cache_push (merge, (Index) {IJKL});
  }
  if (merge.n)
    qsort (merge.a, merge.n, sizeof (Index), by_level_down);
  foreach (p in merge)
    st.coarsened += tree_coarsen (p);

  /* 4. clear the flags */
  morton (TREE_ROOT_INDEX)
    node->flags &= ~(BAMR_NODE_REFINE | BAMR_NODE_COARSEN);
  return st;
}
