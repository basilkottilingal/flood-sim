#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <cwrap/cache.h>

#include "adapt.h"
#include "block-ghosts.h"

/* the wavelet of p for s, scratch a block to prolong into */
static double wavelet_into (Index p, int s, Real * scratch)
{
  BlockShape shape = tree_block_shape ();
  block_prolong_into (p, s, scratch);
  const Real * u = block_data (index_node (p)->s, s);
  /* a select, not fmax (): twice as fast, the same result, a NaN ignored
  .. by both */
  double err = 0;
  block_cells (shape, 0) {
    double d = (double) u [c] - scratch [c];
    d = d < 0 ? -d : d;
    err = d > err ? d : err;
  }
  return err;
}

static Real * scratch_block (void)
{
  Real * b = malloc ((size_t) tree_block_shape ().cells * sizeof (Real));
  if (!b)
    abort ();
  return b;
}

double block_wavelet (Index p, int s)
{
  assert (p.level > 0 && index_is_allocated (p));
  Real * scratch = scratch_block ();
  double e = wavelet_into (p, s, scratch);
  free (scratch);
  return e;
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
  /* a scalar with prolong_none has no wavelet of its own: another's
   * operators fill it (block-ghosts.h). Refused: the run stops, or with
   * ghost_ops_abort 0 nothing is adapted */
  for (int q = 0; q < n; q++)
    if (scalar_get_ops (s [q]).prolongation == prolong_none) {
      fprintf (stderr, "adapt_wavelet: scalar %s (%d): its prolongation is "
               "prolong_none, so it has no wavelet of its own: adapt on the "
               "scalar whose operators fill it\n",
               *scalar_name (s [q]) ? scalar_name (s [q]) : "#unnamed",
               (int) s [q]);
      if (ghost_ops_abort)
        abort ();
      return st;              /* nothing adapted */
    }
  if (maxlevel > tree_get_depth_max ())
    maxlevel = tree_get_depth_max ();
  if (minlevel < 0)
    minlevel = 0;

  /*
  .. The ghosts of every node, parents too, as the wavelet of a parent
  .. decides its coarsening; and of every scalar, as refining prolongs them
  .. all from their ghosts. tree_refine () and tree_coarsen () keep them
  .. valid from then on.
  */
  ghosts_update_scalars (1);
  Real * scratch = scratch_block ();

  /*
  .. 1. flags. The root has no coarser level to be predicted from, so no
  .. wavelet: as a leaf (after tree_coarsen () of it) it is refined, as
  .. tree () makes it, and in 3 no family is coarsened into it.
  */
  cache (Index) refine;
  foreach_leaf () {
    node->flags &= ~(BAMR_NODE_REFINE | BAMR_NODE_COARSEN);
    /* not i: that is the leaf's index, which the block sees */
    int big = level == 0, small = level > 0;
    for (int q = 0; level > 0 && q < n; q++) {
      double e = wavelet_into ((Index) {IJKL}, s [q], scratch);
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
  foreach_node () {
    if (!tree_has_children (node) || level < minlevel || level == 0)
      continue;
    int all = 1;
    children_of ((Index) {IJKL})
      all &= (node->flags & BAMR_NODE_LEAF) && (node->flags & BAMR_NODE_COARSEN);
    for (int q = 0; all && q < n && level + 1 <= maxlevel; q++)
      all &= wavelet_into ((Index) {IJKL}, s [q], scratch) <= max [q];
    if (all)
      cache_push (merge, (Index) {IJKL});
  }
  if (merge.n)
    qsort (merge.a, merge.n, sizeof (Index), by_level_down);
  foreach (p in merge)
    st.coarsened += tree_coarsen (p);

  /* 4. clear the flags */
  foreach_node ()
    node->flags &= ~(BAMR_NODE_REFINE | BAMR_NODE_COARSEN);
  free (scratch);
  return st;
}
