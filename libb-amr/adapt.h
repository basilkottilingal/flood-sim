#ifndef _BAMR_ADAPT_H_
#define _BAMR_ADAPT_H_

  /*
  .. Adapting the tree to the fields it holds, keeping it 2:1 balanced.
  ..
  .. The error of a block is its wavelet: how far each of its cells is from
  .. what the coarser level would predict for it, by the scalar's own
  .. prolongation from the parent's cells (block_prolong_into ()): the
  .. parent's quadrant over the block is its restriction, the layer around
  .. it the parent's other cells and first ghost layer. The error is the
  .. largest |cell - prediction|: zero where the prolongation is exact (a
  .. linear field for prolong_linear (), a flat surface for a lake-at-rest
  .. one), large where the coarse level would lose it. It reads the parent
  .. and its ghosts as they are: ghosts_update_levels () first. The root has
  .. none: adapt_wavelet () refines it as a leaf, and never coarsens a
  .. family into it, so the tree keeps level 1.
  ..
  .. adapt_wavelet (n, s, max, minlevel, maxlevel) then, for scalars s [0..n):
  ..
  ..   0. fills the ghosts of every scalar, at every node
  ..      (ghosts_update_scalars (1)): the values may have changed since the
  ..      last step;
  ..   1. flags every leaf: BAMR_NODE_REFINE where some error > max [i]
  ..      (or the leaf is coarser than minlevel), BAMR_NODE_COARSEN where
  ..      every error < max [i] / 1.5, so a leaf just past the threshold
  ..      does not swing between the two, and wherever the leaf is finer
  ..      than maxlevel;
  ..   2. refines the flagged leaves, coarsest first, below maxlevel:
  ..      tree_refine () refines first any coarser leaf beside them;
  ..   3. coarsens, finest first, each family whose children are all leaves
  ..      flagged, and whose parent is not coarser than minlevel, nor too
  ..      coarse itself (its own wavelet at most max [i]: the error that
  ..      refined it, so a refinement is not undone the next step), unless
  ..      the children are finer than maxlevel:
  ..      tree_coarsen () refuses one that would break the balance;
  ..   4. clears the flags.
  ..
  .. One call moves a leaf by one level at most, as a family can only be
  .. coarsened if its children were leaves when flagged. Call it once per
  .. step, as the fields change.
  ..
  .. New children get values by their scalar's prolongation, and a
  .. coarsened parent its restriction: with the defaults (linear, the mean)
  .. both keep the integral of every scalar. On return, the ghosts of every
  .. leaf are valid.
  */

  #include "tree.h"

  typedef struct
  {
    int refined;      /* leaves divided, those for balance included */
    int coarsened;    /* families merged                             */
  } AdaptStats;

  AdaptStats adapt_wavelet (int n, const int * s, const double * max,
                            int minlevel, int maxlevel);

  /* the wavelet error of the node p, not the root, for scalar s: see
  .. above */
  double block_wavelet (Index p, int s);

#endif
