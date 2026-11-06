#ifndef _BAMR_ADAPT_H_
#define _BAMR_ADAPT_H_

  /*
  .. Adapting the tree to the fields it holds, keeping it 2:1 balanced.
  ..
  .. The error of a block is its wavelet: how far each of its cells is from
  .. what the coarser level would predict for it. The block, gathered with
  .. two ghost layers, is averaged into coarse cells, 2^dimension fine cells
  .. each; each fine cell is then predicted from its coarse cell and that
  .. cell's neighbours, by the same linear prolongation as a refinement
  .. (block_prolong ()). The error is the largest |cell - prediction|: zero
  .. where the field is linear, large where the coarse level would lose it.
  ..
  .. adapt_wavelet (n, s, max, minlevel, maxlevel) then, for scalars s [0..n):
  ..
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
  .. New children get values by linear prolongation, and a coarsened parent
  .. the mean of its children: both keep the integral of every scalar.
  */

  #include "tree.h"

  typedef struct
  {
    int refined;      /* leaves divided, those for balance included */
    int coarsened;    /* families merged                             */
  } AdaptStats;

  AdaptStats adapt_wavelet (int n, const int * s, const double * max,
                            int minlevel, int maxlevel);

  /* the wavelet error of the leaf p for scalar s: see above */
  double block_wavelet (Index p, int s);

#endif
