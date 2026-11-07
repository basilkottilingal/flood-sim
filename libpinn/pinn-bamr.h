#ifndef _PINN_BAMR_H_
#define _PINN_BAMR_H_

  /*
  .. libb-amr <-> libpinn: the leaf cells of the tree as training data for
  .. a Pinn, and the network's values back into scalars. Built by cwrap,
  .. with -dimension=D, as libb-amr is. Python stays behind pinn.h.
  ..
  .. A cell of a leaf block at level l, with m = tree_get_block_m (), is
  .. d = h0 / (m 2^l) wide, and holds a finite volume value: the mean of u
  .. over the cell, not u at its centre. Two ways to fit it:
  ..
  ..   PINN_SAMPLE_CENTRE   one point a cell, its centre: NN (centre) = u,
  ..                        wrong by O (d^2) where u curves
  ..   PINN_SAMPLE_GAUSS2   the 2^D points of the 2-point Gauss rule: the
  ..                        mean of NN over them = u, the cell average to
  ..                        O (d^4). 2^D times the points.
  ..
  .. The same rule gives the predictions written back, so that a scalar
  .. filled by pinn_bamr_predict () is again a cell average.
  ..
  .. An adaptive tree has far more cells where it is refined. Weights:
  ..
  ..   PINN_WEIGHT_AREA     w ~ d^D: the loss is the L2 norm over the domain,
  ..                        every region counts by its area
  ..   PINN_WEIGHT_UNIFORM  w = 1 / n: every cell counts alike, so refined
  ..                        regions, fronts, count the most
  ..
  .. Collocation points, for a residual (pinn.h): a number of uniform random
  .. points in each leaf block. Finer blocks are smaller, so the points
  .. crowd where the tree refined, at no cost: the tree is the sampler. With
  .. PINN_WEIGHT_AREA, each point weighs its block's area over the number,
  .. which keeps the residual term the integral over the domain.
  ..
  .. With time: with_time appends t as the last input, nin = D + 1, so
  .. snapshots at several t train u (x, t).
  */

  #include "pinn.h"

  enum
  {
    PINN_SAMPLE_CENTRE = 0,
    PINN_SAMPLE_GAUSS2 = 1,
  };

  enum
  {
    PINN_WEIGHT_AREA    = 0,
    PINN_WEIGHT_UNIFORM = 1,
  };

  typedef struct
  {
    double x0 [3];    /* the lower corner, as given tree_shift_resize ()  */
    double h0;        /* the side of the domain: [x0, x0 + h0]^D          */
    int    sample;    /* PINN_SAMPLE_*                                    */
    int    weight;    /* PINN_WEIGHT_*                                    */
    int    with_time; /* append t to the inputs                           */
    double t;
  } PinnGrid;

  /* the unit box [0, 1]^D, the default of libb-amr, centre, area */
  PinnGrid pinn_grid (void);

  /* inputs of the grid: dimension, + 1 with time */
  int  pinn_grid_nin (const PinnGrid * g);

  /*
  .. a config for the grid: nin from it, nout outputs, lo and hi the box
  .. of the domain. With time, the box of t is [t0, t1].
  */
  PinnConfig pinn_grid_config (const PinnGrid * g, const char * backend,
                               int nout, double t0, double t1);

  /*
  .. b->x, b->u, b->w: every cell of every leaf, in Morton order of the
  .. blocks and i fastest in a block, with the values of the nf scalars
  .. fields [0..nf) as outputs, in that order (nf = the Pinn's nout).
  .. Allocates them, after freeing what b had; b->nc, xc, wc are kept.
  .. Returns n, the cells, or -1 if out of memory.
  */
  int  pinn_bamr_gather (const PinnGrid * g, int nf, const int * fields,
                         PinnBatch * b);

  /*
  .. b->xc, b->wc: per_block uniform random points in each leaf block, from
  .. seed. Allocates them, after freeing what b had. Returns nc, or -1.
  */
  int  pinn_bamr_collocation (const PinnGrid * g, int per_block,
                              unsigned long seed, PinnBatch * b);

  /* frees the arrays of b, and zeroes it */
  void pinn_batch_free (PinnBatch * b);

  /*
  .. fields [0..nf) = the network in every leaf cell, by the sample rule of
  .. g: the cell average with PINN_SAMPLE_GAUSS2. nf = nout.
  */
  int  pinn_bamr_predict (Pinn * p, const PinnGrid * g, int nf,
                          const int * fields);

  /*
  .. fields [0..nf) = the residual at the centre of every leaf cell,
  .. nf = nout. Its norm in a cell is an indicator of where the network
  .. breaks the physics: for adapt_wavelet (), or a refinement of its own.
  */
  int  pinn_bamr_residual (Pinn * p, const PinnGrid * g, int nf,
                           const int * fields);

#endif
