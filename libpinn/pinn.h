#ifndef _PINN_H_
#define _PINN_H_

  /*
  .. libpinn: a C interface to a neural network trained in Python, by
  .. PyTorch or JAX. The interpreter is embedded (CPython), and each Pinn
  .. is one Python object of the backend module (python/pinn_backend.py),
  .. so the C side never knows which framework runs.
  ..
  .. The network is u = NN (x): nin inputs (x, y [, z] [, t]), nout outputs.
  .. A fit is a few optimiser steps on
  ..
  ..   L = lambda_data * sum_i w_i |mean_q NN (x_iq) - u_i|^2
  ..     + lambda_phys * sum_j wc_j |r (NN, xc_j)|^2
  ..
  .. q runs over the nq points of sample i: nq = 1 for a point value, or the
  .. quadrature points of a cell, whose mean is then the network's cell
  .. average, to compare with a finite volume value. r is the residual
  .. function the config names, written in Python with the backend's own
  .. autodiff. Fits warm start: the weights and the optimiser state stay
  .. from one call to the next, so each new snapshot of a solver reinforces
  .. the same network.
  ..
  .. This header knows nothing of the tree: pinn-bamr.h gathers the cells of
  .. libb-amr into a PinnBatch, and writes predictions back to scalars.
  ..
  .. Every call that can fail returns 0 on success and -1 on error, and
  .. pinn_error () then gives the message (a Python traceback, if Python
  .. raised). Calls are serialised by the GIL: safe from any thread, one at
  .. a time.
  */

  #include <stddef.h>

  #define PINN_MAX_IN    4
  #define PINN_MAX_OUT   16

  typedef struct
  {
    const char * backend;     /* "torch" or "jax"                          */
    int    nin, nout;         /* inputs, outputs                           */
    int    width, depth;      /* hidden layers: depth of width neurons     */
    const char * activation;  /* "tanh" (default), "silu", "sin"           */
    int    fourier;           /* random Fourier features, 0 for none       */
    double fourier_sigma;     /* their scale, in normalised inputs         */
    double lr;                /* Adam learning rate                        */
    double lambda_data;       /* weight of the data term, 0 means 1        */
    double lambda_phys;       /* weight of the residual term               */
    double lo [PINN_MAX_IN];  /* the box of the inputs, mapped to [-1, 1]  */
    double hi [PINN_MAX_IN];  /* inside the network                        */
    const char * residual;    /* "module:function", or NULL for none       */
    const char * device;      /* "cpu" (default), "cuda", ...              */
    const char * dtype;       /* "float32" (default) or "float64"          */
    int    seed;
  } PinnConfig;

  /*
  .. A batch of samples, in arrays the caller owns, row major:
  ..   x   [n * nq * nin]  the points, nq per sample
  ..   u   [n * nout]      the value of each sample: the mean over its points
  ..   w   [n]             its weight, or NULL for 1 / n
  ..   xc  [nc * nin]      collocation points, for the residual; nc may be 0
  ..   wc  [nc]            their weights, or NULL for 1 / nc
  .. Geometry is double, as libb-amr keeps it; fields are float.
  */
  typedef struct
  {
    int n, nq, nc;
    double * x;
    float  * u;
    float  * w;
    double * xc;
    float  * wc;
  } PinnBatch;

  typedef struct
  {
    double total, data, phys;  /* the loss of the last step, by term */
  } PinnLoss;

  typedef struct Pinn Pinn;

  /* a config with defaults: tanh, 4 x 64, lr 1e-3, box [-1, 1]^nin */
  PinnConfig pinn_config (const char * backend, int nin, int nout);

  Pinn * pinn_new (const PinnConfig * cfg);   /* NULL on error: see below */
  void   pinn_free (Pinn * p);

  /*
  .. steps optimiser steps on the whole batch (full batch Adam). The first
  .. fit of a Pinn also sets the scaling of the outputs, from the mean and
  .. deviation of b->u, and keeps it after.
  */
  int    pinn_fit (Pinn * p, const PinnBatch * b, int steps, PinnLoss * loss);

  /* u [n * nout] = NN (x [n * nin]) */
  int    pinn_predict (Pinn * p, int n, const double * x, float * u);

  /* r [n * nout] = the residual at x [n * nin], if the config names one */
  int    pinn_residual (Pinn * p, int n, const double * x, float * r);

  /* weights, optimiser and scaling, to a file and back */
  int    pinn_save (Pinn * p, const char * path);
  int    pinn_load (Pinn * p, const char * path);

  /* the last error of p, or of pinn_new () when p is NULL */
  const char * pinn_error (const Pinn * p);

  /*
  .. The interpreter starts with the first pinn_new () and stays, as torch
  .. and numpy cannot be loaded twice in one process. pinn_shutdown () ends
  .. it, once every Pinn is freed; after that, no more Pinn. If the program
  .. is itself Python (libpinn loaded by ctypes), it is never ended.
  */
  void   pinn_shutdown (void);

#endif
