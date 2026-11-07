/*
.. libpinn on a libb-amr quadtree: gathering cells, fitting two fields, the
.. predictions back into scalars, warm starts, residuals by autodiff, the
.. collocation term, save and load.
..
..   test_pinn torch|jax
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "adapt.h"
#include "block-ghosts.h"
#include "pinn-bamr.h"

static int checks, fails;
#define CHECK(c) do { checks++; if (!(c)) { fails++; \
  printf ("FAIL line %d: %s\n", __LINE__, #c); } } while (0)

scalar u, v, pu, pv, ru, rv;

#define X0  (-1.)
#define H0  2.

/* the fields, on [-1, 1]^D */
static double fu (const double * x)
{
  double r = 0;
  for (int a = 0; a < dimension; a++)
    r += x [a] * x [a];
  return tanh ((sqrt (r) - 0.5) / 0.06);
}

static double fv (const double * x)
{
  return 2. + sin (M_PI * x [0]) * (dimension > 1 ? cos (0.5 * M_PI * x [1]) : 1.);
}

/* u, v at the cell centres */
static void set (void)
{
  BlockShape st = tree_block_shape ();
  int m = tree_get_block_m ();
  foreach_leaf () {
    Real * bu = block_data (sblocks, u), * bv = block_data (sblocks, v);
    double d = H0 / (double) ((long) m << level);
    block_cells (st, 0) {
      double x [3] = { X0 + (i * m + ci + 0.5) * d, 0, 0 };
#if dimension > 1
      x [1] = X0 + (j * m + cj + 0.5) * d;
#endif
#if dimension > 2
      x [2] = X0 + (k * m + ck + 0.5) * d;
#endif
      bu [c] = (Real) fu (x);
      bv [c] = (Real) fv (x);
    }
  }
}

/* the relative L2 error of field p against field s, by area */
static double rel_error (int p, int s)
{
  BlockShape st = tree_block_shape ();
  double e = 0, n = 0;
  foreach_leaf () {
    const Real * a = block_data (sblocks, p), * b = block_data (sblocks, s);
    double w = pow (0.5, dimension * level);
    block_cells (st, 0) {
      e += w * (a [c] - b [c]) * (a [c] - b [c]);
      n += w * b [c] * b [c];
    }
  }
  return sqrt (e / n);
}

static Real max_abs (int s)
{
  Real mx = 0;
  BlockShape st = tree_block_shape ();
  foreach_leaf () {
    const Real * a = block_data (sblocks, s);
    block_cells (st, 0)
      if (fabs (a [c]) > mx)
        mx = fabs (a [c]);
  }
  return mx;
}

#define TRY(call) do { if ((call) < 0) { \
  printf ("FAIL line %d: %s\n%s\n", __LINE__, #call, pinn_error (p)); \
  exit (1); } } while (0)

int main (int argc, char ** argv)
{
  const char * backend = argc > 1 ? argv [1] : "torch";
  tree (5, 4);
  tree_shift_resize (X0, X0, X0, H0);

  /* an adapted tree: the front of u refined, to level 3 (5 in 1D) */
  {
    double max [2] = { dimension == 1 ? 1e-4 : 1e-2, 1. };
    AdaptStats st;
    int its = 0;
    do {
      set ();
      st = adapt_wavelet (2, (int []) { u, v }, max, 1,
                          dimension == 1 ? 5 : 3);
    } while ((st.refined || st.coarsened) && ++its < 20);
    set ();
    CHECK (tree_get_depth () > tree_get_depth_min ());
    printf ("  tree: %d leaves, levels %d to %d\n", tree_get_nleaves (),
            tree_get_depth_min (), tree_get_depth ());
  }

  PinnGrid g = pinn_grid ();
  for (int a = 0; a < 3; a++)
    g.x0 [a] = X0;
  g.h0 = H0;
  g.sample = PINN_SAMPLE_GAUSS2;
  PinnBatch b;
  memset (&b, 0, sizeof (b));

  /* gathering: every cell, its Gauss points around its centre, weights */
  {
    long m = tree_get_block_m (), cells = tree_get_nleaves ();
    for (int a = 0; a < dimension; a++)
      cells *= m;
    int n = pinn_bamr_gather (&g, 2, (int []) { u, v }, &b);
    CHECK (n == cells && b.nq == 1 << dimension);
    /* sum of w d_min^D, d_min the coarsest cell: the area of the domain */
    double dmax = H0 / (double) (m << tree_get_depth_min ()), area = 0;
    for (int s = 0; s < n; s++)
      area += b.w [s] * pow (dmax, dimension);
    CHECK (fabs (area - pow (H0, dimension)) < 1e-9);
    int inside = 1;
    for (long s = 0; s < (long) n * b.nq * dimension; s++)
      inside &= b.x [s] > X0 && b.x [s] < X0 + H0;
    CHECK (inside);
    /* cell 0: the mean of its points is its centre, and u its value */
    double c0 [3] = { 0, 0, 0 };
    for (int q = 0; q < b.nq; q++)
      for (int a = 0; a < dimension; a++)
        c0 [a] += b.x [q * dimension + a] / b.nq;
    CHECK (fabs (b.u [0] - fu (c0)) < 1e-6 && fabs (b.u [1] - fv (c0)) < 1e-6);
  }

  PinnConfig cfg = pinn_grid_config (&g, backend, 2, 0, 0);
  cfg.width = 48;
  cfg.depth = 3;
  cfg.lr = 3e-3;
  cfg.residual = "pinn_test_residual:dudx_torch";
  if (!strcmp (backend, "jax"))
    cfg.residual = "pinn_test_residual:dudx_jax";
  Pinn * p = pinn_new (&cfg);
  if (!p)
  {
    printf ("FAIL pinn_new:\n%s\n", pinn_error (NULL));
    return 1;
  }

  /* a bad config fails, with a message, and does not break the next */
  {
    PinnConfig bad = cfg;
    bad.backend = "tensorflow";
    Pinn * q = pinn_new (&bad);
    CHECK (!q && strstr (pinn_error (NULL), "unknown backend"));
  }

  /*
  .. fitting, on the cell centres (a quarter of the Gauss points in 2D, to
  .. keep the test short): the loss falls, and the predictions, cell
  .. averages by the Gauss rule, come back close
  */
  PinnLoss l0, l1, l2;
  {
    PinnGrid gc = g;
    gc.sample = PINN_SAMPLE_CENTRE;
    CHECK (pinn_bamr_gather (&gc, 2, (int []) { u, v }, &b) > 0 && b.nq == 1);
    TRY (pinn_fit (p, &b, 1, &l0));
    TRY (pinn_fit (p, &b, 1200, &l1));
    CHECK (l1.data < 0.02 * l0.data && l1.phys == 0);
    TRY (pinn_bamr_predict (p, &g, 2, (int []) { pu, pv }));
    double eu = rel_error (pu, u), ev = rel_error (pv, v);
    CHECK (eu < 0.1 && ev < 0.02);
    printf ("  fit: loss %.3g -> %.3g, relative L2 error u %.3g, v %.3g\n",
            l0.data, l1.data, eu, ev);
  }

  /* warm start: a new fit starts from where the last ended */
  {
    TRY (pinn_fit (p, &b, 1, &l2));
    CHECK (l2.data < 2 * l1.data);
  }

  /* the residual, du/dx by autodiff, against differences of predictions */
  {
    double x [2 * 3], eps = 1e-3;
    float up [2 * 2], r [2];
    double x0 [3] = { 0.3, -0.2, 0.1 };
    memcpy (x, x0, sizeof (double) * dimension);
    memcpy (x + dimension, x0, sizeof (double) * dimension);
    x [0] -= eps, x [dimension] += eps;
    TRY (pinn_predict (p, 2, x, up));
    TRY (pinn_residual (p, 1, x0, r));
    for (int k = 0; k < 2; k++)
    {
      double fd = (up [2 + k] - up [k]) / (2 * eps);
      CHECK (fabs (r [k] - fd) < 2e-2 * (1 + fabs (fd)));
    }
    TRY (pinn_bamr_residual (p, &g, 2, (int []) { ru, rv }));
    CHECK (max_abs (ru) > 0.5);     /* u climbs the front, d ~ 1/0.06 */
  }

  /* the collocation term: points in each leaf, trained too */
  {
    int nc = pinn_bamr_collocation (&g, 8, 42, &b);
    CHECK (nc == 8 * tree_get_nleaves ());
    PinnLoss lp;
    TRY (pinn_fit (p, &b, 5, &lp));
    CHECK (lp.phys > 0 && fabs (lp.total - (lp.data + lp.phys)) < 1e-4 * lp.total);
  }

  /* save, train away, load: the predictions of the save come back */
  {
    char path [256];
    snprintf (path, sizeof (path), "%s/pinn-%s.ckpt",
              getenv ("TMPDIR") ? getenv ("TMPDIR") : "/tmp", backend);
    double x [3] = { 0.1, 0.2, 0.3 };
    float a [2], c [2];
    TRY (pinn_predict (p, 1, x, a));
    TRY (pinn_save (p, path));
    b.nc = 0;
    TRY (pinn_fit (p, &b, 50, NULL));
    TRY (pinn_load (p, path));
    TRY (pinn_predict (p, 1, x, c));
    CHECK (a [0] == c [0] && a [1] == c [1]);
    remove (path);
  }

  pinn_batch_free (&b);
  pinn_free (p);
  pinn_shutdown ();
  tree_destroy ();
  printf ("%s: %d checks, %d failed\n", backend, checks, fails);
  return fails != 0;
}
