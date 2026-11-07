/*
.. libb-amr <-> libpinn: see pinn-bamr.h
*/
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "tree.h"
#include "block-ghosts.h"
#include "pinn-bamr.h"

/* the offsets of the 2-point Gauss rule, in cell widths from the centre */
#define GAUSS2  0.28867513459481288225    /* 1 / (2 sqrt (3)) */

PinnGrid pinn_grid (void)
{
  PinnGrid g;
  memset (&g, 0, sizeof (g));
  g.h0 = 1.;
  g.sample = PINN_SAMPLE_CENTRE;
  g.weight = PINN_WEIGHT_AREA;
  return g;
}

int pinn_grid_nin (const PinnGrid * g)
{
  return dimension + (g->with_time ? 1 : 0);
}

PinnConfig pinn_grid_config (const PinnGrid * g, const char * backend,
                             int nout, double t0, double t1)
{
  PinnConfig c = pinn_config (backend, pinn_grid_nin (g), nout);
  for (int a = 0; a < dimension; a++)
    c.lo [a] = g->x0 [a], c.hi [a] = g->x0 [a] + g->h0;
  if (g->with_time)
    c.lo [dimension] = t0, c.hi [dimension] = t1;
  return c;
}

static int points_per_cell (const PinnGrid * g)
{
  return g->sample == PINN_SAMPLE_GAUSS2 ? 1 << dimension : 1;
}

/* the cell width at level */
static double cell_width (const PinnGrid * g, int level)
{
  return g->h0 / (double) ((long) tree_get_block_m () << level);
}

/*
.. the nq points of cell (cell) of the block at (b) [dimension], level,
.. into x [nq * nin]: centre, or the Gauss points, bit a of q choosing the
.. side along axis a.
*/
static void cell_points (const PinnGrid * g, const int * b, const int * cell,
                         int level, double * x)
{
  int nin = pinn_grid_nin (g), nq = points_per_cell (g), m = tree_get_block_m ();
  double d = cell_width (g, level);
  for (int q = 0; q < nq; q++)
  {
    for (int a = 0; a < dimension; a++)
    {
      double off = nq == 1 ? 0. : ((q >> a) & 1 ? GAUSS2 : -GAUSS2);
      x [q * nin + a] = g->x0 [a] + ((double) b [a] * m + cell [a] + 0.5 + off) * d;
    }
    if (g->with_time)
      x [q * nin + dimension] = g->t;
  }
}

/* the block indices of the leaf IJKL, and of a cell, as arrays */
#if dimension == 1
  #define BLOCK_IJK(b)   int b [1] = { i }
  #define CELL_IJK(c)    int c [1] = { ci }
#elif dimension == 2
  #define BLOCK_IJK(b)   int b [2] = { i, j }
  #define CELL_IJK(c)    int c [2] = { ci, cj }
#else
  #define BLOCK_IJK(b)   int b [3] = { i, j, k }
  #define CELL_IJK(c)    int c [3] = { ci, cj, ck }
#endif

static long block_cells_count (void)
{
  long n = 1;
  for (int a = 0; a < dimension; a++)
    n *= tree_get_block_m ();
  return n;
}

void pinn_batch_free (PinnBatch * b)
{
  free (b->x), free (b->u), free (b->w), free (b->xc), free (b->wc);
  memset (b, 0, sizeof (*b));
}

/* the points of every leaf cell, x [n * nq * nin]; n into *n */
static double * all_points (const PinnGrid * g, int * n)
{
  long cells = (long) tree_get_nleaves () * block_cells_count ();
  int nin = pinn_grid_nin (g), nq = points_per_cell (g);
  double * x = malloc ((size_t) (cells ? cells : 1) * nq * nin * sizeof (double));
  if (!x)
    return NULL;
  BlockShape st = tree_block_shape ();
  long s = 0;
  foreach_leaf () {
    BLOCK_IJK (bi);
    block_cells (st, 0) {
      CELL_IJK (cc);
      cell_points (g, bi, cc, level, x + s * nq * nin);
      s++;
    }
  }
  *n = (int) s;
  return x;
}

int pinn_bamr_gather (const PinnGrid * g, int nf, const int * fields,
                      PinnBatch * b)
{
  free (b->x), free (b->u), free (b->w);
  b->x = NULL, b->u = NULL, b->w = NULL, b->n = 0;
  int n = 0;
  double * x = all_points (g, &n);
  float * u = malloc ((size_t) (n ? n : 1) * nf * sizeof (float));
  float * w = malloc ((size_t) (n ? n : 1) * sizeof (float));
  if (!x || !u || !w)
  {
    free (x), free (u), free (w);
    return -1;
  }
  BlockShape st = tree_block_shape ();
  double dmax = cell_width (g, tree_get_depth_min ());
  long s = 0;
  foreach_leaf () {
    /* the area relative to the coarsest leaf: in (0, 1], no underflow */
    float wa = g->weight == PINN_WEIGHT_UNIFORM ? 1.f :
      (float) pow (cell_width (g, level) / dmax, dimension);
    const Real * f [PINN_MAX_OUT];
    for (int k = 0; k < nf; k++)
      f [k] = block_data (sblocks, fields [k]);
    block_cells (st, 0) {
      for (int k = 0; k < nf; k++)
        u [s * nf + k] = (float) f [k][c];
      w [s] = wa;
      s++;
    }
  }
  b->x = x, b->u = u, b->w = w;
  b->n = n;
  b->nq = points_per_cell (g);
  return n;
}

/* splitmix64: a stream of 64 bit numbers from a seed */
static unsigned long long splitmix (unsigned long long * s)
{
  unsigned long long z = (*s += 0x9E3779B97F4A7C15ull);
  z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ull;
  z = (z ^ (z >> 27)) * 0x94D049BB133111EBull;
  return z ^ (z >> 31);
}

static double uniform (unsigned long long * s)   /* in [0, 1) */
{
  return (double) (splitmix (s) >> 11) * 0x1.0p-53;
}

int pinn_bamr_collocation (const PinnGrid * g, int per_block,
                           unsigned long seed, PinnBatch * b)
{
  free (b->xc), free (b->wc);
  b->xc = NULL, b->wc = NULL, b->nc = 0;
  if (per_block < 1)
    return 0;
  int nin = pinn_grid_nin (g), m = tree_get_block_m ();
  long nc = (long) tree_get_nleaves () * per_block;
  double * xc = malloc ((size_t) (nc ? nc : 1) * nin * sizeof (double));
  float * wc = malloc ((size_t) (nc ? nc : 1) * sizeof (float));
  if (!xc || !wc)
  {
    free (xc), free (wc);
    return -1;
  }
  unsigned long long rng = seed;
  double bmax = cell_width (g, tree_get_depth_min ()) * m;
  long s = 0;
  foreach_leaf () {
    BLOCK_IJK (bi);
    double side = cell_width (g, level) * m;
    float wa = g->weight == PINN_WEIGHT_UNIFORM ? 1.f :
      (float) pow (side / bmax, dimension);
    for (int r = 0; r < per_block; r++, s++)
    {
      for (int a = 0; a < dimension; a++)
        xc [s * nin + a] = g->x0 [a] + (bi [a] + uniform (&rng)) * side;
      if (g->with_time)
        xc [s * nin + dimension] = g->t;
      wc [s] = wa;
    }
  }
  b->xc = xc, b->wc = wc;
  b->nc = (int) s;
  return (int) s;
}

/* fields = the mean, over the nq points of each cell, of out [n * nq * nf] */
static void scatter (int nq, int nf, const int * fields, const float * out)
{
  BlockShape st = tree_block_shape ();
  long s = 0;
  foreach_leaf () {
    Real * f [PINN_MAX_OUT];
    for (int k = 0; k < nf; k++)
      f [k] = block_data (sblocks, fields [k]);
    block_cells (st, 0) {
      for (int k = 0; k < nf; k++)
      {
        RealSum sum = 0;
        for (int q = 0; q < nq; q++)
          sum += out [(s * nq + q) * nf + k];
        f [k][c] = (Real) (sum / nq);
      }
      s++;
    }
  }
}

static int eval (Pinn * p, const PinnGrid * g, int nf, const int * fields,
                 int residual)
{
  PinnGrid gc = *g;
  if (residual)
    gc.sample = PINN_SAMPLE_CENTRE;
  if (nf < 1 || nf > PINN_MAX_OUT)
    return -1;
  int n = 0, nq = points_per_cell (&gc);
  double * x = all_points (&gc, &n);
  float * out = malloc ((size_t) (n ? n : 1) * nq * nf * sizeof (float));
  int ret = -1;
  if (x && out)
    ret = residual ? pinn_residual (p, n * nq, x, out) :
                     pinn_predict (p, n * nq, x, out);
  if (!ret)
    scatter (nq, nf, fields, out);
  free (x), free (out);
  return ret;
}

int pinn_bamr_predict (Pinn * p, const PinnGrid * g, int nf,
                       const int * fields)
{
  return eval (p, g, nf, fields, 0);
}

int pinn_bamr_residual (Pinn * p, const PinnGrid * g, int nf,
                        const int * fields)
{
  return eval (p, g, nf, fields, 1);
}
