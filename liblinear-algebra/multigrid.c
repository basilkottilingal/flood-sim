#include <math.h>
#include <assert.h>

#include "multigrid.h"
#include "block-ghosts.h"

/*
.. The neighbours of the cell, face by face: n += alpha_f u_n and
.. d += alpha_f, alpha_f the mean of alpha in the cell and the neighbour. In
.. a stencil loop: the offsets are cwrap's s[i, j, k].
*/
#define MG_FACE(u, al, n, d, ...) do {                                      \
    double w_ = 0.5 * ((double) al[] + al[__VA_ARGS__]);                     \
    n += w_ * u[__VA_ARGS__];                                                \
    d += w_;                                                                 \
  } while (0)

#if dimension == 1
  #define MG_NEIGHBOURS(u, al, n, d) do {                                   \
      MG_FACE (u, al, n, d, 1); MG_FACE (u, al, n, d, -1);                   \
    } while (0)
#elif dimension == 2
  #define MG_NEIGHBOURS(u, al, n, d) do {                                   \
      MG_FACE (u, al, n, d, 1); MG_FACE (u, al, n, d, -1);                   \
      MG_FACE (u, al, n, d, 0, 1); MG_FACE (u, al, n, d, 0, -1);             \
    } while (0)
#else
  #define MG_NEIGHBOURS(u, al, n, d) do {                                   \
      MG_FACE (u, al, n, d, 1); MG_FACE (u, al, n, d, -1);                   \
      MG_FACE (u, al, n, d, 0, 1); MG_FACE (u, al, n, d, 0, -1);             \
      MG_FACE (u, al, n, d, 0, 0, 1); MG_FACE (u, al, n, d, 0, 0, -1);       \
    } while (0)
#endif

static MGParams defaults (MGParams p)
{
  int m = tree_block_shape ().m;
  assert (p.smoother >= MG_SMOOTHER_DEFAULT && p.smoother <= MG_GAUSS_SEIDEL);
  if (p.smoother == MG_SMOOTHER_DEFAULT)
    p.smoother = MG_RED_BLACK;
  if (p.nrelax <= 0)
    p.nrelax = p.smoother == MG_JACOBI ? 4 : 2;
  /* one sweep a level: the cycle diverges in 2D and 3D (measured), each
     smoother */
  assert (p.nrelax >= 2);
  if (p.nrelax_coarse <= 0)
    p.nrelax_coarse = p.smoother == MG_JACOBI ? m * m : m * m / 2;
  if (p.omega <= 0)
    p.omega = 2. * dimension / (2. * dimension + 1.);
  if (p.tolerance <= 0)
    p.tolerance = 1e-3;
  if (p.itermax <= 0)
    p.itermax = 100;
  if (p.minlevel < 0)
    p.minlevel = 0;
  if (p.minlevel > tree_get_depth ())
    p.minlevel = tree_get_depth ();
  return p;
}

void mg_helmholtz_operator (scalar a, scalar alpha, scalar lambda, scalar out)
{
  assert (out != a);
  foreach () {
    double n = 0, d = 0;
    MG_NEIGHBOURS (a, alpha, n, d);
    out[] = (Real) ((n - d * a[]) / (Delta * Delta) + (double) lambda[] * a[]);
  }
}

/* res = b - the operator of a, on the leaves: the largest |res| */
static double residual (scalar a, scalar b, scalar alpha, scalar lambda,
                        scalar res)
{
  double max = 0;
  foreach () {
    double n = 0, d = 0;
    MG_NEIGHBOURS (a, alpha, n, d);
    double r = b[] - ((n - d * a[]) / (Delta * Delta) + (double) lambda[] * a[]);
    res[] = (Real) r;
    r = r < 0 ? -r : r;
    max = r > max ? r : max;
  }
  return max;
}

/*
.. A sweep of weighted Jacobi for the operator of da = res on the grid at
.. depth l, its ghosts made first: dn the new values, then da between the
.. old and the new.
*/
static void relax (scalar da, scalar res, scalar alpha, scalar lambda,
                   scalar dn, int l, double omega)
{
  foreach_level_or_leaf (l) {
    double n = 0, d = 0, h2 = Delta * Delta;
    MG_NEIGHBOURS (da, alpha, n, d);
    dn[] = (Real) ((n / h2 - res[]) / (d / h2 - lambda[]));
  }
  foreach_level_or_leaf (l)
    da[] = (Real) ((1 - omega) * da[] + omega * dn[]);
}

/*
.. A sweep of Gauss-Seidel on one colour of the grid at depth l (BAMR_RED,
.. BAMR_BLACK), or on every cell in order (BAMR_BOTH): in place, reading
.. the other colour, unchanged in the pass (cwrap checks the offsets).
*/
static void relax_colour (scalar da, scalar res, scalar alpha, scalar lambda,
                          int l, int colour)
{
  foreach_level_or_leaf_colour (l, colour) {
    double n = 0, d = 0, h2 = Delta * Delta;
    MG_NEIGHBOURS (da, alpha, n, d);
    da[] = (Real) ((n / h2 - res[]) / (d / h2 - lambda[]));
  }
}

/* one sweep of the smoother on level l, the ghosts made after each pass */
static void sweep (scalar da, scalar res, scalar alpha, scalar lambda,
                   scalar dn, int l, MGParams p)
{
  switch (p.smoother) {
  case MG_JACOBI:
    relax (da, res, alpha, lambda, dn, l, p.omega);
    ghosts_update_level (da, l);
    break;
  case MG_RED_BLACK:
    relax_colour (da, res, alpha, lambda, l, BAMR_RED);
    ghosts_update_level (da, l);
    relax_colour (da, res, alpha, lambda, l, BAMR_BLACK);
    ghosts_update_level (da, l);
    break;
  default:
    relax_colour (da, res, alpha, lambda, l, BAMR_BOTH);
    ghosts_update_level (da, l);
  }
}

/* one cycle: the residual res on the leaves, a corrected */
static void cycle (scalar a, scalar res, scalar da, scalar dn, scalar alpha,
                   scalar lambda, MGParams p)
{
  tree_restrict (res);                  /* the residual on every level */
  int depth = tree_get_depth ();
  for (int l = p.minlevel; l <= depth; l++) {
    if (l == p.minlevel)
      foreach_level_or_leaf (l)
        da[] = 0;
    else {
      /* level l from level l - 1, by da's prolongation */
      const TreeCaches * tc = tree_caches ();
      for (int q = tc->level_first [l - 1]; q < tc->level_first [l]; q++)
        if (tree_has_children (tc->by_level [q].node))
          block_prolong (tc->by_level [q].p, da);
    }
    ghosts_update_level (da, l);
    int n = l == p.minlevel ? p.nrelax_coarse : p.nrelax;
    for (int k = 0; k < n; k++)
      sweep (da, res, alpha, lambda, dn, l, p);
  }
  foreach ()
    a[] += da[];
}

/* the cycles, dn Jacobi's new values (none for Gauss-Seidel) */
static MGStats solve (scalar a, scalar b, scalar alpha, scalar lambda,
                      scalar res, scalar da, scalar dn, MGParams p)
{
  /* the correction's walls: a's, which must be homogeneous */
  for (int f = 0; f < GHOST_FACES; f++)
    ghost_ops_set (da, GHOST_BOUNDARY, f, scalar_get_boundary (a, f), NULL, 0);
  /* alpha and lambda on every level; alpha's ghosts at every node */
  tree_restrict (lambda);
  ghosts_update_levels (alpha);
  MGStats st = { 0, 0, 0 };
  st.resb = st.resa = residual (a, b, alpha, lambda, res);
  while (st.i < p.itermax && st.resa > p.tolerance) {
    cycle (a, res, da, dn, alpha, lambda, p);
    st.resa = residual (a, b, alpha, lambda, res);
    st.i++;
  }
  return st;
}

MGStats mg_helmholtz (scalar a, scalar b, scalar alpha, scalar lambda,
                      MGParams p)
{
  p = defaults (p);
  scalar res, da;
  if (p.smoother == MG_JACOBI) {
    scalar dn;
    return solve (a, b, alpha, lambda, res, da, dn, p);
  }
  return solve (a, b, alpha, lambda, res, da, -1, p);
}

MGStats mg_poisson (scalar a, scalar b, MGParams p)
{
  scalar alpha, lambda;
  foreach () {
    alpha[] = 1;
    lambda[] = 0;
  }
  return mg_helmholtz (a, b, alpha, lambda, p);
}
