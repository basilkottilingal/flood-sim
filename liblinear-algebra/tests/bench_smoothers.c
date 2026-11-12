/* The smoothers of multigrid.h against each other: an adaptive tree,
 * div ((1 + x) grad a) - a = b, b from a known field, a from zero to a
 * tolerance; cycles, time, and the mean residual ratio per cycle, for each
 * smoother and number of sweeps. Prints JSON lines. */
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <cwrap/cache.h>
#include "adapt.h"
#include "block-ghosts.h"
#include "multigrid.h"

@boundary bc_odd;

scalar a, b, e, alpha, lambda;

static double now (void)
{
  struct timespec t;
  clock_gettime (CLOCK_MONOTONIC, &t);
  return t.tv_sec + 1e-9 * t.tv_nsec;
}

int main (void)
{
  ghost_ops (a, boundary.all (bc_odd));
  ghost_ops (e, boundary.all (bc_odd));
  int deep = dimension == 3 ? 3 : dimension == 2 ? 5 : 7;
  tree (deep + 1, 4);
  for (int pass = 0; pass < 3 * deep; pass++) {
    cache (Index) todo;
    foreach_leaf () {
      double n = 1 << level, c0 = (i + 0.5) / n - 0.3, r2 = c0 * c0;
#if dimension > 1
      double c1 = (j + 0.5) / n - 0.3;
      r2 += c1 * c1;
#endif
#if dimension > 2
      double c2 = (k + 0.5) / n - 0.3;
      r2 += c2 * c2;
#endif
      if (level < deep && sqrt (r2) < 0.25)
        cache_push (todo, (Index) {IJKL});
    }
    if (!todo.n)
      break;
    foreach (p in todo)
      tree_refine (p);
  }
  foreach () {
    double v = sin (M_PI * x);
#if dimension > 1
    v *= sin (M_PI * y);
#endif
#if dimension > 2
    v *= sin (M_PI * z);
#endif
    e[] = (Real) v;
    alpha[] = (Real) (1 + x);
    lambda[] = -1;
  }
  mg_helmholtz_operator (e, alpha, lambda, b);
  double dmin = tree_get_size () / (tree_block_shape ().m << tree_get_depth ());
  double tol = 100 * REAL_EPSILON * 4 * dimension * 2 / (dmin * dmin);
  int smoothers [3] = { MG_JACOBI, MG_RED_BLACK, MG_GAUSS_SEIDEL };
  const char * names [3] = { "jacobi", "red_black", "gauss_seidel" };
  for (int s = 0; s < 3; s++)
    for (int nr = 1; nr <= 4; nr++) {
      double best = 1e30;
      MGStats st = { 0, 0, 0 };
      for (int rep = 0; rep < 3; rep++) {
        foreach ()
          a[] = 0;
        MGParams p = { .smoother = smoothers [s], .nrelax = nr, .tolerance = tol,
                       .itermax = 200 };
        double t0 = now ();
        st = mg_helmholtz (a, b, alpha, lambda, p);
        best = fmin (best, now () - t0);
      }
      double rate = st.i ? pow (st.resa / st.resb, 1. / st.i) : 0;
      printf ("{\"dimension\": %d, \"smoother\": \"%s\", \"nrelax\": %d, \"cycles\": %d, "
              "\"ms\": %.2f, \"rate\": %.4f, \"converged\": %d, \"leaves\": %d, "
              "\"levels\": [%d, %d]}\n", dimension, names [s], nr, st.i, 1e3 * best, rate,
              st.resa <= tol, tree_get_nleaves (), tree_get_depth_min (), tree_get_depth ());
      fflush (stdout);
    }
  tree_destroy ();
  return 0;
}
