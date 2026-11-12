/* the cost of a smoothing sweep's parts on the finest grid of an adaptive
 * tree: one relaxation pass (R), Jacobi's mixing pass (M), a level ghost
 * update (G); the best of several runs, us */
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <cwrap/cache.h>
#include "adapt.h"
#include "block-ghosts.h"

scalar da, dn, res, alpha, lambda;

static double now (void)
{
  struct timespec t;
  clock_gettime (CLOCK_MONOTONIC, &t);
  return t.tv_sec + 1e-9 * t.tv_nsec;
}

#if dimension == 2
#define NB(u, n, d) do { double w_;                                         \
    w_ = 0.5 * ((double) alpha[] + alpha[1]);  n += w_ * u[1];  d += w_;    \
    w_ = 0.5 * ((double) alpha[] + alpha[-1]); n += w_ * u[-1]; d += w_;    \
    w_ = 0.5 * ((double) alpha[] + alpha[0,1]);  n += w_ * u[0,1];  d += w_; \
    w_ = 0.5 * ((double) alpha[] + alpha[0,-1]); n += w_ * u[0,-1]; d += w_; \
  } while (0)
#else
#define NB(u, n, d) do { double w_;                                         \
    w_ = 0.5 * ((double) alpha[] + alpha[1]);  n += w_ * u[1];  d += w_;    \
    w_ = 0.5 * ((double) alpha[] + alpha[-1]); n += w_ * u[-1]; d += w_;    \
    w_ = 0.5 * ((double) alpha[] + alpha[0,1]);  n += w_ * u[0,1];  d += w_; \
    w_ = 0.5 * ((double) alpha[] + alpha[0,-1]); n += w_ * u[0,-1]; d += w_; \
    w_ = 0.5 * ((double) alpha[] + alpha[0,0,1]);  n += w_ * u[0,0,1];  d += w_; \
    w_ = 0.5 * ((double) alpha[] + alpha[0,0,-1]); n += w_ * u[0,0,-1]; d += w_; \
  } while (0)
#endif

static void relax_pass (int l)
{
  foreach_level_or_leaf (l) {
    double n = 0, d = 0, h2 = Delta * Delta;
    NB (da, n, d);
    dn[] = (Real) ((n / h2 - res[]) / (d / h2 - lambda[]));
  }
}

static void mix_pass (int l)
{
  foreach_level_or_leaf (l)
    da[] = (Real) (0.2 * da[] + 0.8 * dn[]);
}

#define BEST(var, reps, body) do {                                       \
    double best_ = 1e30;                                                 \
    for (int r_ = 0; r_ < 7; r_++) {                                     \
      double t0_ = now ();                                               \
      for (int k_ = 0; k_ < (reps); k_++) { body; }                      \
      best_ = fmin (best_, (now () - t0_) / (reps));                     \
    }                                                                    \
    var = 1e6 * best_;                                                   \
  } while (0)

int main (void)
{
  int deep = dimension == 3 ? 3 : 5;
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
    da[] = (Real) (x * x), dn[] = 0, res[] = (Real) y, alpha[] = 1, lambda[] = 0;
  }
  ghosts_update_levels (alpha);
  tree_restrict (res), tree_restrict (lambda);
  int l = tree_get_depth (), cells = 0;
  BlockShape sh = tree_block_shape ();
  foreach_level_or_leaf_block (l) cells += (int) (sh.m * sh.m * (dimension > 2 ? sh.m : 1));
  double R, M, G;
  BEST (R, 20, relax_pass (l));
  BEST (M, 20, mix_pass (l));
  BEST (G, 20, ghosts_update_level (da, l));
  printf ("{\"dimension\": %d, \"levels\": [%d, %d], \"blocks\": %d, \"cells\": %d, "
          "\"R\": %.2f, \"M\": %.2f, \"G\": %.2f}\n", dimension,
          tree_get_depth_min (), l, tree_get_nleaves (), cells, R, M, G);
  tree_destroy ();
  return 0;
}
