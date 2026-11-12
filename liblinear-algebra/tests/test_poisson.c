/* The multigrid solver of multigrid.h: exact for its own discretisation on
 * an adaptive tree (variable alpha, lambda < 0), its rate per cycle, second
 * order for the Laplacian against the continuous solution, and the
 * parameters honoured. Walls: zero value (bc_odd). */
#include <stdio.h>
#include <math.h>
#include <cwrap/cache.h>
#include "adapt.h"
#include "block-ghosts.h"
#include "multigrid.h"

@boundary bc_odd;

scalar a, b, e, alpha, lambda;

static int checks, fails;
#define CHECK(c) do {                                                      \
    checks++;                                                              \
    if (!(c)) {                                                            \
      fails++;                                                             \
      fprintf (stderr, "FAIL line %d: %s\n", __LINE__, #c);                \
    }                                                                      \
  } while (0)

static const double tol = 1e3 * REAL_EPSILON;

/* the solution: zero on the walls of the unit square (cube) */
static double exact (double x, double y, double z)
{
  double v = sin (M_PI * x);
#if dimension > 1
  v *= sin (M_PI * y);
#endif
#if dimension > 2
  v *= sin (M_PI * z);
#endif
  (void) y, (void) z;
  return v;
}

/* refine the leaves coarser than level 'to' within r of the point (0.25,
 * 0.25, 0.25), until none is left (balance refines others) */
static void refine_near (int to, double r)
{
  for (int pass = 0; pass < 2 * to + 2; pass++) {
    cache (Index) todo;
    BlockShape shape = tree_block_shape ();
    foreach_leaf () {
      double n = 1 << level, d = 0;
      double c [3] = { (i + 0.5) / n - 0.25, 0, 0 };
#if dimension > 1
      c [1] = (j + 0.5) / n - 0.25;
#endif
#if dimension > 2
      c [2] = (k + 0.5) / n - 0.25;
#endif
      for (int q = 0; q < dimension; q++)
        d += c [q] * c [q];
      if (level < to && sqrt (d) < r)
        cache_push (todo, (Index) {IJKL});
    }
    (void) shape;
    if (!todo.n)
      break;
    foreach (p in todo)
      tree_refine (p);
  }
}

static void refine_all (int to)
{
  for (int l = 0; l < to; l++) {
    cache (Index) todo;
    foreach_leaf ()
      if (level < to)
        cache_push (todo, (Index) {IJKL});
    foreach (p in todo)
      tree_refine (p);
  }
}

static double max_diff (scalar u, scalar v)
{
  double m = 0;
  foreach () {
    double d = fabs ((double) u[] - v[]);
    m = d > m ? d : m;
  }
  return m;
}

/*
.. The residual's floor: round-off in an operator scaled by 1 / Delta^2,
.. alpha at most amax, |u| at most umax, on the finest leaves
*/
static double floor_of (double amax, double umax)
{
  double dmin = tree_get_size () / (tree_block_shape ().m << tree_get_depth ());
  return REAL_EPSILON * 4 * dimension * amax * umax / (dmin * dmin);
}

/* the tolerance: above the floor, by less in single precision, whose
 * floor is close to the residuals a cycle reaches */
static double tolerance_of (double amax, double umax)
{
  return (sizeof (Real) < 8 ? 2 : 100) * floor_of (amax, umax);
}

/* the Laplacian's error against the continuous solution, uniform at level l */
static double poisson_error (int l)
{
  tree (l + 1, 4);
  refine_all (l);
  foreach () {
    e[] = (Real) exact (x, y, z);
    b[] = (Real) (- dimension * M_PI * M_PI * e[]);
    a[] = 0;
  }
  MGParams p = { .tolerance = tolerance_of (1, 1) };
  MGStats st = mg_poisson (a, b, p);
  double err = max_diff (a, e);
  printf ("  level %d: %d cycles, residual %.3g to %.3g, error %.3g\n",
          l, st.i, st.resb, st.resa, err);
  CHECK (st.resa <= p.tolerance);
  tree_destroy ();
  return err;
}

int main (void)
{
  ghost_ops (a, boundary.all (bc_odd));
  ghost_ops (e, boundary.all (bc_odd));

  /* 1. exact for its own discretisation: b from the operator of e, on an
   * adaptive tree, alpha = 1 + x, lambda = -1; a from zero */
  int deep = dimension == 3 ? 2 : 3;
  tree (deep + 1, 4);
  refine_near (deep, 0.2);
  CHECK (tree_get_depth () == deep && tree_get_depth_min () < deep);
  foreach () {
    e[] = (Real) exact (x, y, z);
    alpha[] = (Real) (1 + x);
    lambda[] = -1;
    a[] = 0;
  }
  mg_helmholtz_operator (e, alpha, lambda, b);
  double floor = floor_of (2, 1);
  /* each smoother, the default (red-black) first */
  const int smoothers [3] = { MG_SMOOTHER_DEFAULT, MG_JACOBI, MG_GAUSS_SEIDEL };
  const char * names [3] = { "red-black", "Jacobi", "Gauss-Seidel" };
  for (int q = 0; q < 3; q++) {
    foreach ()
      a[] = 0;
    MGParams p = { .smoother = smoothers [q], .tolerance = tolerance_of (2, 1) };
    MGStats st = mg_helmholtz (a, b, alpha, lambda, p);
    double err = max_diff (a, e);
    printf ("  adaptive, levels %d to %d, %d leaves, %s: %d cycles, residual "
            "%.3g to %.3g, |a - e| %.3g\n", tree_get_depth_min (),
            tree_get_depth (), tree_get_nleaves (), names [q], st.i, st.resb,
            st.resa, err);
    CHECK (st.resa <= p.tolerance && st.i < 40);
    CHECK (err < 100 * tol);
  }

  /* 2. the rate: the residual divided by about the same each cycle, while
   * well above its floor (in single precision, one or two cycles reach it:
   * nothing measured) */
  {
    foreach ()
      a[] = 0;
    double worst = 0;
    int measured = 0;
    MGParams one = { .itermax = 1, .tolerance = 1e-300 };
    for (int c = 0; c < 6; c++) {
      MGStats s1 = mg_helmholtz (a, b, alpha, lambda, one);
      if (s1.resa < 1e3 * floor)
        break;
      worst = fmax (worst, s1.resa / s1.resb);
      measured++;
    }
    printf ("  rate: the worst of %d cycles %.3f\n", measured, worst);
    CHECK (worst < 0.5);
    CHECK (measured >= 3 || sizeof (Real) < 8);
  }

  /* 3. itermax honoured */
  {
    foreach ()
      a[] = 0;
    MGParams two = { .itermax = 2, .tolerance = 1e-300 };
    CHECK (mg_helmholtz (a, b, alpha, lambda, two).i == 2);
  }
  tree_destroy ();

  /* 4. second order: the Laplacian against the continuous solution,
   * uniform trees one level apart */
  {
    double e0 = poisson_error (1), e1 = poisson_error (2);
    printf ("  order: error ratio %.2f\n", e0 / e1);
    CHECK (e0 / e1 > 3.5);
  }

  printf ("dimension %d: %d checks, %d failed\n", dimension, checks, fails);
  return fails != 0;
}
