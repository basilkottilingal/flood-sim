/*
.. hm_from_function (), hm_tree_from_function (): the interface f = 0, in 2D
.. and 3D (built with cwrap, linked with libb-amr):
..   every point on f = 0 (round-off); closed interfaces closed (no boundary
..   half); every normal along grad f (out of f < 0); the Euler
..   characteristic: a circle, two circles (0), a sphere (2), a torus (0),
..   two spheres (4); the enclosed area, volume against the exact, second
..   order as the lattice is refined
..   f = 0 through lattice nodes (a line, a plane): no flat simplex, the
..   exact length, area
..   an interface reaching the box: it ends on the box's faces
..   on the tree: every edge shorter than its cell (checked directly, and
..   by hm_long_edges () with -DEULERIAN_LAGRANGIAN); hm_vof () of the
..   interface against its own enclosed volume, to round-off
..   refused: m not D - 1, n < 1
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "tree.h"
#include "block-ghosts.h"
#include "hmesh-tree.h"
#include "hmesh-vof.h"

static int fails = 0;

#define CHECK(c, ...) do {                                        \
    if (!(c)) {                                                   \
      fprintf (stderr, "%s:%d: ", __FILE__, __LINE__);            \
      fprintf (stderr, __VA_ARGS__);                              \
      fputc ('\n', stderr);                                       \
      fails++;                                                    \
    }                                                             \
  } while (0)

#define D HM_DIMENSION

/* -------- the functions -------- */

static const double R = 0.31, CX = 0.52, CY = 0.47;
#if D == 3
static const double CZ = 0.49;
#endif

#if D == 2
static double circle (double x, double y)
{
  return sqrt ((x - CX)*(x - CX) + (y - CY)*(y - CY)) - R;
}

static double two_circles (double x, double y)
{
  double a = sqrt ((x - 0.28)*(x - 0.28) + (y - 0.3)*(y - 0.3)) - 0.17;
  double b = sqrt ((x - 0.7)*(x - 0.7) + (y - 0.66)*(y - 0.66)) - 0.2;
  return fmin (a, b);
}

static double line_on_nodes (double x, double y)
{
  (void) y;
  return x - 0.5;
}

static double water (double x, double y)
{
  return y - 0.37 - 0.05*sin (2*M_PI*x);
}

static HmFunction fn;
static double F (const double * x)
{
  return fn (x[0], x[1]);
}
#else
static double sphere (double x, double y, double z)
{
  return sqrt ((x - CX)*(x - CX) + (y - CY)*(y - CY) + (z - CZ)*(z - CZ))
    - R;
}

static double torus (double x, double y, double z)
{
  double q = sqrt ((x - 0.5)*(x - 0.5) + (y - 0.5)*(y - 0.5)) - 0.27;
  return sqrt (q*q + (z - 0.5)*(z - 0.5)) - 0.1;
}

static double two_spheres (double x, double y, double z)
{
  double a = sqrt ((x - 0.3)*(x - 0.3) + (y - 0.3)*(y - 0.3) +
                   (z - 0.3)*(z - 0.3)) - 0.18;
  double b = sqrt ((x - 0.7)*(x - 0.7) + (y - 0.68)*(y - 0.68) +
                   (z - 0.66)*(z - 0.66)) - 0.2;
  return fmin (a, b);
}

static double plane_on_nodes (double x, double y, double z)
{
  (void) x, (void) y;
  return z - 0.5;
}

static double water (double x, double y, double z)
{
  return z - 0.37 - 0.05*sin (2*M_PI*x)*cos (2*M_PI*y);
}

static HmFunction fn;
static double F (const double * x)
{
  return fn (x[0], x[1], x[2]);
}
#endif

/* -------- what an interface must be -------- */

typedef struct {
  long   chi;
  double measure, volume;       /* its length (area), what it encloses */
  long   open;                  /* boundary halves */
} Shape;

static Shape check_interface (Hmesh * M, HmFunction f, const char * what)
{
  fn = f;
  int m = M->dim;
  Shape S = { 0, 0, 0, 0 };
  CHECK (hm_check (M) == 0, "%s: invariants", what);
  for (int k = 0; k <= m; k++)
    S.chi += (k & 1 ? -1 : 1) * (long) hm_count (M, k);
  /* the points on f = 0 */
  double fmax_ = 0;
  hm_foreach (M, 0, v) {
    double x[3];
    for (int d = 0; d < D; d++)
      x[d] = hm_x (M, d)[v];
    fmax_ = fmax (fmax_, fabs (F (x)));
  }
  CHECK (fmax_ < 1e-12, "%s: a point at f = %g", what, fmax_);
  /* the normals along grad f; the measure, the volume (divergence: the
     sum of x . A / D) */
  long wrong = 0, flat = 0;
  hm_foreach (M, m, s) {
    double x[3][3] = { { 0 } }, A[3] = { 0, 0, 0 }, c[3] = { 0, 0, 0 };
    for (int i = 0; i <= m; i++)
      for (int d = 0; d < D; d++)
        x[i][d] = hm_x (M, d)[hm_vertex (M, m, s, i)];
#if D == 2
    A[0] = x[1][1] - x[0][1], A[1] = - (x[1][0] - x[0][0]);
#else
    double e1[3], e2[3];
    for (int d = 0; d < 3; d++)
      e1[d] = x[1][d] - x[0][d], e2[d] = x[2][d] - x[0][d];
    A[0] = 0.5*(e1[1]*e2[2] - e1[2]*e2[1]);
    A[1] = 0.5*(e1[2]*e2[0] - e1[0]*e2[2]);
    A[2] = 0.5*(e1[0]*e2[1] - e1[1]*e2[0]);
#endif
    double g[3] = { 0, 0, 0 }, dot = 0, a = 0;
    for (int d = 0; d < D; d++)
      for (int i = 0; i <= m; i++)
        c[d] += x[i][d] / (m + 1);
    for (int d = 0; d < D; d++) {
      double p[3] = { c[0], c[1], c[2] }, q[3] = { c[0], c[1], c[2] };
      p[d] += 1e-6, q[d] -= 1e-6;
      g[d] = (F (p) - F (q)) / 2e-6;
      dot += g[d]*A[d], a += A[d]*A[d];
      S.volume += c[d]*A[d] / D;
    }
    S.measure += sqrt (a);
    wrong += dot <= 0;
    flat += a == 0;
    for (int i = 0; i <= m; i++)
      S.open += hm_half_boundary (M, m, hm_half (m, s, i));
  }
  CHECK (!wrong && !flat, "%s: %ld normals against grad f, %ld flat", what,
         wrong, flat);
  return S;
}

static void test_closed (void)
{
  double lo[3] = { 0, 0, 0 };
#if D == 2
  double want = M_PI*R*R, err[3];
  for (int r = 0; r < 3; r++) {
    Hmesh * M = hmesh_new (1);
    int n = 32 << r;
    long added = hm_from_function (M, circle, lo, 1., n);
    Shape S = check_interface (M, circle, "circle");
    err[r] = fabs (S.volume - want);
    CHECK (added == (long) hm_count (M, 1) && S.chi == 0 && S.open == 0 &&
           err[r] < 2e-3, "circle, n %d: %ld edges, Euler %ld, %ld open, "
           "area off by %g", n, added, S.chi, S.open, err[r]);
    hmesh_free (M);
  }
  printf ("circle: area off by %.2g, %.2g, %.2g (n 32, 64, 128): order "
          "%.2f\n", err[0], err[1], err[2], log2 (err[1]/err[2]));
  CHECK (log2 (err[1]/err[2]) > 1.7, "circle: not second order");

  Hmesh * M = hmesh_new (1);
  hm_from_function (M, two_circles, lo, 1., 64);
  Shape S = check_interface (M, two_circles, "two circles");
  CHECK (S.chi == 0 && S.open == 0 &&
         fabs (S.volume - M_PI*(0.17*0.17 + 0.2*0.2)) < 2e-3,
         "two circles: Euler %ld, %ld open, area %g", S.chi, S.open,
         S.volume);
  hmesh_free (M);
#else
  double want = 4*M_PI*R*R*R/3, err[3];
  for (int r = 0; r < 3; r++) {
    Hmesh * M = hmesh_new (2);
    int n = 16 << r;
    long added = hm_from_function (M, sphere, lo, 1., n);
    Shape S = check_interface (M, sphere, "sphere");
    err[r] = fabs (S.volume - want);
    CHECK (added == (long) hm_count (M, 2) && S.chi == 2 && S.open == 0 &&
           err[r] < 2e-3, "sphere, n %d: %ld triangles, Euler %ld, %ld open, "
           "volume off by %g", n, added, S.chi, S.open, err[r]);
    hmesh_free (M);
  }
  printf ("sphere: volume off by %.2g, %.2g, %.2g (n 16, 32, 64): order "
          "%.2f\n", err[0], err[1], err[2], log2 (err[1]/err[2]));
  CHECK (log2 (err[1]/err[2]) > 1.7, "sphere: not second order");

  Hmesh * M = hmesh_new (2);
  hm_from_function (M, torus, lo, 1., 48);
  Shape S = check_interface (M, torus, "torus");
  double vt = 2*M_PI*M_PI*0.27*0.1*0.1;
  CHECK (S.chi == 0 && S.open == 0 && fabs (S.volume - vt) < 2e-3,
         "torus: Euler %ld, %ld open, volume %g, not %g", S.chi, S.open,
         S.volume, vt);
  hmesh_free (M);

  M = hmesh_new (2);
  hm_from_function (M, two_spheres, lo, 1., 40);
  S = check_interface (M, two_spheres, "two spheres");
  CHECK (S.chi == 4 && S.open == 0, "two spheres: Euler %ld, %ld open",
         S.chi, S.open);
  hmesh_free (M);
#endif
}

/* f = 0 through the lattice's nodes, and an interface reaching the box */
static void test_open (void)
{
  double lo[3] = { 0, 0, 0 };
  Hmesh * M = hmesh_new (D - 1);
#if D == 2
  hm_from_function (M, line_on_nodes, lo, 1., 32);
  Shape S = check_interface (M, line_on_nodes, "on nodes");
#else
  hm_from_function (M, plane_on_nodes, lo, 1., 16);
  Shape S = check_interface (M, plane_on_nodes, "on nodes");
#endif
  CHECK (fabs (S.measure - 1) < 1e-13, "on nodes: measure %.17g, not 1",
         S.measure);
  hm_set_domain (M, lo, 1.);
  long on = 0;
  hm_foreach (M, D - 1, s)
    for (int i = 0; i < D; i++)
      on += hm_half_on_domain (M, D - 1, hm_half (D - 1, s, i)) != 0;
  CHECK (on == S.open && on > 0, "on nodes: %ld of %ld open halves on the "
         "box", on, S.open);
  hmesh_free (M);

  /* the water: its interface ends on the box's sides */
  M = hmesh_new (D - 1);
  hm_from_function (M, water, lo, 1., 40);
  S = check_interface (M, water, "water");
  hm_set_domain (M, lo, 1.);
  on = 0;
  hm_foreach (M, D - 1, s)
    for (int i = 0; i < D; i++)
      on += hm_half_on_domain (M, D - 1, hm_half (D - 1, s, i)) != 0;
  CHECK (on == S.open && on > 0 && S.chi == (D == 2 ? 1 : 1),
         "water: %ld of %ld open halves on the box, Euler %ld", on, S.open,
         S.chi);
  hmesh_free (M);

  /* refused */
  M = hmesh_new (D == 2 ? 2 : 1);
  CHECK (hm_from_function (M, water, lo, 1., 8) == -1, "a wrong dimension");
  hmesh_free (M);
  M = hmesh_new (D - 1);
  CHECK (hm_from_function (M, water, lo, 1., 0) == -1, "n 0");
  hmesh_free (M);
}

/* -------- on the tree -------- */

static void refine_half (int to)
{
  for (int pass = 0; pass < 4 * to; pass++) {
    const TreeCaches * tc = tree_caches ();
    int n = 0;
    Index * todo = malloc ((size_t) tc->nleaves * sizeof (Index));
    for (int q = 0; q < tc->nleaves; q++) {
      Index p = tc->leaves[q].p;
      if (p.level < to && (p.i + 0.5)/(1 << p.level) < 0.5)
        todo[n++] = p;
    }
    for (int a = 0; a < n; a++)
      tree_refine (todo[a]);
    free (todo);
    if (!n)
      break;
  }
}

static void test_tree (void)
{
  tree (3, 4);
  refine_half (2);
  scalar c = scalar_new ();
#if D == 2
  HmFunction f = circle;
#else
  HmFunction f = sphere;
#endif
  Hmesh * M = hmesh_new (D - 1);
  long added = hm_tree_from_function (M, f);
  Shape S = check_interface (M, f, "tree");
  CHECK (added > 0 && S.open == 0 && M->domain, "tree: %ld added, %ld open",
         added, S.open);
  /* every edge shorter than its cells, directly */
  double worst = 0;
  hm_foreach (M, 1, e) {
    double a[3], b[3], l = 0;
    for (int d = 0; d < D; d++) {
      a[d] = hm_x (M, d)[hm_vertex (M, 1, e, 0)];
      b[d] = hm_x (M, d)[hm_vertex (M, 1, e, 1)];
      l += (b[d] - a[d])*(b[d] - a[d]);
    }
    double h = fmin (hm_tree_cell_size (M, a), hm_tree_cell_size (M, b));
    worst = fmax (worst, sqrt (l) / h);
  }
  CHECK (worst < 1, "tree: an edge %g of its cell", worst);
#ifdef EULERIAN_LAGRANGIAN
  CHECK (hm_long_edges (M) == 0, "tree: %u long edges", hm_long_edges (M));
#endif
  /* the volume fractions: the interface's own volume */
  long left = hm_vof (M, c);
  double sum = 0;
  foreach (reduction (+:sum))
    sum += c[] * pow (Delta, D);
  double tol = sizeof (Real) == sizeof (float) ? 1e-6 : 1e-12;
  CHECK (left == 0 && fabs (sum - S.volume) < tol, "tree: vof %ld left, "
         "total %.15g, the interface encloses %.15g", left, sum, S.volume);
  printf ("tree: %ld simplices, the longest edge %.3f of its cell, vof "
          "total off by %.2g\n", added, worst, fabs (sum - S.volume));
  hmesh_free (M);
  scalar_free (c);
  tree_destroy ();
}

int main (void)
{
  test_closed ();
  test_open ();
  test_tree ();
  if (fails)
    printf ("%d failed\n", fails);
  else
    printf ("passed, %dD\n", D);
  return fails != 0;
}
