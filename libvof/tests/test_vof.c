/* vof_advection () (vof.h): the plane geometry; exact translation of a slab;
 * c conserved to round-off and kept in [0, 1] through a solenoidal vortex
 * on an adaptive tree, level jumps across the interface, and the shape
 * back after the flow is reversed, closer at a finer grid; the order of
 * the sweeps turned at every call; prolong_fraction: refinement that keeps
 * c, within [0, 1], exact for a plane along an axis, sharper than minmod. */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cwrap/cache.h>
#include "adapt.h"
#include "block-ghosts.h"
#include "vof.h"
#include "vof-ghosts.h"

scalar c, c0;
face_vector uf;

static int checks, fails;
#define CHECK(c) do {                                                      \
    checks++;                                                              \
    if (!(c)) {                                                            \
      fails++;                                                             \
      fprintf (stderr, "FAIL line %d: %s\n", __LINE__, #c);                \
    }                                                                      \
  } while (0)

static const int single = sizeof (Real) < 8;

static double urand (void)
{
  return rand () / (RAND_MAX + 1.);
}

#if dimension > 1
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

#endif

/* refine the leaves coarser than 'to' whose centre is below x = 0.5 (and
 * y = 0.5): level jumps across the middle of the domain */
static void refine_half (int to)
{
  for (int pass = 0; pass < 3 * to; pass++) {
    cache (Index) todo;
    foreach_leaf () {
      double n = 1 << level;
      int in = (i + 0.5) / n < 0.5;
#if dimension > 1
      in = in && (j + 0.5) / n < 0.5 + 0.25;
#endif
      if (level < to && in)
        cache_push (todo, (Index) {IJKL});
    }
    if (!todo.n)
      break;
    foreach (p in todo)
      tree_refine (p);
  }
}

static double finest (void)
{
  return tree_get_size () / (tree_block_shape ().m << tree_get_depth ());
}

static double integral (scalar s)
{
  double sum = 0;
  foreach () {
    double v = 1;
    for (int a = 0; a < dimension; a++)
      v *= Delta;
    sum += s[] * v;
  }
  return sum;
}

/* ---------------------------------------------------------------------- */
/* 1. the geometry                                                          */

static coord C (const double * a)
{
  return (coord) { a [0], a [1], a [2] };
}

static void random_normal (double * n)
{
  for (int a = 0; a < dimension; a++) {
    double r = urand ();
    n [a] = r < 0.1 ? 0 : r < 0.2 ? 1e-9 * (urand () - 0.5)
          : 2 * urand () - 1;
  }
  if (n [0] == 0)
    n [0] = 0.3;
}

/* the fraction of [lo, hi] below the plane, by S^d sample points */
static double sampled (const double * n, double alpha, const double * lo,
                       const double * hi, int S)
{
  long in = 0, all = 0;
  int s [3] = { 0, 0, 0 }, top [3] = { 1, 1, 1 };
  for (int a = 0; a < dimension; a++)
    top [a] = S;
  for (s [2] = 0; s [2] < top [2]; s [2]++)
    for (s [1] = 0; s [1] < top [1]; s [1]++)
      for (s [0] = 0; s [0] < top [0]; s [0]++) {
        double d = 0;
        for (int a = 0; a < dimension; a++)
          d += n [a] * (lo [a] + (s [a] + 0.5) / S * (hi [a] - lo [a]));
        in += d <= alpha;
        all++;
      }
  return (double) in / all;
}

static void test_geometry (void)
{
  double lo [3] = { 0, 0, 0 }, hi [3] = { 1, 1, 1 };
  /* known values */
  double n1 [3] = { 1, 0, 0 };
  CHECK (fabs (vof_plane_volume (C (n1), 0.3, C (lo), C (hi)) - 0.3) < 1e-15);
  CHECK (fabs (vof_plane_alpha (C (n1), 0.3) - 0.3) < 1e-15);
#if dimension > 1
  double n2 [3] = { 1, 1, 0 };
  CHECK (fabs (vof_plane_volume (C (n2), 0.5, C (lo), C (hi)) - 0.125) < 1e-15);
  CHECK (fabs (vof_plane_volume (C (n2), 1.5, C (lo), C (hi)) - 0.875) < 1e-15);
  double n3 [3] = { -1, 1, 0 };                 /* y <= x: half */
  CHECK (fabs (vof_plane_volume (C (n3), 0, C (lo), C (hi)) - 0.5) < 1e-15);
#endif
#if dimension > 2
  double n4 [3] = { 1, 1, 1 };
  CHECK (fabs (vof_plane_volume (C (n4), 0.5, C (lo), C (hi)) - 0.125 / 6) < 1e-15);
  CHECK (fabs (vof_plane_volume (C (n4), 1.5, C (lo), C (hi)) - 0.5) < 1e-15);
  CHECK (fabs (vof_plane_alpha (C (n4), 0.5) - 1.5) < 1e-14);
#endif
  /* alpha and the volume inverse to each other; a box cut in two adds up
   * to the cell; the volume against points sampled */
  double worst = 0, worst_split = 0, worst_sampled = 0;
  for (int t = 0; t < 20000; t++) {
    double n [3], f = urand ();
    if (t % 10 == 0)
      f = t % 20 ? 1e-12 : 1 - 1e-12;
    random_normal (n);
    double alpha = vof_plane_alpha (C (n), f);
    double e = fabs (vof_plane_volume (C (n), alpha, C (lo), C (hi)) - f);
    worst = e > worst ? e : worst;
    double s = urand (), l2 [3] = { s, 0, 0 }, h1 [3] = { s, 1, 1 };
    double split = s * vof_plane_volume (C (n), alpha, C (lo), C (h1))
                 + (1 - s) * vof_plane_volume (C (n), alpha, C (l2), C (hi));
    e = fabs (split - f);
    worst_split = e > worst_split ? e : worst_split;
    if (t < 20) {
      int S = dimension == 3 ? 60 : dimension == 2 ? 400 : 100000;
      e = fabs (sampled (n, alpha, l2, hi, S)
                - vof_plane_volume (C (n), alpha, C (l2), C (hi)));
      worst_sampled = e > worst_sampled ? e : worst_sampled;
    }
  }
  printf ("  geometry: inverse %.2g, split %.2g, sampled %.2g\n",
          worst, worst_split, worst_sampled);
  CHECK (worst < 1e-12);
  CHECK (worst_split < 1e-12);
  CHECK (worst_sampled < (dimension == 3 ? 0.05 : 0.01));
}

/* ---------------------------------------------------------------------- */
/* 2. a slab carried along x by a uniform flow, exact                           */

static double slug (double xl, double xr, double a, double b)
{
  double l = xl > a ? xl : a, r = xr < b ? xr : b;
  return r > l ? (r - l) / (xr - xl) : 0;
}

static void test_translation (void)
{
  int l = dimension == 3 ? 1 : 2;
  tree (l + 1, 4);
  refine_half (l);                              /* level jumps */
  double a = 0.1234, b = 0.3456, dist = 0;
  foreach () {
    c[] = (Real) slug (x - Delta / 2, x + Delta / 2, a, b);
    c0[] = c[];
  }
  foreach_face ()
    uf.x[] = 0;
  foreach_face_x ()
    uf.x[] = 1;
  double dt = 0.4 * finest (), m0 = integral (c);
  int n = (int) (0.4 / dt);
  for (int s = 0; s < n; s++) {
    vof_advection (c, uf, dt);
    dist += dt;
  }
  double e = 0;
  foreach () {
    double d = fabs (c[] - slug (x - Delta / 2, x + Delta / 2, a + dist,
                                  b + dist));
    e = d > e ? d : e;
  }
  double dm = fabs (integral (c) - m0) / m0;
  printf ("  translation across a level jump: error %.2g, mass %.2g\n",
          e, dm);
  CHECK (e < (single ? 1e-4 : 1e-10));
  CHECK (dm < (single ? 1e-5 : 1e-13));
  tree_destroy ();
}

/* ---------------------------------------------------------------------- */
/* 2b. prolong_fraction: refining a plane                                    */

/* the fraction of the cell centred at (x, y, z), size D, below n . x <= a */
static double plane_cell (coord n, double a, double x, double y, double z,
                          double D)
{
  coord lo = { 0, 0, 0 }, hi = { 1, 1, 1 };
  double o = n.x * (x - D / 2);
#if dimension > 1
  o += n.y * (y - D / 2);
#endif
#if dimension > 2
  o += n.z * (z - D / 2);
#endif
  (void) y, (void) z;
  return vof_plane_volume (n, (a - o) / D, lo, hi);
}

/* c (prolong_fraction) and c0 (prolong_minmod) from the plane, refined
 * near the middle: what is kept, the bounds, and the error of each against
 * the plane */
static void refined_plane (coord n, double a, double * ef, double * em)
{
  tree (3, 4);
  foreach () {
    c[] = (Real) plane_cell (n, a, x, y, z, Delta);
    c0[] = c[];
  }
  double m0 = integral (c), m1 = integral (c0);
  ghosts_update (c, c0);        /* refinement reads the parents' ghosts */
  refine_half (2);
  double dm = fabs (integral (c) - m0), dm0 = fabs (integral (c0) - m1);
  int bad = 0;
  *ef = *em = 0;
  foreach () {
    double e = plane_cell (n, a, x, y, z, Delta);
    bad += !(c[] >= 0 && c[] <= 1);
    *ef = fmax (*ef, fabs (c[] - e));
    *em = fmax (*em, fabs (c0[] - e));
  }
  printf ("  refined: c changed by %.2g of %.3g (minmod %.2g)\n", dm, m0, dm0);
  CHECK (dm < (single ? 1e-6 : 1e-13) * m0);
  CHECK (dm0 < (single ? 1e-6 : 1e-13) * m1);
  CHECK (bad == 0);
  tree_destroy ();
}

static void test_prolong_fraction (void)
{
  double ef, em;
  coord along = { 1, 0, 0 };
  refined_plane (along, 0.3 + 1. / 7, &ef, &em);
  printf ("  prolong_fraction, a plane along an axis: error %.2g "
          "(minmod %.2g)\n", ef, em);
  CHECK (ef < (single ? 1e-6 : 1e-12));
#if dimension > 1
  coord oblique = { 1, 0.6, dimension > 2 ? 0.3 : 0 };
  refined_plane (oblique, 0.7, &ef, &em);
  printf ("  prolong_fraction, an oblique plane: error %.2g (minmod %.2g)\n",
          ef, em);
  CHECK (ef < em);
#endif
}

/* ---------------------------------------------------------------------- */
/* 3. a vortex, reversed                                                    */

#if dimension > 1
static double psi (double a, double b)
{
  return sin (M_PI * a) * sin (M_PI * a) * sin (M_PI * b) * sin (M_PI * b)
         / M_PI;
}

/* the vector potential: (psi (y, z), 0, psi (x, y)) in 3D, psi (x, y) in
 * 2D; each face's flux its circulation around the face's edges, so that
 * uf is divergence free face by face, and zero through the walls */
static double A_x (double x, double y, double z)
{
  (void) x;
  return dimension > 2 ? psi (y, z) : 0;
}

static double A_z (double x, double y, double z)
{
  (void) z;
  return psi (x, y);
}

static void set_vortex (double scale)
{
  foreach_face_x () {
    double h = Delta / 2;
    uf.x[] = (Real) (scale * (A_z (x, y + h, z) - A_z (x, y - h, z))
                     / Delta);
  }
  foreach_face_y () {
    double h = Delta / 2;
    uf.y[] = (Real) (scale * ((A_x (x, y, z + h) - A_x (x, y, z - h))
                              - (A_z (x + h, y, z) - A_z (x - h, y, z)))
                     / Delta);
  }
#if dimension > 2
  foreach_face_z () {
    double h = Delta / 2;
    uf.z[] = (Real) (- scale * (A_x (x, y + h, z) - A_x (x, y - h, z))
                     / Delta);
  }
#endif
}

/* the fraction of a sphere (circle) of radius r at p in the cell, by S^d
 * sample points */
static double ball (double x, double y, double z, double D, const double * p,
                    double r)
{
  int S = 8, in = 0, all = 0;
  for (int a = 0; a < S; a++)
    for (int b = 0; b < S; b++)
      for (int e = 0; e < (dimension > 2 ? S : 1); e++) {
        double q [3] = { x + ((a + 0.5) / S - 0.5) * D,
                         y + ((b + 0.5) / S - 0.5) * D,
                         z + ((e + 0.5) / S - 0.5) * D }, d = 0;
        for (int k = 0; k < dimension; k++)
          d += (q [k] - p [k]) * (q [k] - p [k]);
        in += d < r * r;
        all++;
      }
  return (double) in / all;
}

/* the L1 error of the shape after the flow is reversed: the tree uniform
 * at 'level', refined to 'level' + 1 below the middle when 'jumps' */
static double vortex (int level, int jumps, double T)
{
  tree (level + 2, 4);
  refine_all (level);
  if (jumps)
    refine_half (level + 1);
  double p [3] = { 0.5, 0.7, 0.5 }, r = 0.15;
  foreach () {
    c[] = (Real) ball (x, y, z, Delta, p, r);
    c0[] = c[];
  }
  double m0 = integral (c), dt = 0.2 * finest (), t = 0;
  int n = (int) ceil (T / dt), bad = 0;
  dt = T / n;
  double worst = 0;
  for (int s = 0; s < n; s++) {
    set_vortex (cos (M_PI * (t + dt / 2) / T));
    vof_advection (c, uf, dt);
    t += dt;
    double dm = fabs (integral (c) - m0) / m0;
    worst = dm > worst ? dm : worst;
    foreach ()
      bad += !(c[] >= 0 && c[] <= 1);
  }
  double e = 0;
  foreach () {
    double v = 1;
    for (int a = 0; a < dimension; a++)
      v *= Delta;
    e += fabs ((double) c[] - c0[]) * v;
  }
  printf ("  vortex, %s at level %d: %d steps, mass %.2g, L1 %.3g\n",
          jumps ? "level jumps" : "uniform", tree_get_depth (), n, worst, e);
  CHECK (worst < (single ? 1e-4 : 1e-12));
  CHECK (bad == 0);
  /* c's own operators given back after each call */
  CHECK (scalar_get_ops (c).prolongation == prolong_fraction);
  tree_destroy ();
  return e;
}

static void test_vortex (void)
{
  /* tree () starts at level 1 (32^dimension cells); in 3D the finer
   * uniform grid (64^3) is left out, too slow under the sanitizers */
  int l = dimension == 2 ? 2 : 1;
  double T = dimension == 2 ? 2 : 1;
  double coarse = vortex (l, 0, T);
  double jump = vortex (l, 1, T);
  CHECK (coarse < (dimension == 2 ? 2e-3 : 1e-2));
  /* refined where it goes (2D); in 3D the refined half holds only part of
   * the ball's path: the jumps must not make it worse, by much */
  CHECK (jump < (dimension == 2 ? 1 : 1.1) * coarse);
#if dimension == 2
  double fine = vortex (l + 1, 0, T);
  CHECK (fine < coarse / 1.8);              /* at least first order */
#endif
}
#endif

/* ---------------------------------------------------------------------- */
/* 4. the order of the sweeps, turned at every call                          */

#if dimension > 1
/*
.. Problem A: c0 (x) = F (x), u = (gx, gy, gz) (x). Problem B is A turned as
.. vof.c turns the axes (x to y to z; in 2D x and y swapped): B (x, y, z) =
.. A (y, z, x), its y velocity A's x one there, ... (2D: B (x, y) = A (y,
.. x)). A step of B whose sweeps start one axis later is then a step of A
.. turned, bit for bit: B after call n + 1 is A after call n, turned. After
.. call n + dimension, whose order is A's, it is not.
*/
static double F (double x, double y, double z)
{
  double d = (x - 0.45) * (x - 0.45) / 0.04 + (y - 0.55) * (y - 0.55) / 0.02;
  if (dimension > 2)
    d += (z - 0.5) * (z - 0.5) / 0.03;
  return d < 1;
}

static double g (int axis, double x, double y, double z)
{
  switch (axis) {
  case 0: return 0.3 + 0.2 * sin (3 * y + z);
  case 1: return -0.25 + 0.1 * cos (2 * x - z);
  default: return 0.15 * sin (x + 2 * y);
  }
}

/* the fraction of F at the cell centred at (x, y, z), turned if 'turn' */
static double frac (double x, double y, double z, double D, int turn)
{
  int S = 6, in = 0, all = 0;
  double q [3];
  for (int a = 0; a < S; a++)
    for (int b = 0; b < S; b++)
      for (int e = 0; e < (dimension > 2 ? S : 1); e++) {
        q [0] = x + ((a + 0.5) / S - 0.5) * D;
        q [1] = y + ((b + 0.5) / S - 0.5) * D;
        q [2] = dimension > 2 ? z + ((e + 0.5) / S - 0.5) * D : 0;
        if (!turn)
          in += F (q [0], q [1], q [2]) > 0;
        else if (dimension == 2)
          in += F (q [1], q [0], 0) > 0;
        else
          in += F (q [1], q [2], q [0]) > 0;
        all++;
      }
  return (double) in / all;
}

/* A's velocities, or B's */
static double gu (int axis, double x, double y, double z, int turn)
{
  if (!turn)
    return g (axis, x, y, z);
  if (dimension == 2)
    return g (1 - axis, y, x, 0);
  return g ((axis + 2) % 3, y, z, x);
}

static void set_turned (int turn)
{
  foreach ()
    c[] = (Real) frac (x, y, z, Delta, turn);
  foreach_face_x ()
    uf.x[] = (Real) gu (0, x, y, z, turn);
  foreach_face_y ()
    uf.y[] = (Real) gu (1, x, y, z, turn);
#if dimension > 2
  foreach_face_z ()
    uf.z[] = (Real) gu (2, x, y, z, turn);
#endif
}

static void test_rotation (void)
{
  tree (2, 4);                                     /* level 1, uniform */
  int N = tree_block_shape ().m << tree_get_depth ();
  double dt = 0.4 / N;
  Real * a = calloc ((size_t) N * N * (dimension > 2 ? N : 1), sizeof (Real));
#define AT(i, j, k) a [((size_t) (k) * N + (j)) * N + (i)]
  set_turned (0);
  vof_advection (c, uf, dt);                       /* call n             */
  foreach () {
    int i = (int) (x * N), j = (int) (y * N), k = (int) (z * N);
    AT (i, j, k) = c[];
  }
  int same [2] = { 0, 0 }, changed = 0;
  for (int pass = 0; pass < 2; pass++) {
    set_turned (1);
    vof_advection (c, uf, dt);                     /* n + 1, n + dimension */
    foreach () {
      int i = (int) (x * N), j = (int) (y * N), k = (int) (z * N);
      Real v = dimension == 2 ? AT (j, i, 0) : AT (j, k, i);
      same [pass] += c[] == v;
      changed += pass == 0 && c[] != (Real) frac (x, y, z, Delta, 1);
    }
    /* calls to bring the next to n + dimension: a step of nothing */
    for (int s = 0; s < dimension - 2; s++)
      vof_advection (c0, uf, 0);
  }
#undef AT
  int cells = N * N * (dimension > 2 ? N : 1);
  printf ("  rotation: %d of %d cells equal turned, %d with A's order, "
          "%d changed\n", same [0], cells, same [1], changed);
  CHECK (changed > 0);
  CHECK (same [0] == cells);
  CHECK (same [1] < cells);
  free (a);
  tree_destroy ();
}
#endif

int main (void)
{
  srand (1);
  ghost_ops (c, prolongation (prolong_fraction));
  ghost_ops (c0, prolongation (prolong_minmod));
  test_geometry ();
  test_translation ();
  test_prolong_fraction ();
#if dimension > 1
  test_vortex ();
  test_rotation ();
#endif
  printf ("%d checks, %d failed\n", checks, fails);
  return fails != 0;
}
