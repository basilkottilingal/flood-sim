/*
.. manifold_advect (): the manifold moved by a face velocity, its volume
.. fractions brought up to date, in 2D and 3D (built with cwrap, linked with
.. libb-amr). The fields are linear, so that their interpolation is exact:
..   one step of a translation and a rotation: every point moved by dt u (x)
..   to round-off; c equal to hm_vof () of the moved mesh, in every cell
..   steps on a tree with level jumps: c equal to hm_vof () after each
..   order: a steady rotation, and one whose speed grows in time, against
..   the exact paths: first order without uf_prev, second with
..   a thread in one cell column (c = 0.4) widened by 0.35 of a cell on both
..   sides: the column full, from its cut neighbours' faces (rounding its
..   old c would say empty)
..   points on the box move along it
..   a step of 0.6 of a cell: refused, nothing changed
..   a cell far from the interface: neither read nor written
..   the edges longer than their cells, counted
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

static const double tol_c = sizeof (Real) == sizeof (float) ? 3e-7 : 1e-12;

/* -------- the velocity: a translation, a rotation about z, a stretch
   along x, all linear, times g -------- */

static double U0, V0, W0, OMEGA, SX, XC = 0.5;

static double vx (double x, double y)
{
  return U0 - OMEGA*(y - 0.5) + SX*(x - XC);
}

static double vy (double x, double y)
{
  (void) y;
  return V0 + OMEGA*(x - 0.5);
}

static void field (double g)
{
  (void) g;
}

static void fill (face_vector uf, double g)
{
  field (g);
#if D == 2
  foreach_face_x ()
    uf.x[] = (Real) (g * vx (x, y));
  foreach_face_y ()
    uf.y[] = (Real) (g * vy (x, y));
#else
  foreach_face_x ()
    uf.x[] = (Real) (g * vx (x, y));
  foreach_face_y ()
    uf.y[] = (Real) (g * vy (x, y));
  foreach_face_z ()
    uf.z[] = (Real) (g * W0);
#endif
}

static void set_field (double u, double v, double w, double omega, double sx)
{
  U0 = u, V0 = v, W0 = w, OMEGA = omega, SX = sx;
}

/* -------- the interfaces -------- */

#if D == 2
static double ball (double x, double y)
{
  return sqrt ((x - 0.5)*(x - 0.5) + (y - 0.47)*(y - 0.47)) - 0.25;
}

static double thread (double x, double y)
{
  (void) y;
  return fabs (x - XC) - 0.2 / 32;
}
#else
static double ball (double x, double y, double z)
{
  return sqrt ((x - 0.5)*(x - 0.5) + (y - 0.47)*(y - 0.47) +
               (z - 0.52)*(z - 0.52)) - 0.25;
}

static double thread (double x, double y, double z)
{
  (void) y, (void) z;
  return fabs (x - XC) - 0.2 / 32;
}
#endif

/* -------- checks -------- */

/* c against hm_vof () of the mesh as it is: the largest difference */
static double against_vof (Hmesh * M, scalar c, scalar c2)
{
  hm_vof (M, c2);
  double e = 0;
  foreach (reduction (max:e))
    e = fmax (e, fabs (c[] - c2[]));
  return e;
}

static double * positions (const Hmesh * M)
{
  double * p = malloc ((size_t) hm_end (M, 0) * 3 * sizeof (double));
  hm_foreach (M, 0, v)
    for (int d = 0; d < D; d++)
      p[3*v + d] = hm_x (M, d)[v];
  return p;
}

static double cell_min (void)
{
  return tree_get_size () / (tree_get_block_m () << tree_get_depth ());
}

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

/* -------- one step, then many -------- */

static void test_steps (void)
{
  tree (3, 4);                         /* level 1: 32 cells a side */
  face_vector uf;
  scalar c = scalar_new (), c2 = scalar_new ();
  set_field (0.3, -0.2, 0.1, 1., 0.2);   /* u_x varies along x too */
  fill (uf, 1.);
  Hmesh * M = hmesh_new (D - 1);
  hm_tree_from_function (M, ball);
  hm_vof (M, c);
  double * p0 = positions (M);
  double dt = 0.3 * cell_min () / 0.8;
  int r = manifold_advect (M, uf, c, dt, uf, 0.);
  double perr = 0;
  hm_foreach (M, 0, v) {
    double x = p0[3*v], y = p0[3*v + 1], want[3] = {
      x + dt*vx (x, y), y + dt*vy (x, y), p0[3*v + 2] + dt*W0 };
    for (int d = 0; d < D; d++)
      perr = fmax (perr, fabs (hm_x (M, d)[v] - want[d]));
  }
  double tol_p = sizeof (Real) == sizeof (float) ? 1e-7 : 1e-14;
  double e = against_vof (M, c, c2);
  CHECK (r == 0 && perr < tol_p && e < tol_c, "one step: %d, points off by "
         "%g, c off hm_vof () by %g", r, perr, e);
  printf ("one step: points off by %.2g, c off hm_vof () by %.2g\n", perr, e);
  free (p0);
  hmesh_free (M);
  tree_destroy ();

  /* level jumps: the lower half in x one level finer */
  tree (3, 4);
  refine_half (2);
  fill (uf, 1.);
  M = hmesh_new (D - 1);
  hm_tree_from_function (M, ball);
  hm_vof (M, c);
  dt = 0.3 * cell_min () / 0.8;
  double emax = 0;
  int steps = D == 2 ? 40 : 6, bad = 0;
  for (int s = 0; s < steps; s++) {
    bad += manifold_advect (M, uf, c, dt, uf, 0.) < 0;
    emax = fmax (emax, against_vof (M, c, c2));
  }
  CHECK (!bad && emax < tol_c, "%d steps, level jumps: %d refused, c off "
         "hm_vof () by %g", steps, bad, emax);
  printf ("%d steps, level jumps: c off hm_vof () by %.2g at most\n", steps,
          emax);
  hmesh_free (M);
  scalar_free (c), scalar_free (c2);
  tree_destroy ();
}

/* -------- order -------- */

/* the error of the points after a quarter turn in n steps, against the
   exact rotation of their starting positions; the speed 1 + growth t */
static int refusals;

static double quarter_turn (int n, int second, double growth)
{
  face_vector uf, ufp;
  scalar c = scalar_new ();
  Hmesh * M = hmesh_new (D - 1);
  double lo[3] = { 0, 0, 0 };
  hm_from_function (M, ball, lo, 1., 16);   /* coarse: the paths only */
  hm_tree_bind (M);
  hm_vof (M, c);
  double * p0 = positions (M);
  double w = 2*M_PI, T = 0.25, dt = T / n, t = 0;
  /* the angle at T: w (T + growth T^2 / 2) */
  double th = w*(T + 0.5*growth*T*T);
  set_field (0, 0, 0, w, 0);
  fill (ufp, 1.);
  for (int s = 0; s < n; s++) {
    fill (uf, 1 + growth*t);
    if (s > 0)
      fill (ufp, 1 + growth*(t - dt));
    refusals += manifold_advect (M, uf, c, dt, ufp, second && s > 0 ? dt :
                                 0.) < 0;
    t += dt;
  }
  /* a first step of order 1 in a second order run: its error is O(dt^2),
     the run's O(dt^2) too */
  double e = 0;
  hm_foreach (M, 0, v) {
    double x = p0[3*v] - 0.5, y = p0[3*v + 1] - 0.5;
    double ex = 0.5 + cos (th)*x - sin (th)*y, ey = 0.5 + sin (th)*x
      + cos (th)*y;
    e = fmax (e, fmax (fabs (hm_x (M, 0)[v] - ex), fabs (hm_x (M, 1)[v] - ey)));
  }
  free (p0);
  hmesh_free (M);
  scalar_free (c);
  return e;
}

static void test_order (void)
{
  tree (3, 4);
  for (int growth = 0; growth < 2; growth++)
    for (int second = 0; second < 2; second++) {
      /* 48 steps: at most 0.4 of a cell at the fastest (1.25 2 pi r) */
      refusals = 0;
      double e1 = quarter_turn (48, second, growth),
        e2 = quarter_turn (96, second, growth);
      double order = log2 (e1 / e2);
      printf ("%s rotation, %s: error %.2g, %.2g (48, 96 steps): order "
              "%.2f\n", growth ? "accelerating" : "steady",
              second ? "second order" : "first order", e1, e2, order);
      CHECK (!refusals, "order: %d steps refused", refusals);
      if (sizeof (Real) == sizeof (double))
        CHECK (second ? order > 1.8 : order > 0.8 && order < 1.3,
               "%s, %s: order %.2f", growth ? "accelerating" : "steady",
               second ? "second" : "first", order);
    }
  tree_destroy ();
}

/* -------- a thread widened from both sides; points on the box -------- */

static void test_thread (void)
{
  tree (3, 4);                         /* 32 cells: XC a cell's centre */
  XC = (16 + 0.5) / 32;
  face_vector uf;
  scalar c = scalar_new (), c2 = scalar_new ();
  Hmesh * M = hmesh_new (D - 1);
  hm_tree_from_function (M, thread);
  hm_vof (M, c);
  double before = 0;
  foreach (reduction (max:before))
    if (fabs (x - XC) < 0.25 / 32)
      before = fmax (before, c[]);
  /* x' = s (x - XC), s dt = 1.75: the sides from 0.2 to 0.55 of a cell */
  double dt = 1;
  set_field (0, 0, 0, 0, 1.75);
  fill (uf, 1.);
  int r = manifold_advect (M, uf, c, dt, uf, 0.);
  double e = against_vof (M, c, c2), col = 1;
  foreach (reduction (min:col))
    if (fabs (x - XC) < 0.25 / 32)
      col = fmin (col, c[]);
  CHECK (r == 0 && fabs (before - 0.4) < 1e-6 && col == 1 && e < tol_c,
         "thread: %d, the column %g before, %g after (full), c off hm_vof () "
         "by %g", r, before, col, e);
  printf ("thread: the column %.2f, widened from both sides: %.2f, c off "
          "hm_vof () by %.2g\n", before, col, e);

  /* along y too, by 0.3 of a cell: the points on the box move along it */
  set_field (0, 0.3 / 32, 0, 0, 0);
  fill (uf, 1.);
  r = manifold_advect (M, uf, c, dt, uf, 0.);
  long on = 0, off = 0;
  hm_foreach (M, 0, v) {
    double y = hm_x (M, 1)[v];
    on += y == 0 || y == 1;
    off += y < 0 || y > 1;
  }
  e = against_vof (M, c, c2);
  CHECK (r == 0 && on > 0 && !off && e < tol_c && hm_check (M) == 0,
         "along y: %d, %ld points on y = 0, 1, %ld off the box, c off by %g",
         r, on, off, e);
  hmesh_free (M);
  scalar_free (c), scalar_free (c2);
  XC = 0.5;
  tree_destroy ();
}

/* -------- refused; untouched; long edges -------- */

static void test_limits (void)
{
  tree (3, 4);
  face_vector uf;
  scalar c = scalar_new ();
  Hmesh * M = hmesh_new (D - 1);
  hm_tree_from_function (M, ball);
  hm_vof (M, c);
  double h = cell_min ();

  /* 0.6 of a cell in one step: refused, nothing changed */
  set_field (1, 0, 0, 0, 0);
  fill (uf, 1.);
  double * p0 = positions (M), sum0 = 0, sum1 = 0;
  foreach (reduction (+:sum0))
    sum0 += c[];
  int r = manifold_advect (M, uf, c, 0.6 * h, uf, 0.);
  long moved = 0;
  hm_foreach (M, 0, v)
    moved += hm_x (M, 0)[v] != p0[3*v];
  foreach (reduction (+:sum1))
    sum1 += c[];
  CHECK (r == -1 && !moved && sum0 == sum1, "CFL 0.6: %d, %ld points moved",
         r, moved);

  /* a cell far from the interface: a value it could not have, kept */
  foreach ()
    if (x < h && y < h)
      c[] = 0.25;
  r = manifold_advect (M, uf, c, 0.3 * h, uf, 0.);
  double corner = 0;
  foreach (reduction (max:corner))
    if (x < h && y < h)
      corner = c[];
  CHECK (r == 0 && corner == (Real) 0.25, "a far cell written: %g", corner);
  foreach ()
    if (x < h && y < h)
      c[] = 0;

  /* stretched along x until edges outgrow their cells: the count */
  set_field (0, 0, 0, 0, 1.);
  fill (uf, 1.);
  r = 0;
  for (int s = 0; s < 400 && r == 0; s++)
    r = manifold_advect (M, uf, c, 0.4 * h, uf, 0.);
  int count = 0;
  hm_foreach (M, 1, e) {
    double a[3], b[3], l = 0;
    for (int d = 0; d < D; d++) {
      a[d] = hm_x (M, d)[hm_vertex (M, 1, e, 0)];
      b[d] = hm_x (M, d)[hm_vertex (M, 1, e, 1)];
      l += (b[d] - a[d])*(b[d] - a[d]);
    }
    count += !(sqrt (l) < fmin (hm_tree_cell_size (M, a),
                                hm_tree_cell_size (M, b)));
  }
  CHECK (r > 0 && r == count, "stretched: %d long edges returned, %d "
         "counted", r, count);
  printf ("stretched: %d edges outgrew their cells\n", r);

  /* refused: a curve in 3D, a surface in 2D */
  Hmesh * W = hmesh_new (D == 2 ? 2 : 1);
  CHECK (manifold_advect (W, uf, c, 0.1 * h, uf, 0.) == -2,
         "a wrong dimension");
  hmesh_free (W);
  free (p0);
  hmesh_free (M);
  scalar_free (c);
  tree_destroy ();
}

int main (void)
{
  test_steps ();
  test_order ();
  test_thread ();
  test_limits ();
  if (fails)
    printf ("%d failed\n", fails);
  else
    printf ("passed, %dD\n", D);
  return fails != 0;
}
