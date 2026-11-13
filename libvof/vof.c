#include <stdio.h>
#include <math.h>
#include <assert.h>

#include "vof.h"
#include "vof-ghosts.h"

/*
.. The plane in the unit cell. With every component of n made positive (the
.. axes of the negative ones flipped) and n scaled to a sum of 1, the fluid
.. is m . x <= a, and its volume V (a) is that of the cube cut by a plane,
.. as Scardovelli and Zaleski (2000); V (1 - a) = 1 - V (a), so only a <=
.. 1/2 is computed, with m sorted, m1 <= m2 (<= m3).
..
.. Written for SIMD lanes and GPU threads: no branch that depends on the
.. data, no loop of a varying count, no array indexed as it runs. The
.. components are sorted by a network of min and max; each case of V is
.. computed and the right one selected; the cubic of 3D is solved by a fixed
.. number of Newton steps. A divisor that may be zero is kept from it
.. (fmax (m, tiny)) in the cases where it is not selected.
..
.. The 3D cases are written as vof_utils.cpp's volume_of_cube (), but
.. without dividing by a small m1 where it can be helped: with b = a - m1,
.. a^3 - b^3 = m1 (a^2 + a b + b^2), so only (a - m2)^3 / m1 and
.. (a - m3)^3 / m1 remain, where a - m2 < m1, a - m3 < m1.
*/

#define VOF_TINY 1e-300

/* a, b, c in increasing order: a network of min and max */
static inline void sort3 (double * a, double * b, double * c)
{
  double lo = fmin (*a, *b), hi = fmax (*a, *b);
  double t = fmax (lo, *c);
  *a = fmin (lo, *c);
  *b = fmin (hi, t);
  *c = fmax (hi, t);
}

static inline void sort2 (double * a, double * b)
{
  double lo = fmin (*a, *b);
  *b = fmax (*a, *b);
  *a = lo;
}

/* The normal made positive, sorted and scaled to a sum of 1: m1 <= m2 <=
.. m3 (those beyond the dimension 0); *shift the sum of its negative
.. components, which the flips move alpha by; the sum. Both sums are taken in
.. sorted order, so that a normal with its components permuted (the same
.. plane, the axes turned) gives the same alpha and volumes, bit for bit. */
typedef struct { double m1, m2, m3, shift, sum, i1, i23; } Unit;

static inline Unit unit_of (coord n)
{
  Unit u;
#if dimension == 1
  u.m1 = u.m2 = 0, u.m3 = fabs (n.x);
  u.shift = fmin (n.x, 0);
  u.sum = u.m3;
#elif dimension == 2
  double a = fabs (n.x), b = fabs (n.y), p = fmin (n.x, 0), q = fmin (n.y, 0);
  sort2 (&a, &b);
  sort2 (&p, &q);
  u.m1 = 0, u.m2 = a, u.m3 = b;
  u.shift = p + q;
  u.sum = a + b;
#else
  double a = fabs (n.x), b = fabs (n.y), c = fabs (n.z);
  double p = fmin (n.x, 0), q = fmin (n.y, 0), r = fmin (n.z, 0);
  sort3 (&a, &b, &c);
  sort3 (&p, &q, &r);
  u.m1 = a, u.m2 = b, u.m3 = c;
  u.shift = (p + q) + r;
  u.sum = (a + b) + c;
#endif
  double inv = 1 / fmax (u.sum, VOF_TINY);
  u.m1 *= inv, u.m2 *= inv, u.m3 *= inv;
  /* the reciprocals the cases divide by, once: 1 / m1, 1 / (2 m2 m3) */
  u.i1 = 1 / fmax (u.m1, VOF_TINY);
  u.i23 = 1 / fmax (2 * u.m2 * u.m3, VOF_TINY);
  return u;
}

#if dimension == 3
/* the volume for m2 <= a < m1 + m2 (two corners or three cut), and dV/da */
static inline double mid_volume (Unit u, double a, double * dv)
{
  double b = a - u.m1, q = a * a + a * b + b * b;
  double e = fmax (a - u.m2, 0), g = fmax (a - u.m3, 0);
  *dv = ((2 * a - u.m1) - (e * e + g * g) * u.i1) * u.i23;
  return (q - (e * e * e + g * g * g) * u.i1) * u.i23 * (1. / 3);
}
#endif

/* V (a) for 0 <= a <= 1/2: each case, the right one selected */
static inline double unit_volume (Unit u, double a)
{
#if dimension == 1
  (void) u;
  return a;
#elif dimension == 2
  double m2 = fmax (u.m2, VOF_TINY);
  double corner = (a / m2) * a / (2 * u.m3), slab = (a - u.m2 / 2) / u.m3;
  return a < u.m2 ? corner : slab;
#else
  double dv, i3 = 1 / u.m3;
  double corner = (a * u.i1) * a * a * u.i23 * (1. / 3);
  double b = a - u.m1;
  double edge = (a * a + a * b + b * b) * u.i23 * (1. / 3);
  double mid = mid_volume (u, a, &dv);
  double slab = (a - (u.m1 + u.m2) / 2) * i3;
  return a < u.m1 ? corner : a < u.m2 ? edge : a < u.m1 + u.m2 ? mid : slab;
#endif
}

/* a, 0 <= a <= 1/2, such that V (a) = v, 0 <= v <= 1/2 */
static inline double unit_alpha (Unit u, double v)
{
#if dimension == 1
  (void) u;
  return v;
#elif dimension == 2
  double m2 = fmax (u.m2, VOF_TINY);
  double corner = sqrt (2 * u.m2 * u.m3 * v), slab = u.m3 * v + u.m2 / 2;
  return v < u.m2 / (2 * u.m3) && m2 > VOF_TINY ? corner : slab;
#else
  double m1 = u.m1, m3 = u.m3, m12 = m1 + u.m2, i3 = 1 / m3;
  double v1 = m1 * m1 * u.i23 * (1. / 3), v2 = v1 + (u.m2 - m1) * 0.5 * i3;
  double corner = cbrt (6 * m1 * u.m2 * m3 * v);
  double edge = m1 / 2 + sqrt (fmax (2 * u.m2 * m3 * v - m1 * m1 / 12, 0));
  double slab = m3 * v + m12 / 2;
  /* the cubic, between m2 and min (m12, 1/2): Newton from the upper end,
   * where V is above v; V is convex there, so the steps come down to the
   * root, and six reach round-off (a test of two million planes) */
  double lo = u.m2, hi = fmin (m12, 0.5), a = hi, dv;
  for (int k = 0; k < 6; k++) {
    double f = mid_volume (u, a, &dv) - v;
    a = fmin (fmax (a - f / fmax (dv, VOF_TINY), lo), hi);
  }
  int in_slab = m12 <= m3 && v >= m12 * 0.5 * i3;
  return v < v1 ? corner : v < v2 ? edge : in_slab ? slab : a;
#endif
}

/* the fraction of the unit cell where n . x <= alpha */
static inline double cell_volume (coord n, double alpha)
{
  Unit u = unit_of (n);
  double a = (alpha - u.shift) / fmax (u.sum, VOF_TINY);
  double h = fmin (fmax (a, 0), 1), low = fmin (h, 1 - h);
  double v = unit_volume (u, low);
  double r = h <= 0.5 ? v : 1 - v;
  r = a <= 0 ? 0 : a >= 1 ? 1 : r;
  return u.sum > 0 ? r : alpha >= 0;
}

double vof_plane_alpha (coord n, double f)
{
  Unit u = unit_of (n);
  f = fmin (fmax (f, 0), 1);
  double low = fmin (f, 1 - f), a = unit_alpha (u, low);
  a = f <= 0.5 ? a : 1 - a;
  return u.sum > 0 ? a * u.sum + u.shift : 0;
}

double vof_plane_volume (coord n, double alpha, coord lo, coord hi)
{
  /* x = lo + (hi - lo) X: the box is the unit cell of X */
  coord nb = { 0, 0, 0 };
  foreach_dimension () {
    nb.x = n.x * (hi.x - lo.x);
    alpha -= n.x * lo.x;
  }
  return cell_volume (nb, alpha);
}

/*
.. The interface on every leaf: n, Youngs' normal (calc_slope () in
.. vof_utils.cpp), minus the gradient of c, and alpha, the plane
.. n . x <= alpha holding c in the unit cell. Two loops: the normal, a
.. stencil with no branch, on every cell (vectorised along rows); then
.. alpha, on the interface cells only, behind one test that is the same for
.. whole regions (the bulk on either side), so that SIMD lanes and GPU warps
.. mostly agree on it, the geometry itself without a branch. Elsewhere alpha
.. is 0: the faces of a bulk cell do not read it.
*/
static void reconstruct (scalar c, vector n, scalar alpha)
{
  foreach ()
    foreach_dimension ()
#if dimension == 1
      n.x[] = (Real) ((double) c[-1] - c[1]);
#elif dimension == 2
      n.x[] = (Real) (((double) c[-1,1] + 2 * c[-1] + c[-1,-1])
                      - ((double) c[1,1] + 2 * c[1] + c[1,-1]));
#else
      n.x[] = (Real) ((4 * (double) c[-1]
                       + 2 * ((double) c[-1,1] + c[-1,-1] + c[-1,0,1] + c[-1,0,-1])
                       + ((double) c[-1,1,1] + c[-1,1,-1] + c[-1,-1,1] + c[-1,-1,-1]))
                      - (4 * (double) c[1]
                         + 2 * ((double) c[1,1] + c[1,-1] + c[1,0,1] + c[1,0,-1])
                         + ((double) c[1,1,1] + c[1,1,-1] + c[1,-1,1] + c[1,-1,-1])));
#endif
  const double eps = VOF_EPS;
  foreach () {
    double f = c[], a = 0;
    if (f > eps && f < 1 - eps) {
      coord m = { 0, 0, 0 };
      foreach_dimension ()
        m.x = n.x[];
      a = fabs (m.x) + fabs (m.y) + fabs (m.z) == 0 ? 0 : vof_plane_alpha (m, f);
    }
    alpha[] = (Real) a;
  }
}

/*
.. The fraction of fluid in the slab of a cell, plane m . x <= a, that
.. crosses its face along x in dt: lambda wide (over the cell's width), at
.. its upper side if up, its lower if not.
*/
foreach_dimension ()
static double slab_x (coord m, double a, double lambda, int up)
{
  if (lambda > 1)
    lambda = 1;
  coord lo = { 0, 0, 0 }, hi = { 1, 1, 1 };
  if (up)
    lo.x = 1 - lambda;
  else
    hi.x = lambda;
  return vof_plane_volume (m, a, lo, hi);
}

/*
.. A sweep along x, written once and turned for y and z:
..   - the volume through each face in dt, over its area, from the upwind
..     cell (-1 if the flow is up x, the cell's own if not), cut by its
..     plane;
..   - what the faces bring in, and cc times what the flow expands.
*/
foreach_dimension ()
static void sweep_x (scalar c, scalar cc, vector n, scalar alpha,
                     face_vector uf, face_vector flux, double dt)
{
  const double eps = VOF_EPS;
  foreach_face_x () {
    double u = uf.x[], phi = 0;
    int up = u > 0;
    double f = up ? (double) c[-1] : (double) c[];
    if (u != 0 && f > eps) {
      if (f >= 1 - eps)
        phi = u * dt;
      else {
        /* the upwind cell's plane: every component, read across x */
        coord m = { 0, 0, 0 };
        m.x = up ? n.x[-1] : n.x[];
#if dimension > 1
        m.y = up ? n.y[-1] : n.y[];
#endif
#if dimension > 2
        m.z = up ? n.z[-1] : n.z[];
#endif
        double a = up ? alpha[-1] : alpha[];
        phi = u * dt * (fabs (m.x) + fabs (m.y) + fabs (m.z) == 0 ? f
                        : slab_x (m, a, fabs (u) * dt / Delta, up));
      }
    }
    flux.x[] = (Real) phi;
  }
  foreach () {
    double v = c[] + ((double) flux.x[] - flux.x[1]
                      + cc[] * ((double) uf.x[1] - uf.x[]) * dt) / Delta;
    c[] = (Real) (v < 0 ? 0 : v > 1 ? 1 : v);
  }
}

/* the calls so far: the first axis of the next */
static unsigned vof_calls;
enum { axis_x, axis_y, axis_z };

void vof_advection (scalar c, face_vector uf, double dt)
{
  static int warned;
  if (!warned && scalar_get_ops (c).prolongation == prolong_linear) {
    warned = 1;
    fprintf (stderr, "vof_advection: %s has the linear prolongation, which "
             "loses c at level jumps when refining: bind prolong_fraction "
             "(vof-ghosts.h)\n", scalar_name (c));
  }
  scalar cc, alpha;
  vector n;
  face_vector flux;
  /*
  .. For the call, c's operators fill n and alpha with c, at level jumps
  .. and on parents (prolong_plane, restriction_plane); theirs do nothing.
  .. So the one ghost update that c, n and alpha get, as the sweep reads
  .. them, makes the three, the ghosts at level jumps cut by the coarse
  .. cell's plane. c's own operators are given back at the end.
  */
  ScalarOps kept = scalar_get_ops (c);
  ScalarOps none = { .prolongation = prolong_none,
                     .restriction = restriction_none };
  foreach_dimension ()
    scalar_set_ops (n.x, none);
  scalar_set_ops (alpha, none);
  ghost_ops (c, prolongation (prolong_plane, n, alpha),
             restriction (restriction_plane, n, alpha));
  foreach ()
    cc[] = c[] > R (0.5);
  int first = (int) (vof_calls++ % dimension);
  for (int k = 0; k < dimension; k++) {
    int axis = (first + k) % dimension;
    reconstruct (c, n, alpha);
    foreach_dimension ()
      if (axis == axis_x)
        sweep_x (c, cc, n, alpha, uf, flux, dt);
  }
  scalar_restore_ops (c, kept);
}
