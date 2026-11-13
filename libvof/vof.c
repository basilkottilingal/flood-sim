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
.. The 3D cases are written as vof_utils.cpp's volume_of_cube (), but
.. without dividing by a small m1 where it can be helped: with b = a - m1,
.. a^3 - b^3 = m1 (a^2 + a b + b^2), so only (a - m2)^3 / m1 and
.. (a - m3)^3 / m1 remain, where a - m2 < m1, a - m3 < m1.
*/

static void sort_up (double * m, int n)
{
  for (int i = 1; i < n; i++)
    for (int j = i; j > 0 && m [j] < m [j - 1]; j--) {
      double t = m [j];
      m [j] = m [j - 1];
      m [j - 1] = t;
    }
}

#if dimension == 3
/* V (a) and dV/da for a <= 1/2, m sorted, sum 1 */
static double unit_volume3 (const double * m, double a, double * dv)
{
  double m1 = m [0], m2 = m [1], m3 = m [2], m12 = m1 + m2;
  if (a <= 0) {
    if (dv)
      *dv = 0;
    return 0;
  }
  if (a < m1) {                              /* a corner                    */
    double p = (a / m1) * (a / m2);
    if (dv)
      *dv = p / (2 * m3);
    return p * a / (6 * m3);
  }
  double b = a - m1, q = a * a + a * b + b * b;
  if (a < m2) {                              /* an edge                     */
    if (dv)
      *dv = (2 * a - m1) / (2 * m2 * m3);
    return q / (6 * m2 * m3);
  }
  if (a < m12) {
    double e = a - m2, ee = e * e * (e / m1), dd = e * (e / m1);
    if (a >= m3) {                           /* past the third corner too   */
      double g = a - m3;
      ee += g * g * (g / m1);
      dd += g * (g / m1);
    }
    if (dv)
      *dv = ((2 * a - m1) - dd) / (2 * m2 * m3);
    return (q - ee) / (6 * m2 * m3);
  }
  if (dv)                                    /* a slab across m3            */
    *dv = 1 / m3;
  return (a - m12 / 2) / m3;
}
#endif

/* V (a) for a <= 1/2, m sorted, sum 1 */
static double unit_volume (const double * m, double a)
{
  if (a <= 0)
    return 0;
#if dimension == 1
  (void) m;
  return a;
#elif dimension == 2
  if (a < m [0])
    return (a / m [0]) * a / (2 * m [1]);
  return (a - m [0] / 2) / m [1];
#else
  return unit_volume3 (m, a, 0);
#endif
}

/* a <= 1/2 such that V (a) = v, v <= 1/2, m sorted, sum 1 */
static double unit_alpha (const double * m, double v)
{
  if (v <= 0)
    return 0;
#if dimension == 1
  (void) m;
  return v;
#elif dimension == 2
  double m1 = m [0], m2 = m [1];
  if (v < m1 / (2 * m2))
    return sqrt (2 * m1 * m2 * v);
  return m2 * v + m1 / 2;
#else
  double m1 = m [0], m2 = m [1], m3 = m [2], m12 = m1 + m2;
  double v1 = m1 > 0 ? m1 * (m1 / (6 * m2 * m3)) : 0;
  if (v < v1)
    return cbrt (6 * m1 * m2 * m3 * v);
  if (v < v1 + (m2 - m1) / (2 * m3)) {
    double d = 2 * m2 * m3 * v - m1 * m1 / 12;
    return m1 / 2 + sqrt (d > 0 ? d : 0);
  }
  if (m12 <= m3 && v >= m12 / (2 * m3))
    return m3 * v + m12 / 2;
  /* a cubic, between m2 and min (m12, 1/2): Newton, bisection where it
     would leave the bracket (V is increasing) */
  double lo = m2, hi = m12 < 0.5 ? m12 : 0.5, a = 0.5 * (lo + hi);
  for (int it = 0; it < 100; it++) {
    double dv, f = unit_volume3 (m, a, &dv) - v;
    if (f > 0)
      hi = a;
    else
      lo = a;
    double next = dv > 0 ? a - f / dv : 0.5 * (lo + hi);
    if (!(next > lo && next < hi))
      next = 0.5 * (lo + hi);
    if (fabs (next - a) <= 4 * DBL_EPSILON * a || hi - lo <= 4 * DBL_EPSILON)
      return next;
    a = next;
  }
  return a;
#endif
}

/*
.. m = |n| / sum, sorted; shift = the sum of the negative components; the
.. sum. Both sums are taken in sorted order, so that a normal with its
.. components permuted (the same plane, the axes turned) gives the same
.. alpha and volumes, bit for bit.
*/
static double normalise (const double * n, double * m, double * shift)
{
  double neg [3];
  int k = 0;
  for (int i = 0; i < dimension; i++) {
    m [i] = fabs (n [i]);
    if (n [i] < 0)
      neg [k++] = n [i];
  }
  sort_up (m, dimension);
  sort_up (neg, k);
  double s = 0;
  for (int i = 0; i < dimension; i++)
    s += m [i];
  *shift = 0;
  for (int i = 0; i < k; i++)
    *shift += neg [i];
  if (s > 0)
    for (int i = 0; i < dimension; i++)
      m [i] /= s;
  return s;
}

/* the fraction of the unit cell where n . x <= alpha */
static double cell_volume (const double * n, double alpha)
{
  double m [3], shift, s = normalise (n, m, &shift);
  if (s == 0)
    return alpha >= 0;
  double a = (alpha - shift) / s;
  if (a <= 0)
    return 0;
  if (a >= 1)
    return 1;
  return a <= 0.5 ? unit_volume (m, a) : 1 - unit_volume (m, 1 - a);
}

static void components (coord c, double * a)
{
  a [0] = c.x;
#if dimension > 1
  a [1] = c.y;
#endif
#if dimension > 2
  a [2] = c.z;
#endif
}

double vof_plane_alpha (coord n, double f)
{
  double nn [3], m [3], shift;
  components (n, nn);
  double s = normalise (nn, m, &shift);
  if (s == 0)
    return 0;
  f = f < 0 ? 0 : f > 1 ? 1 : f;
  double a = f <= 0.5 ? unit_alpha (m, f) : 1 - unit_alpha (m, 1 - f);
  return a * s + shift;
}

double vof_plane_volume (coord n, double alpha, coord lo, coord hi)
{
  /* x = lo + (hi - lo) X: the box is the unit cell of X */
  double nn [3], l [3], h [3], nb [3];
  components (n, nn), components (lo, l), components (hi, h);
  for (int i = 0; i < dimension; i++) {
    nb [i] = nn [i] * (h [i] - l [i]);
    alpha -= nn [i] * l [i];
  }
  return cell_volume (nb, alpha);
}

/*
.. The interface on every leaf: n, Youngs' normal (calc_slope () in
.. vof_utils.cpp), minus the gradient of c, and alpha, the plane
.. n . x <= alpha holding c in the unit cell. All the branching of the
.. reconstruction is here, once per cell; the sweep's faces only cut the
.. plane of their upwind cell.
*/
static void reconstruct (scalar c, vector n, scalar alpha)
{
  foreach () {
    coord m = { 0, 0, 0 };
    foreach_dimension ()
#if dimension == 1
      m.x = (double) c[-1] - c[1];
#elif dimension == 2
      m.x = ((double) c[-1,1] + 2 * c[-1] + c[-1,-1])
          - ((double) c[1,1] + 2 * c[1] + c[1,-1]);
#else
      m.x = (4 * (double) c[-1]
             + 2 * ((double) c[-1,1] + c[-1,-1] + c[-1,0,1] + c[-1,0,-1])
             + ((double) c[-1,1,1] + c[-1,1,-1] + c[-1,-1,1] + c[-1,-1,-1]))
          - (4 * (double) c[1]
             + 2 * ((double) c[1,1] + c[1,-1] + c[1,0,1] + c[1,0,-1])
             + ((double) c[1,1,1] + c[1,1,-1] + c[1,-1,1] + c[1,-1,-1]));
#endif
    foreach_dimension ()
      n.x[] = (Real) m.x;
    alpha[] = (Real) (fabs (m.x) + fabs (m.y) + fabs (m.z) == 0 ? 0
                      : vof_plane_alpha (m, c[]));
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
