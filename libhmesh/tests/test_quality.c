/*
.. hm_shape () and hm_quality () (hmesh-quality.h), in the euclidean
.. dimension it is built with, against references computed another way:
..   closed forms: the equilateral and the right isosceles triangle, the
..   regular tetrahedron, the corner tetrahedron (circumcentre (1/2, 1/2,
..   1/2), dihedral angles 90 and acos (1/sqrt 3)), a regular hexagon as a
..   curve (turning 60), a fan of six triangles (valences 6 and 3), a chain
..   on a line (edge ratios at its points)
..   random triangles and tetrahedra: the radius ratio from the
..   circumcentre solved as a linear system, the dihedral angles from the
..   faces' outward normals, the angles of a triangle summing to 180
..   the statistics of a known mesh: min, mean, max and their ids
..   the same shape moved, turned and scaled: the same ratios and angles
..   degenerate simplices counted and left out; refusals
..   -DEULERIAN_LAGRANGIAN: each edge over its cells (a cell size set)
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "hmesh-quality.h"

#define D HM_DIMENSION
#define DEG (180. / 3.14159265358979323846)

static int fails = 0, checks = 0;

#define CHECK(c, ...) do {                                        \
    checks++;                                                     \
    if (!(c)) {                                                   \
      fprintf (stderr, "%s:%d: ", __FILE__, __LINE__);            \
      fprintf (stderr, __VA_ARGS__);                              \
      fputc ('\n', stderr);                                       \
      fails++;                                                    \
    }                                                             \
  } while (0)

/* relatively near; a value not finite only equal to itself */
static inline int near (double a, double b, double tol)
{
  if (!isfinite (a) || !isfinite (b))
    return a == b;
  return fabs (a - b) <= tol * fmax (1., fmax (fabs (a), fabs (b)));
}

#define NEAR(a, b, tol, what) \
  CHECK (near (a, b, tol), "%s: %.17g, expected %.17g", what, \
         (double) (a), (double) (b))

static uint64_t seed = 88172645463325252ull;

static inline double urand (void)
{
  seed ^= seed << 13;
  seed ^= seed >> 7;
  seed ^= seed << 17;
  return (double) (seed >> 11) / 9007199254740992.;
}

/* a point of R^D from three coordinates, the ones past D dropped */
static uint32_t point (Hmesh * M, const double * p)
{
  double x[3] = { p[0], p[1], p[2] };
  return hm_vertex_add (M, x);
}

static uint32_t simplex (Hmesh * M, int k, const double (* p)[3])
{
  uint32_t v[4];
  for (int i = 0; i <= k; i++)
    v[i] = point (M, p[i]);
  return hm_add (M, k, v);
}

static double dot3 (const double * a, const double * b)
{
  return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

static void cross3 (const double * a, const double * b, double * c)
{
  c[0] = a[1]*b[2] - a[2]*b[1];
  c[1] = a[2]*b[0] - a[0]*b[2];
  c[2] = a[0]*b[1] - a[1]*b[0];
}

static double dist3 (const double * a, const double * b)
{
  double d[3] = { a[0] - b[0], a[1] - b[1], a[2] - b[2] };
  return sqrt (dot3 (d, d));
}

/* the solution of the 3x3 system A x = b, by Cramer's rule */
static void solve3 (double A[3][3], const double * b, double * x)
{
  double c0[3], c1[3], c2[3];
  for (int i = 0; i < 3; i++)
    c0[i] = A[i][0], c1[i] = A[i][1], c2[i] = A[i][2];
  double t[3];
  cross3 (c1, c2, t);
  double det = dot3 (c0, t);
  double bc[3];
  cross3 (c1, c2, bc);
  x[0] = dot3 (b, bc) / det;
  cross3 (b, c2, bc);
  x[1] = dot3 (c0, bc) / det;
  cross3 (c1, b, bc);
  x[2] = dot3 (c0, bc) / det;
}

/* --- references, computed another way --- */

/* a triangle's radius ratio: R from the circumcentre in its plane (solved
   with the normal as the third equation), r = 2 A / P */
static inline double ref_tri_radius (const double (* p)[3])
{
  double u[3], w[3], n[3];
  for (int d = 0; d < 3; d++)
    u[d] = p[1][d] - p[0][d], w[d] = p[2][d] - p[0][d];
  cross3 (u, w, n);
  double A[3][3] = { { 2*u[0], 2*u[1], 2*u[2] },
                     { 2*w[0], 2*w[1], 2*w[2] },
                     { n[0], n[1], n[2] } };
  double b[3] = { dot3 (u, u), dot3 (w, w), 0. }, c[3];
  solve3 (A, b, c);
  double R = sqrt (dot3 (c, c));
  double area = 0.5 * sqrt (dot3 (n, n));
  double P = dist3 (p[0], p[1]) + dist3 (p[1], p[2]) + dist3 (p[0], p[2]);
  return R / (2. * (2. * area / P));
}

/* a tetrahedron's radius ratio: R from the circumcentre, a linear system;
   r = 3 V / S */
static inline double ref_tet_radius (const double (* p)[3])
{
  double A[3][3], b[3], c[3], e[3][3];
  for (int i = 0; i < 3; i++) {
    for (int d = 0; d < 3; d++)
      e[i][d] = p[i + 1][d] - p[0][d], A[i][d] = 2. * e[i][d];
    b[i] = dot3 (e[i], e[i]);
  }
  solve3 (A, b, c);
  double R = sqrt (dot3 (c, c)), n[3];
  cross3 (e[1], e[2], n);
  double V = fabs (dot3 (e[0], n)) / 6., S = 0.;
  for (int o = 0; o < 4; o++) {
    const double * a = p[o == 0 ? 1 : 0], * q = p[o <= 1 ? 2 : 1];
    const double * s = p[o <= 2 ? 3 : 2];
    double u[3], w[3], m[3];
    for (int d = 0; d < 3; d++)
      u[d] = q[d] - a[d], w[d] = s[d] - a[d];
    cross3 (u, w, m);
    S += 0.5 * sqrt (dot3 (m, m));
  }
  return R / (3. * (3. * V / S));
}

/* the dihedral angles of a tetrahedron from its faces' outward normals:
   along the edge ij, 180 minus the angle between the normals of the faces
   opposite k and l */
static inline void ref_dihedral (const double (* p)[3], double * mn,
                                 double * mx)
{
  double nrm[4][3];
  for (int o = 0; o < 4; o++) {
    int f[3], n = 0;
    for (int i = 0; i < 4; i++)
      if (i != o)
        f[n++] = i;
    double u[3], w[3], c[3], t[3];
    for (int d = 0; d < 3; d++)
      u[d] = p[f[1]][d] - p[f[0]][d], w[d] = p[f[2]][d] - p[f[0]][d],
        t[d] = p[o][d] - p[f[0]][d];
    cross3 (u, w, c);
    double s = dot3 (c, t) > 0. ? -1. : 1., l = sqrt (dot3 (c, c));
    for (int d = 0; d < 3; d++)
      nrm[o][d] = s * c[d] / l;
  }
  *mn = INFINITY, *mx = -INFINITY;
  for (int k = 0; k < 4; k++)
    for (int l = k + 1; l < 4; l++) {
      double c = fmax (-1., fmin (1., dot3 (nrm[k], nrm[l])));
      double t = 180. - acos (c) * DEG;
      *mn = fmin (*mn, t), *mx = fmax (*mx, t);
    }
}

/* --- the cases --- */

static void refusals (void)
{
  Hmesh * M = hmesh_new (1);
  HmShape q;
  HmQuality Q;
  double p[2][3] = { { 0, 0, 0 }, { 1, 0, 0 } };
  uint32_t e = simplex (M, 1, p);
  CHECK (hm_shape (M, 1, e, &q) == 0 && q.measure == 1. &&
         q.edge_ratio == 1. && isnan (q.angle_min), "an edge's shape");
  CHECK (hm_shape (M, 1, e + 1, &q) == -1, "a dead id refused");
  CHECK (hm_shape (M, 2, 0, &q) == -1, "k > m refused");
  CHECK (hm_shape (M, 0, 0, &q) == -1, "k = 0 refused");
  CHECK (hm_shape (NULL, 1, 0, &q) == -1, "no mesh refused");
  CHECK (hm_quality (NULL, &Q) == -1, "hm_quality (NULL) refused");
  hmesh_free (M);
  /* points only */
  M = hmesh_new (0);
  double x[3] = { 0.5, 0, 0 };
  hm_vertex_add (M, x), hm_vertex_add (M, x);
  CHECK (hm_quality (M, &Q) == 0 && Q.points == 2 && Q.simplices == 2 &&
         Q.edge.n == 0 && Q.valence.n == 0, "points only");
  hmesh_free (M);
}

/* a chain on the x axis, 0 1 3 6: edges 1 2 3; ratios 2 and 1.5 at the
   inner points, no turning; valences 1 2 2 1 */
static void chain (void)
{
  Hmesh * M = hmesh_new (1);
  double x[4] = { 0, 1, 3, 6 };
  uint32_t v[4];
  for (int i = 0; i < 4; i++) {
    double p[3] = { x[i], 0, 0 };
    v[i] = point (M, p);
  }
  for (int i = 0; i < 3; i++) {
    uint32_t e[2] = { v[i], v[i + 1] };
    hm_add (M, 1, e);
  }
  HmQuality Q;
  CHECK (hm_quality (M, &Q) == 0, "chain");
  CHECK (Q.edge.n == 3 && Q.edge.min == 1. && Q.edge.max == 3. &&
         Q.edge.mean == 2., "chain: edges %g %g %g", Q.edge.min,
         Q.edge.mean, Q.edge.max);
  CHECK (Q.measure.n == 3 && Q.measure.mean == 2., "chain: lengths");
  CHECK (Q.edge_ratio.n == 2 && Q.edge_ratio.min == 1.5 &&
         Q.edge_ratio.max == 2. && Q.edge_ratio.argmax == v[1] &&
         Q.edge_ratio.argmin == v[2], "chain: ratios at the points");
  CHECK (Q.turning.n == 2 && Q.turning.max == 0., "chain: straight");
  CHECK (Q.valence.min == 1. && Q.valence.max == 2. &&
         Q.valence.mean == 1.5, "chain: valences");
  CHECK (Q.aspect.n == 0 && Q.degenerate == 0, "chain: no triangles");
  hmesh_free (M);
}

#if D >= 2
/* a regular hexagon as a closed curve: turning 60 at every point */
static void hexagon (void)
{
  Hmesh * M = hmesh_new (1);
  uint32_t v[6];
  for (int i = 0; i < 6; i++) {
    double a = i * 3.14159265358979323846 / 3.;
    double p[3] = { 2. + cos (a), -1. + sin (a), 0 };
    v[i] = point (M, p);
  }
  for (int i = 0; i < 6; i++) {
    uint32_t e[2] = { v[i], v[(i + 1) % 6] };
    hm_add (M, 1, e);
  }
  HmQuality Q;
  hm_quality (M, &Q);
  NEAR (Q.turning.min, 60., 1e-12, "hexagon: turning min");
  NEAR (Q.turning.max, 60., 1e-12, "hexagon: turning max");
  NEAR (Q.edge_ratio.max, 1., 1e-12, "hexagon: edge ratio");
  CHECK (Q.turning.n == 6 && Q.valence.min == 2. && Q.valence.max == 2.,
         "hexagon: every point");
  hmesh_free (M);
}

/* closed forms of two triangles, and the statistics of the mesh of them */
static void triangles (void)
{
  const double s3 = sqrt (3.);
  double eq[3][3] = { { 0, 0, 0 }, { 1, 0, 0 }, { 0.5, s3 / 2., 0 } };
  double ri[3][3] = { { 0, 0, 0 }, { 1, 0, 0 }, { 0, 1, 0 } };
  Hmesh * M = hmesh_new (2);
  uint32_t t0 = simplex (M, 2, eq), t1 = simplex (M, 2, ri);
  HmShape q;
  CHECK (hm_shape (M, 2, t0, &q) == 0, "equilateral");
  NEAR (q.measure, s3 / 4., 1e-15, "equilateral: area");
  NEAR (q.edge_ratio, 1., 1e-15, "equilateral: edge ratio");
  NEAR (q.aspect, 1., 1e-14, "equilateral: aspect");
  NEAR (q.radius_ratio, 1., 1e-14, "equilateral: radius ratio");
  NEAR (q.angle_min, 60., 1e-13, "equilateral: min angle");
  NEAR (q.angle_max, 60., 1e-13, "equilateral: max angle");
  /* right isosceles: aspect sqrt 2 (2 + sqrt 2) / (4 sqrt 3 / 2) =
     (1 + sqrt 2) / sqrt 3; R = sqrt 2 / 2, r = 1 / (2 + sqrt 2) */
  const double asp = (1. + sqrt (2.)) / s3, rad = (1. + sqrt (2.)) / 2.;
  CHECK (hm_shape (M, 2, t1, &q) == 0, "right isosceles");
  NEAR (q.measure, 0.5, 1e-15, "right: area");
  NEAR (q.edge_ratio, sqrt (2.), 1e-15, "right: edge ratio");
  NEAR (q.aspect, asp, 1e-14, "right: aspect");
  NEAR (q.radius_ratio, rad, 1e-14, "right: radius ratio");
  NEAR (q.angle_min, 45., 1e-13, "right: min angle");
  NEAR (q.angle_max, 90., 1e-13, "right: max angle");
  /* obtuse at its third point: 45 at (0,0), atan (1/3) at (4,0), the
     rest at (1,1) */
  {
    Hmesh * O = hmesh_new (2);
    double ob[3][3] = { { 0, 0, 0 }, { 4, 0, 0 }, { 1, 1, 0 } };
    uint32_t to = simplex (O, 2, ob);
    hm_shape (O, 2, to, &q);
    NEAR (q.angle_min, atan (1. / 3.) * DEG, 1e-13, "obtuse: min angle");
    NEAR (q.angle_max, 135. - atan (1. / 3.) * DEG, 1e-13,
          "obtuse: max angle");
    NEAR (q.measure, 2., 1e-15, "obtuse: area");
    hmesh_free (O);
  }
  /* the statistics of the two */
  HmQuality Q;
  hm_quality (M, &Q);
  CHECK (Q.simplices == 2 && Q.aspect.n == 2 && Q.degenerate == 0,
         "two triangles: counts");
  NEAR (Q.aspect.min, 1., 1e-14, "stats: aspect min");
  NEAR (Q.aspect.max, asp, 1e-14, "stats: aspect max");
  NEAR (Q.aspect.mean, (1. + asp) / 2., 1e-14, "stats: aspect mean");
  CHECK (Q.aspect.argmin == t0 && Q.aspect.argmax == t1,
         "stats: aspect's ids");
  NEAR (Q.angle_min.min, 45., 1e-13, "stats: smallest angle");
  NEAR (Q.angle_min.max, 60., 1e-13, "stats: largest smallest angle");
  NEAR (Q.angle_max.max, 90., 1e-13, "stats: largest angle");
  CHECK (Q.angle_min.argmin == t1 && Q.angle_max.argmax == t1,
         "stats: angles' ids");
  NEAR (Q.measure.mean, (s3 / 4. + 0.5) / 2., 1e-15, "stats: mean area");
  CHECK (Q.edge.n == 6, "stats: six edges");
  /* a flat triangle: counted, left out */
  double flat[3][3] = { { 0, 0, 0 }, { 1, 0, 0 }, { 2, 0, 0 } };
  uint32_t t2 = simplex (M, 2, flat);
  CHECK (hm_shape (M, 2, t2, &q) == 0 && isinf (q.aspect) &&
         isinf (q.radius_ratio) && q.measure == 0., "flat: infinite");
  NEAR (q.angle_max, 180., 1e-13, "flat: max angle 180");
  hm_quality (M, &Q);
  CHECK (Q.degenerate == 1 && Q.aspect.n == 2 && Q.measure.n == 3,
         "flat: counted, left out of the ratios");
  NEAR (Q.aspect.max, asp, 1e-14, "flat: ratios unchanged");
  hm_quality_print (stdout, &Q);
  hmesh_free (M);
}

/* a fan of six equilateral triangles about a point: valence 6 there, 3 on
   the rim */
static void fan (void)
{
  Hmesh * M = hmesh_new (2);
  double c[3] = { 0, 0, 0 };
  uint32_t o = point (M, c), r[6];
  for (int i = 0; i < 6; i++) {
    double a = i * 3.14159265358979323846 / 3.;
    double p[3] = { cos (a), sin (a), 0 };
    r[i] = point (M, p);
  }
  for (int i = 0; i < 6; i++) {
    uint32_t t[3] = { o, r[i], r[(i + 1) % 6] };
    hm_add (M, 2, t);
  }
  HmQuality Q;
  hm_quality (M, &Q);
  CHECK (Q.valence.n == 7 && Q.valence.max == 6. && Q.valence.min == 3. &&
         Q.valence.argmax == o, "fan: valences");
  NEAR (Q.valence.mean, 24. / 7., 1e-15, "fan: mean valence");
  NEAR (Q.edge.min, 1., 1e-15, "fan: edges");
  NEAR (Q.edge.max, 1., 1e-15, "fan: edges");
  NEAR (Q.angle_min.min, 60., 1e-12, "fan: angles");
  NEAR (Q.radius_ratio.max, 1., 1e-12, "fan: radius ratio");
  CHECK (Q.edge.n == 12 && Q.simplices == 6, "fan: counts");
  hmesh_free (M);
}

/* random triangles against the circumcentre; the angles sum to 180 */
static void random_triangles (void)
{
  Hmesh * M = hmesh_new (2);
  int bad_r = 0, bad_sum = 0;
  for (int t = 0; t < 1000; t++) {
    double p[3][3] = { { 0 } };
    for (int i = 0; i < 3; i++)
      for (int d = 0; d < D; d++)
        p[i][d] = urand ();
    uint32_t s = simplex (M, 2, p);
    HmShape q;
    hm_shape (M, 2, s, &q);
    if (isinf (q.aspect))
      continue;
    double a0 = 0., a1 = 0., a2 = 0.;
    {
      double u[3], w[3], n[3];
      for (int k = 0; k < 3; k++) {
        const double * o = p[k], * b = p[(k + 1) % 3], * c = p[(k + 2) % 3];
        for (int d = 0; d < 3; d++)
          u[d] = b[d] - o[d], w[d] = c[d] - o[d];
        cross3 (u, w, n);
        double ang = atan2 (sqrt (dot3 (n, n)), dot3 (u, w)) * DEG;
        if (k == 0) a0 = ang;
        else if (k == 1) a1 = ang;
        else a2 = ang;
      }
    }
    if (!near (a0 + a1 + a2, 180., 1e-12) ||
        !near (q.angle_min, fmin (a0, fmin (a1, a2)), 1e-12))
      bad_sum++;
    double ref = ref_tri_radius (p);
    if (q.measure > 1e-3 && !near (q.radius_ratio, ref, 1e-8))
      bad_r++;
  }
  CHECK (!bad_r, "random triangles: %d radius ratios off", bad_r);
  CHECK (!bad_sum, "random triangles: %d angle sums off", bad_sum);
  hmesh_free (M);
}
#endif

#if D == 3
/* the triangle turned, moved and scaled: the same measures */
static void invariance (void)
{
  double ri[3][3] = { { 0, 0, 0 }, { 1, 0, 0 }, { 0, 1, 0 } }, mv[3][3];
  double a = 0.7, b = -1.3, sc = 250., off[3] = { 3, -40, 7 };
  for (int i = 0; i < 3; i++) {
    double x = ri[i][0], y = ri[i][1], z = ri[i][2];
    double x1 = cos (a)*x - sin (a)*y, y1 = sin (a)*x + cos (a)*y, z1 = z;
    double y2 = cos (b)*y1 - sin (b)*z1, z2 = sin (b)*y1 + cos (b)*z1;
    mv[i][0] = sc*x1 + off[0], mv[i][1] = sc*y2 + off[1];
    mv[i][2] = sc*z2 + off[2];
  }
  Hmesh * M = hmesh_new (2);
  uint32_t s0 = simplex (M, 2, ri), s1 = simplex (M, 2, mv);
  HmShape q0, q1;
  hm_shape (M, 2, s0, &q0);
  hm_shape (M, 2, s1, &q1);
  NEAR (q1.aspect, q0.aspect, 1e-12, "moved: aspect");
  NEAR (q1.radius_ratio, q0.radius_ratio, 1e-12, "moved: radius ratio");
  NEAR (q1.angle_min, q0.angle_min, 1e-11, "moved: min angle");
  NEAR (q1.angle_max, q0.angle_max, 1e-11, "moved: max angle");
  NEAR (q1.measure, sc * sc * q0.measure, 1e-12, "moved: area");
  hmesh_free (M);
}

static void tetrahedra (void)
{
  Hmesh * M = hmesh_new (3);
  double reg[4][3] = { { 1, 1, 1 }, { 1, -1, -1 }, { -1, 1, -1 },
                       { -1, -1, 1 } };
  double cor[4][3] = { { 0, 0, 0 }, { 1, 0, 0 }, { 0, 1, 0 }, { 0, 0, 1 } };
  uint32_t t0 = simplex (M, 3, reg), t1 = simplex (M, 3, cor);
  HmShape q;
  hm_shape (M, 3, t0, &q);
  NEAR (q.measure, 8. / 3., 1e-14, "regular: volume");
  NEAR (q.edge_ratio, 1., 1e-15, "regular: edge ratio");
  NEAR (q.aspect, 1., 1e-14, "regular: aspect");
  NEAR (q.radius_ratio, 1., 1e-14, "regular: radius ratio");
  NEAR (q.angle_min, acos (1. / 3.) * DEG, 1e-13, "regular: dihedral");
  NEAR (q.angle_max, acos (1. / 3.) * DEG, 1e-13, "regular: dihedral");
  /* the corner: R = sqrt 3 / 2 (the centre (1/2, 1/2, 1/2)), r = 3 V / S
     = 1/2 / (3/2 + sqrt 3 / 2): both ratios (1 + sqrt 3) / 2 */
  const double cr = (1. + sqrt (3.)) / 2.;
  hm_shape (M, 3, t1, &q);
  NEAR (q.measure, 1. / 6., 1e-15, "corner: volume");
  NEAR (q.edge_ratio, sqrt (2.), 1e-15, "corner: edge ratio");
  NEAR (q.aspect, cr, 1e-14, "corner: aspect");
  NEAR (q.radius_ratio, cr, 1e-14, "corner: radius ratio");
  NEAR (q.angle_min, acos (1. / sqrt (3.)) * DEG, 1e-13,
        "corner: smallest dihedral");
  NEAR (q.angle_max, 90., 1e-13, "corner: largest dihedral");
  /* random ones, against the circumcentre and the outward normals */
  int bad_r = 0, bad_a = 0;
  for (int t = 0; t < 1000; t++) {
    double p[4][3];
    for (int i = 0; i < 4; i++)
      for (int d = 0; d < 3; d++)
        p[i][d] = urand ();
    uint32_t s = simplex (M, 3, p);
    hm_shape (M, 3, s, &q);
    if (q.measure < 1e-3)
      continue;
    if (!near (q.radius_ratio, ref_tet_radius (p), 1e-8))
      bad_r++;
    double mn, mx;
    ref_dihedral (p, &mn, &mx);
    if (!near (q.angle_min, mn, 1e-9) || !near (q.angle_max, mx, 1e-9))
      bad_a++;
  }
  CHECK (!bad_r, "random tetrahedra: %d radius ratios off", bad_r);
  CHECK (!bad_a, "random tetrahedra: %d dihedral angles off", bad_a);
  /* a flat one */
  double flat[4][3] = { { 0, 0, 0 }, { 1, 0, 0 }, { 0, 1, 0 },
                        { 1, 1, 0 } };
  uint32_t t2 = simplex (M, 3, flat);
  hm_shape (M, 3, t2, &q);
  CHECK (isinf (q.aspect) && isinf (q.radius_ratio), "flat tetrahedron");
  HmQuality Q;
  hm_quality (M, &Q);
  CHECK (Q.degenerate == 1 && Q.aspect.n == Q.simplices - 1,
         "tetrahedra: the flat one left out");
  CHECK (Q.aspect.argmin == t0, "tetrahedra: the regular one the best");
  hmesh_free (M);
}
#endif

#ifdef EULERIAN_LAGRANGIAN
static double cell_two (const double * x, void * ctx)
{
  (void) x, (void) ctx;
  return 2.;
}

static void edge_cell (void)
{
  Hmesh * M = hmesh_new (1);
  double p[2][3] = { { 0, 0, 0 }, { 1.5, 0, 0 } };
  simplex (M, 1, p);
  HmQuality Q;
  hm_quality (M, &Q);
  CHECK (Q.edge_cell.n == 0, "no cell size: no edge / cell");
  hm_set_cell_size (M, cell_two, NULL);
  hm_quality (M, &Q);
  CHECK (Q.edge_cell.n == 1 && Q.edge_cell.max == 0.75,
         "edge / cell: %g", Q.edge_cell.max);
  hmesh_free (M);
}
#endif

int main (void)
{
  refusals ();
  chain ();
#if D >= 2
  hexagon ();
  triangles ();
  fan ();
  random_triangles ();
#endif
#if D == 3
  invariance ();
  tetrahedra ();
#endif
#ifdef EULERIAN_LAGRANGIAN
  edge_cell ();
#endif
  if (fails)
    printf ("%d of %d checks failed\n", fails, checks);
  else
    printf ("%d checks passed\n", checks);
  return fails != 0;
}
