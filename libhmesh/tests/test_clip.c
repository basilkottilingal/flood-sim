/*
.. hm_clip (): manifolds crossing the box, in the dimension it is built
.. with. Against what the cut must give exactly:
..   a straight curve (its length inside), a circle crossing every side of
..   the square (four arcs: eight ends on the box)
..   a rotated triangulation covering the square (area 1, its boundary the
..   square's sides, length 4), one half out of it, a tilted plane in the
..   cube cut by its sides and by z = 0, 1 (3D)
..   a rotated cube of tetrahedra covering the unit cube (volume 1, boundary
..   area 6, 3D)
.. and for each: hm_check (), orientations kept (signed measures positive),
.. Euler characteristic, d d = 0, twins of opposite signs inside, every
.. boundary half on the box where the mesh covered it (hm_half_on_domain
.. ()), fields carried (a linear field on the points interpolated exactly, a
.. constant on the top simplices kept); points near a plane snapped onto it;
.. isolated points and dangling edges; clipping twice changes nothing.
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "hmesh.h"

static int fails = 0;

#define CHECK(c, ...) do {                                        \
    if (!(c)) {                                                   \
      fprintf (stderr, "%s:%d: ", __FILE__, __LINE__);            \
      fprintf (stderr, __VA_ARGS__);                              \
      fputc ('\n', stderr);                                       \
      fails++;                                                    \
    }                                                             \
  } while (0)

static const double O[3] = { 0, 0, 0 };

static uint32_t point (Hmesh * M, double x, double y, double z)
{
  double p[3] = { x, y, z };
  return hm_vertex_add (M, p);
}

static long euler (const Hmesh * M)
{
  long chi = 0;
  for (int k = 0; k <= M->dim; k++)
    chi += (k & 1 ? -1 : 1) * (long) hm_count (M, k);
  return chi;
}

static int boundary_of_boundary (const Hmesh * M)
{
  for (int k = 2; k <= M->dim; k++)
    hm_foreach (M, k, s) {
      uint32_t e[16];
      int sum[16], n = 0;
      for (int i = 0; i <= k; i++) {
        uint32_t h = hm_half (k, s, i), f = hm_half_facet (M, k, h);
        for (int j = 0; j < k; j++) {
          uint32_t g = hm_half (k - 1, f, j), x = hm_half_facet (M, k - 1, g);
          int sg = hm_half_sign (M, k, h)*hm_half_sign (M, k - 1, g), a = 0;
          while (a < n && e[a] != x)
            a++;
          if (a == n) {
            e[n] = x;
            sum[n++] = 0;
          }
          sum[a] += sg;
        }
      }
      for (int a = 0; a < n; a++)
        if (sum[a])
          return -1;
    }
  return 0;
}

/* what a clipped mesh must satisfy whatever it is; its measure, and the
   measures of its boundary on the box and off it */
typedef struct {
  double measure, on, off;
  uint32_t non, noff;
} Summary;

static Summary summary (const Hmesh * M, const char * what, int oriented)
{
  int m = M->dim;
  Summary S = { 0, 0, 0, 0, 0 };
  CHECK (hm_check (M) == 0, "%s: invariants", what);
  CHECK (boundary_of_boundary (M) == 0, "%s: d d is not 0", what);
  int flipped = 0, same = 0;
  hm_foreach (M, m, s) {
    S.measure += hm_measure (M, m, s);
    if (oriented && hm_volume_signed (M, s) <= 0)
      flipped++;
    for (int i = 0; i <= m; i++) {
      uint32_t h = hm_half (m, s, i), t = hm_half_twin (M, m, h);
      double f = m == 1 ? 1. : hm_measure (M, m - 1, hm_half_facet (M, m, h));
      if (t == h) {
        if (hm_half_on_domain (M, m, h))
          S.on += f, S.non++;
        else
          S.off += f, S.noff++;
      }
      else if (hm_half_twin (M, m, t) == h &&
               hm_half_sign (M, m, h) == hm_half_sign (M, m, t))
        same++;
    }
  }
  CHECK (!flipped, "%s: %d simplices of the wrong orientation", what,
         flipped);
  CHECK (!same, "%s: %d twins of the same sign", what, same);
  /* every point in the box */
  double x[3];
  int out = 0;
  hm_foreach (M, 0, v) {
    for (int d = 0; d < HM_DIMENSION; d++)
      x[d] = hm_x (M, d)[v];
    out += !hm_inside (M, x);
  }
  CHECK (!out, "%s: %d points outside", what, out);
  return S;
}

/* the fields: a linear one on the points (their coordinates), a constant on
   the top simplices */
static void fields_set (Hmesh * M, int * fx, int * fc)
{
  *fx = hm_vector_new (M, 0, "x");
  *fc = hm_scalar_new (M, M->dim, "c");
  hm_foreach (M, 0, v)
    for (int d = 0; d < HM_DIMENSION; d++)
      hm_val (M, *fx, d)[v] = (HmReal) hm_x (M, d)[v];
  hm_foreach (M, M->dim, s)
    hm_val (M, *fc, 0)[s] = 3;
}

static void fields_check (const Hmesh * M, int fx, int fc, const char * what)
{
  double tol = sizeof (HmReal) == sizeof (float) ? 1e-5 : 1e-12;
  int bad = 0;
  hm_foreach (M, 0, v)
    for (int d = 0; d < HM_DIMENSION; d++)
      bad += fabs (hm_val (M, fx, d)[v] - hm_x (M, d)[v]) > tol*4;
  hm_foreach (M, M->dim, s)
    bad += hm_val (M, fc, 0)[s] != 3;
  CHECK (!bad, "%s: %d field values not carried", what, bad);
}

/* -------- curves -------- */

static void test_curve (void)
{
  /* a line from out to out: its length inside */
  Hmesh * M = hmesh_new (1);
  int n = 23;
  double a[3] = { -0.4, 0.2, 0.5 }, b[3] = { 1.3, 0.7, 0.5 };
  for (int i = 0; i <= n; i++) {
    double t = (double) i / n;
    point (M, a[0] + t*(b[0] - a[0]), a[1] + t*(b[1] - a[1]),
           a[2] + t*(b[2] - a[2]));
  }
  for (uint32_t i = 0; i < (uint32_t) n; i++)
    hm_add (M, 1, (uint32_t []) { i, i + 1 });
  int fx, fc;
  fields_set (M, &fx, &fc);
  CHECK (hm_clip (M, O, 1.) == 0, "curve: clip");
  Summary S = summary (M, "curve", HM_DIMENSION == 1);
  double slope = HM_DIMENSION == 1 ? 0 : (b[1] - a[1]) / (b[0] - a[0]);
  double want = sqrt (1 + slope*slope);
  CHECK (fabs (S.measure - want) < 1e-12, "curve: length %.15g, not %.15g",
         S.measure, want);
  CHECK (S.non == 2 && S.noff == 0, "curve: %u ends on the box, %u off",
         S.non, S.noff);
  CHECK (euler (M) == 1, "curve: Euler %ld", euler (M));
  fields_check (M, fx, fc, "curve");
  /* the ends: on x = 0 and x = 1 */
  int faces = 0;
  hm_foreach (M, 1, s)
    for (int i = 0; i < 2; i++)
      faces |= hm_half_on_domain (M, 1, hm_half (1, s, i));
  CHECK (faces == 3, "curve: ends on faces %x, not x = 0 and x = 1", faces);
  hmesh_free (M);

#if HM_DIMENSION >= 2
  /* a circle of radius 0.7 about the square's centre: four arcs */
  M = hmesh_new (1);
  n = 200;
  for (int i = 0; i < n; i++) {
    double t = 2*M_PI*(i + 0.5)/n;
    point (M, 0.5 + 0.7*cos (t), 0.5 + 0.7*sin (t), 0.5);
  }
  for (uint32_t i = 0; i < (uint32_t) n; i++)
    hm_add (M, 1, (uint32_t []) { i, (i + 1) % (uint32_t) n });
  CHECK (hm_clip (M, O, 1.) == 0, "circle: clip");
  S = summary (M, "circle", 0);
  CHECK (S.non == 8 && S.noff == 0, "circle: %u ends on the box, %u off",
         S.non, S.noff);
  CHECK (euler (M) == 4, "circle: Euler %ld, not 4 arcs", euler (M));
  hmesh_free (M);
#endif
}

/* -------- surfaces -------- */

#if HM_DIMENSION >= 2
/* an n x n triangulation of [lo, hi]^2, rotated by angle about the
   square's centre, lifted to z = z0 + zx x */
static Hmesh * grid (int n, double lo, double hi, double angle, double z0,
                     double zx)
{
  Hmesh * M = hmesh_new (2);
  double c = cos (angle), s = sin (angle);
  for (int j = 0; j <= n; j++)
    for (int i = 0; i <= n; i++) {
      double u = lo + (hi - lo)*i/n - 0.5, v = lo + (hi - lo)*j/n - 0.5;
      double x = 0.5 + c*u - s*v, y = 0.5 + s*u + c*v;
      point (M, x, y, z0 + zx*x);
    }
  uint32_t w = (uint32_t) n + 1;
  for (uint32_t j = 0; j < (uint32_t) n; j++)
    for (uint32_t i = 0; i < (uint32_t) n; i++) {
      uint32_t a = j*w + i, b = a + 1, d = a + w + 1, e = a + w;
      hm_add (M, 2, (uint32_t []) { a, b, d });
      hm_add (M, 2, (uint32_t []) { a, d, e });
    }
  return M;
}
#endif

static void test_surface (void)
{
#if HM_DIMENSION >= 2
  /* covering the square, rotated: area 1, the boundary the square's */
  Hmesh * M = grid (17, -0.6, 1.6, 0.5236, 0.5, 0);
  int fx, fc;
  fields_set (M, &fx, &fc);
  CHECK (hm_clip (M, O, 1.) == 0, "covering: clip");
  Summary S = summary (M, "covering", HM_DIMENSION == 2);
  CHECK (fabs (S.measure - 1) < 1e-12, "covering: area %.15g", S.measure);
  CHECK (fabs (S.on - 4) < 1e-12 && S.noff == 0, "covering: boundary %.15g "
         "on the box, %u halves off", S.on, S.noff);
  CHECK (euler (M) == 1, "covering: Euler %ld", euler (M));
  fields_check (M, fx, fc, "covering");
  uint32_t counts[3] = { hm_count (M, 0), hm_count (M, 1), hm_count (M, 2) };

  /* twice: nothing changes */
  CHECK (hm_clip (M, O, 1.) == 0, "twice: clip");
  CHECK (hm_count (M, 0) == counts[0] && hm_count (M, 1) == counts[1] &&
         hm_count (M, 2) == counts[2], "twice: the counts changed");
  /* the box holds now */
  CHECK (point (M, 1.5, 0.5, 0.5) == HM_NONE, "a point outside, after");
  hmesh_free (M);

  /* half out: [0.5, 1.5]^2, area 1/4, half its boundary on the box */
  M = grid (8, 0.5, 1.5, 0, 0.5, 0);
  S = summary (M, "half: before", HM_DIMENSION == 2);
  CHECK (hm_clip (M, O, 1.) == 0, "half: clip");
  S = summary (M, "half", HM_DIMENSION == 2);
  CHECK (fabs (S.measure - 0.25) < 1e-12, "half: area %.15g", S.measure);
  CHECK (fabs (S.on - 1) < 1e-12 && fabs (S.off - 1) < 1e-12,
         "half: boundary %.15g on the box, %.15g off", S.on, S.off);
  hmesh_free (M);

  /* a point 1e-13 inside x = 0: snapped onto it, no sliver */
  M = hmesh_new (2);
  uint32_t a = point (M, 1e-13, 0.5, 0.5), b = point (M, 0.5, 0.2, 0.5),
    c = point (M, 0.5, 0.8, 0.5), d = point (M, -0.5, 0.5, 0.5);
  hm_add (M, 2, (uint32_t []) { b, c, a });
  hm_add (M, 2, (uint32_t []) { a, c, d });
  CHECK (hm_clip (M, O, 1.) == 0, "snap: clip");
  summary (M, "snap", HM_DIMENSION == 2);
  /* (a, c, d) leaves one triangle, a on x = 0: not a quadrangle in two,
     one a sliver 1e-13 wide */
  double amin = 1;
  hm_foreach (M, 2, t)
    amin = fmin (amin, hm_measure (M, 2, t));
  CHECK (hm_x (M, 0)[a] == 0 && hm_count (M, 2) == 2 && amin > 1e-3,
         "snap: x = %g, %u triangles, the smallest %g", hm_x (M, 0)[a],
         hm_count (M, 2), amin);

  /* isolated points and a dangling edge */
  hmesh_free (M);
  M = grid (4, 0.2, 0.8, 0, 0.5, 0);
  uint32_t pin = point (M, 0.1, 0.1, 0.5), pout = point (M, 1.1, 0.1, 0.5);
  uint32_t e0 = point (M, 0.9, 0.9, 0.5), e1 = point (M, 1.4, 1.2, 0.5);
  hm_add (M, 1, (uint32_t []) { e0, e1 });
  CHECK (hm_clip (M, O, 1.) == 0, "isolated: clip");
  summary (M, "isolated", HM_DIMENSION == 2);
  CHECK (hm_alive (M, 0, pin) && !hm_alive (M, 0, pout) &&
         hm_alive (M, 0, e0) && !hm_alive (M, 0, e1),
         "isolated: the points in kept, out gone");
  uint32_t dangling = 0, ends = 0;
  hm_foreach (M, 1, e)
    if (hm_use (M, 1, e) == HM_NONE) {
      dangling++;
      ends += hm_half_on_domain (M, 1, hm_half (1, e, 0)) != 0 ||
        hm_half_on_domain (M, 1, hm_half (1, e, 1)) != 0;
    }
  CHECK (dangling == 1 && ends == 1, "isolated: %u dangling edges, %u ending "
         "on the box", dangling, ends);
  hmesh_free (M);
#endif

#if HM_DIMENSION == 3
  /* the plane z = 1.5 x - 0.25 over y: cut by x, y and z = 0, 1;
     x in [1/6, 5/6] */
  M = grid (12, -0.5, 1.5, 0.3, -0.25, 1.5);
  fields_set (M, &fx, &fc);
  CHECK (hm_clip (M, O, 1.) == 0, "tilted: clip");
  S = summary (M, "tilted", 0);
  double slant = sqrt (1 + 2.25), area = 2./3*slant;
  CHECK (fabs (S.measure - area) < 1e-12, "tilted: area %.15g, not %.15g",
         S.measure, area);
  CHECK (fabs (S.on - (2 + 2*area)) < 1e-12 && S.noff == 0,
         "tilted: boundary %.15g, not %.15g", S.on, 2 + 2*area);
  CHECK (euler (M) == 1, "tilted: Euler %ld", euler (M));
  fields_check (M, fx, fc, "tilted");
  hmesh_free (M);
#endif
}

/* -------- volumes -------- */

#if HM_DIMENSION == 3
static void rotate (double * p, int a, int b, double angle)
{
  double c = cos (angle), s = sin (angle);
  double u = p[a] - 0.5, v = p[b] - 0.5;
  p[a] = 0.5 + c*u - s*v;
  p[b] = 0.5 + s*u + c*v;
}

/* the Kuhn tetrahedra of [-1, 2]^3 in n^3 cubes, rotated about z then x */
static Hmesh * cube (int n)
{
  Hmesh * M = hmesh_new (3);
  uint32_t w = (uint32_t) n + 1;
  for (int k = 0; k <= n; k++)
    for (int j = 0; j <= n; j++)
      for (int i = 0; i <= n; i++) {
        double p[3] = { -1 + 3.*i/n, -1 + 3.*j/n, -1 + 3.*k/n };
        rotate (p, 0, 1, 0.35);
        rotate (p, 1, 2, 0.2);
        point (M, p[0], p[1], p[2]);
      }
  static const int perm[6][3] = {{0,1,2},{0,2,1},{1,0,2},{1,2,0},{2,0,1},
                                 {2,1,0}};
  for (uint32_t k = 0; k < (uint32_t) n; k++)
    for (uint32_t j = 0; j < (uint32_t) n; j++)
      for (uint32_t i = 0; i < (uint32_t) n; i++)
        for (int p = 0; p < 6; p++) {
          uint32_t c[3] = { i, j, k }, v[4];
          v[0] = (k*w + j)*w + i;
          for (int q = 0; q < 3; q++) {
            c[perm[p][q]]++;
            v[q + 1] = (c[2]*w + c[1])*w + c[0];
          }
          uint32_t t = hm_add (M, 3, v);
          if (hm_volume_signed (M, t) < 0) {
            hm_remove (M, 3, t);
            uint32_t x = v[0];
            v[0] = v[1], v[1] = x;
            hm_add (M, 3, v);
          }
        }
  return M;
}
#endif

static void test_volume (void)
{
#if HM_DIMENSION == 3
  Hmesh * M = cube (7);
  int fx, fc;
  fields_set (M, &fx, &fc);
  CHECK (hm_clip (M, O, 1.) == 0, "cube: clip");
  Summary S = summary (M, "cube", 1);
  CHECK (fabs (S.measure - 1) < 1e-12, "cube: volume %.15g", S.measure);
  CHECK (fabs (S.on - 6) < 1e-12 && S.noff == 0, "cube: boundary %.15g on "
         "the box, %u halves off", S.on, S.noff);
  CHECK (euler (M) == 1, "cube: Euler %ld", euler (M));
  fields_check (M, fx, fc, "cube");
  /* the faces of the box the boundary halves are on: all six */
  int faces = 0;
  hm_foreach (M, 3, s)
    for (int i = 0; i < 4; i++)
      faces |= hm_half_on_domain (M, 3, hm_half (3, s, i));
  CHECK (faces == 63, "cube: boundary on faces %x", faces);
  printf ("cube: %u tetrahedra inside\n", hm_count (M, 3));
  hmesh_free (M);
#endif
}

int main (void)
{
  test_curve ();
  test_surface ();
  test_volume ();
  if (fails)
    printf ("%d failed\n", fails);
  else
    printf ("passed, %dD\n", HM_DIMENSION);
  return fails != 0;
}
