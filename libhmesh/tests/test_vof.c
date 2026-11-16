/*
.. hm_vof (): volume fractions from the manifold (Front2VOF), in 2D and 3D
.. (built with cwrap, linked with libb-amr):
..   one box against the exact fraction under a plane: a large element in
..   the plane m . x = alpha (the paper's 3.1), m along an axis, two, three;
..   both orientations; alpha through the corners (1/3, 1/2, 2/3), 0 and 1
..   (the box wholly in or out: unknown, then hm_vof_inside ()); a slab
..   between two parallel elements (3.2)
..   on a tree refined over half its domain, every cell against a fraction
..   computed another way: 2D a convex polygon, a star (several pieces in a
..   cell), the same reversed (1 - c), a water line cut to the domain (the
..   phase along the bottom); 3D a tilted plane cut to the domain (the plane
..   formula, cell by cell), an octahedron (its tetrahedra cut to each cell
..   by hm_clip ()); totals against the exact areas and volumes, a sphere
..   of triangles (its volume by the divergence theorem)
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "tree.h"
#include "block-ghosts.h"
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

/* a fraction stored in c: Real */
static const double tol_c = sizeof (Real) == sizeof (float) ? 2e-7 : 1e-12;

/* -------- the exact fraction under a plane -------- */

/*
.. The fraction of the unit box [0, 1]^D with m . u <= alpha: for k
.. positive components, 1 / (k! prod m) sum over the corners v of the
.. k-cube of (-1)^|v| max (alpha - m . v, 0)^k; a negative component
.. turned (u -> 1 - u), a zero one left out.
*/
static double plane_fraction (const double * m0, double alpha)
{
  double m[3];
  int k = 0;
  for (int d = 0; d < D; d++) {
    double c = m0[d];
    if (c < 0)
      alpha -= c, c = -c;
    if (c > 0)
      m[k++] = c;
  }
  if (!k)
    return alpha >= 0 ? 1 : 0;
  /* outside [0, sum m], 0 or 1; above the middle, by the box's symmetry
     (u -> 1 - u): the alternating sum cancels badly far from 0 */
  double sum = 0;
  for (int i = 0; i < k; i++)
    sum += m[i];
  if (alpha <= 0)
    return 0;
  if (alpha >= sum)
    return 1;
  int upper = alpha > sum / 2;
  if (upper)
    alpha = sum - alpha;
  double prod = 1, fact = 1, s = 0;
  for (int i = 0; i < k; i++)
    prod *= m[i], fact *= i + 1;
  for (int v = 0; v < 1 << k; v++) {
    double mv = 0;
    int bits = 0;
    for (int i = 0; i < k; i++)
      if (v >> i & 1)
        mv += m[i], bits++;
    double t = alpha - mv;
    if (t > 0)
      s += (bits & 1 ? -1 : 1) * pow (t, k);
  }
  double f = s / (fact * prod);
  if (upper)
    f = 1 - f;
  return f < 0 ? 0 : f > 1 ? 1 : f;
}

/* -------- one box -------- */

/* the element in the plane m . x = alpha about c0, large; its normal +m
   (the phase m . x <= alpha) or -m */
static uint32_t plane_element (Hmesh * M, const double * m, const double * c0,
                               int flip)
{
  double n = 0, mm[3] = { 0, 0, 0 };
  for (int d = 0; d < D; d++)
    n += m[d]*m[d];
  for (int d = 0; d < D; d++)
    mm[d] = m[d] / sqrt (n);
  uint32_t v[3];
#if D == 2
  double t[2] = { -mm[1], mm[0] };   /* normal (t_y, -t_x) = +m */
  for (int i = 0; i < 2; i++) {
    double x[2], s = i ? 8 : -8;
    for (int d = 0; d < 2; d++)
      x[d] = c0[d] + s*t[d];
    v[flip ? 1 - i : i] = hm_vertex_add (M, x);
  }
  return hm_add (M, 1, v);
#else
  /* e1 normal to m, e2 = m x e1: counterclockwise about m */
  double e1[3], e2[3], a[3] = { 1, 0, 0 };
  if (fabs (mm[0]) > 0.9)
    a[0] = 0, a[1] = 1;
  e1[0] = a[1]*mm[2] - a[2]*mm[1];
  e1[1] = a[2]*mm[0] - a[0]*mm[2];
  e1[2] = a[0]*mm[1] - a[1]*mm[0];
  double l = sqrt (e1[0]*e1[0] + e1[1]*e1[1] + e1[2]*e1[2]);
  for (int d = 0; d < 3; d++)
    e1[d] /= l;
  e2[0] = mm[1]*e1[2] - mm[2]*e1[1];
  e2[1] = mm[2]*e1[0] - mm[0]*e1[2];
  e2[2] = mm[0]*e1[1] - mm[1]*e1[0];
  for (int i = 0; i < 3; i++) {
    double th = 2*M_PI*i/3, x[3], r = 8/sqrt (3.);
    for (int d = 0; d < 3; d++)
      x[d] = c0[d] + r*(cos (th)*e1[d] + sin (th)*e2[d]);
    v[flip ? 2 - i : i] = hm_vertex_add (M, x);
  }
  return hm_add (M, 2, v);
#endif
}

static void test_box (void)
{
  double ms[3][3] = {
    { 1, 0, 0 },
    { 0.5, 0.5, 0 },
    { D == 2 ? 0.3 : 1./3, D == 2 ? 0.7 : 1./3, D == 2 ? 0 : 1./3 },
  };
  double lo[3] = { 0, 0, 0 }, emax = 0;
  int cases = 0;
  for (int k = 0; k < 3; k++)
    for (int flip = 0; flip < 2; flip++)
      for (int i = 0; i <= 24; i++) {
        double alpha = i <= 20 ? i / 20. : (double []) { 1./3, 0.5, 2./3,
                                                          0.55 } [i - 21];
        Hmesh * M = hmesh_new (D - 1);
        double c0[3] = { alpha, alpha, alpha };
        uint32_t e = plane_element (M, ms[k], c0, flip);
        double want = plane_fraction (ms[k], alpha);
        if (flip)
          want = 1 - want;
        int known;
        double f = hm_vof_box (M, lo, 1., &e, 1, &known);
        if (!known) {
          /* wholly in or out: the element only touches the box */
          double ctr[3] = { 0.5, 0.5, 0.5 };
          f = hm_vof_inside (M, ctr, &e, 1);
          CHECK (fabs (want - round (want)) < 1e-12, "box: m %d alpha %g "
                 "unknown, the fraction %g", k, alpha, want);
        }
        emax = fmax (emax, fabs (f - want));
        CHECK (fabs (f - want) < 1e-12, "box: m %d, %s, alpha %g: %.17g, "
               "not %.17g", k, flip ? "-m" : "+m", alpha, f, want);
        cases++;
        hmesh_free (M);
      }

  /* a slab between two parallel elements: plane beta, normal -m (the phase
     above), plane alpha, normal +m */
  for (int i = 0; i <= 20; i++) {
    double alpha = 0.1 + 0.05*i, beta = alpha - 0.35;
    Hmesh * M = hmesh_new (D - 1);
    double ca[3] = { alpha, alpha, alpha }, cb[3] = { beta, beta, beta };
    uint32_t e[2] = { plane_element (M, ms[2], ca, 0),
                      plane_element (M, ms[2], cb, 1) };
    double want = plane_fraction (ms[2], alpha) - plane_fraction (ms[2], beta);
    int known;
    double f = hm_vof_box (M, lo, 1., e, 2, &known);
    if (!known)
      f = hm_vof_inside (M, (double []) { 0.5, 0.5, 0.5 }, e, 2);
    emax = fmax (emax, fabs (f - want));
    CHECK (fabs (f - want) < 1e-12, "slab: alpha %g: %.17g, not %.17g",
           alpha, f, want);
    cases++;
    hmesh_free (M);
  }
  printf ("one box: %d cases, the largest error %.2g\n", cases, emax);
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

/* the tree: level 1 (32 cells a side), the lower half in x at level 2 */
static void grid (void)
{
  tree (3, 4);
  refine_half (2);
}

#if D == 2
/* the area of the polygon p [0 ... n - 1] (any, counterclockwise) inside
   the box [x0, x0 + h] x [y0, y0 + h]: clipped by its four sides (a convex
   window: right for any polygon), then the shoelace */
static double polygon_in_box (const double (* p)[2], int n, double x0,
                              double y0, double h)
{
  double a[512][2], b[512][2];
  int na = n;
  for (int i = 0; i < n; i++)
    a[i][0] = p[i][0], a[i][1] = p[i][1];
  for (int side = 0; side < 4; side++) {
    int ax = side / 2, up = side % 2, nb = 0;
    double w = (ax ? y0 : x0) + (up ? h : 0);
    for (int i = 0; i < na; i++) {
      double * P = a[i], * Q = a[(i + 1) % na];
      double dp = up ? w - P[ax] : P[ax] - w, dq = up ? w - Q[ax] : Q[ax] - w;
      if (dp >= 0)
        b[nb][0] = P[0], b[nb][1] = P[1], nb++;
      if ((dp > 0 && dq < 0) || (dp < 0 && dq > 0)) {
        double t = dp / (dp - dq);
        b[nb][0] = P[0] + t*(Q[0] - P[0]);
        b[nb][1] = P[1] + t*(Q[1] - P[1]);
        nb++;
      }
    }
    for (int i = 0; i < nb; i++)
      a[i][0] = b[i][0], a[i][1] = b[i][1];
    na = nb;
  }
  double s = 0;
  for (int i = 0; i < na; i++)
    s += a[i][0]*a[(i + 1) % na][1] - a[(i + 1) % na][0]*a[i][1];
  return 0.5*s;
}

/* the manifold of the closed polygon p, counterclockwise (or not) */
static Hmesh * polygon_mesh (const double (* p)[2], int n, int reverse)
{
  Hmesh * M = hmesh_new (1);
  for (int i = 0; i < n; i++)
    hm_vertex_add (M, p[i]);
  for (uint32_t i = 0; i < (uint32_t) n; i++) {
    uint32_t a = i, b = (i + 1) % (uint32_t) n;
    hm_add (M, 1, reverse ? (uint32_t []) { b, a } : (uint32_t []) { a, b });
  }
  return M;
}

/* c against the polygon's area in each cell; returns the largest error */
static double against_polygon (scalar c, const double (* p)[2], int n,
                               int reverse, double * total)
{
  double emax = 0, sum = 0;
  foreach (reduction (max:emax) reduction (+:sum)) {
    double h = Delta;
    double want = polygon_in_box (p, n, x - h/2, y - h/2, h) / (h*h);
    if (reverse)
      want = 1 - want;
    emax = fmax (emax, fabs (c[] - want));
    sum += c[]*h*h;
  }
  *total = sum;
  return emax;
}

static void test_tree2 (void)
{
  grid ();
  scalar c = scalar_new ();
  double p[400][2], total;
  int n;

  /* a convex polygon, 97 points about (0.52, 0.47) */
  n = 97;
  for (int i = 0; i < n; i++) {
    double t = 2*M_PI*i/n;
    p[i][0] = 0.52 + 0.31*cos (t), p[i][1] = 0.47 + 0.31*sin (t);
  }
  double area = polygon_in_box ((const double (*)[2]) p, n, 0, 0, 1);
  for (int reverse = 0; reverse < 2; reverse++) {
    Hmesh * M = polygon_mesh ((const double (*)[2]) p, n, reverse);
    long left = hm_vof (M, c);
    double e = against_polygon (c, (const double (*)[2]) p, n, reverse,
                                &total);
    double want = reverse ? 1 - area : area;
    CHECK (left == 0 && e < tol_c && fabs (total - want) < 1e-6,
           "circle%s: %ld left, error %g, total %.15g, not %.15g",
           reverse ? " reversed" : "", left, e, total, want);
    printf ("2D circle%s: error %.2g a cell, total off by %.2g\n",
            reverse ? " reversed" : "", e, fabs (total - want));
    hmesh_free (M);
  }

  /* a star, seven points: several pieces in a cell, cells between arms */
  n = 280;
  for (int i = 0; i < n; i++) {
    double t = 2*M_PI*i/n, r = 0.28 + 0.12*cos (7*t);
    p[i][0] = 0.5 + r*cos (t), p[i][1] = 0.5 + r*sin (t);
  }
  area = polygon_in_box ((const double (*)[2]) p, n, 0, 0, 1);
  Hmesh * M = polygon_mesh ((const double (*)[2]) p, n, 0);
  long left = hm_vof (M, c);
  double e = against_polygon (c, (const double (*)[2]) p, n, 0, &total);
  CHECK (left == 0 && e < tol_c && fabs (total - area) < 1e-6,
         "star: %ld left, error %g, total %.15g, not %.15g", left, e, total,
         area);
  printf ("2D star: error %.2g a cell, total off by %.2g\n", e,
          fabs (total - area));
  hmesh_free (M);

  /* a water line from x = 1.2 to -0.2, the phase below (the normal up),
     cut to the domain: the region under it in the square */
  M = hmesh_new (1);
  n = 60;
  for (int i = 0; i <= n; i++) {
    double x = 1.2 - 1.4*i/n, y = 0.4 + 0.1*sin (2*M_PI*1.3*x);
    hm_vertex_add (M, (double []) { x, y });
  }
  for (uint32_t i = 0; i < (uint32_t) n; i++)
    hm_add (M, 1, (uint32_t []) { i, i + 1 });
  CHECK (hm_clip (M, (double []) { 0, 0 }, 1.) == 0, "water: clip");
  /* the polygon under it: the line within [0, 1], then the bottom */
  int k = 0;
  for (int i = 0; i <= n; i++) {
    double x = 1.2 - 1.4*i/n;
    if (x >= 0 && x <= 1)
      p[k][0] = x, p[k][1] = 0.4 + 0.1*sin (2*M_PI*1.3*x), k++;
    double xn = 1.2 - 1.4*(i + 1)/n;
    if (i < n && ((x > 1 && xn < 1) || (x > 0 && xn < 0))) {
      double w = x > 1 && xn < 1 ? 1 : 0, t = (w - x) / (xn - x);
      double y0 = 0.4 + 0.1*sin (2*M_PI*1.3*x),
        y1 = 0.4 + 0.1*sin (2*M_PI*1.3*xn);
      p[k][0] = w, p[k][1] = y0 + t*(y1 - y0), k++;
    }
  }
  p[k][0] = 0, p[k][1] = 0, k++;
  p[k][0] = 1, p[k][1] = 0, k++;
  area = polygon_in_box ((const double (*)[2]) p, k, 0, 0, 1);
  left = hm_vof (M, c);
  e = against_polygon (c, (const double (*)[2]) p, k, 0, &total);
  CHECK (left == 0 && e < tol_c && fabs (total - area) < 1e-6,
         "water: %ld left, error %g, total %.15g, not %.15g", left, e,
         total, area);
  printf ("2D water line, clipped: error %.2g a cell, total off by %.2g\n",
          e, fabs (total - area));
  hmesh_free (M);

  /* no manifold: every cell left, at 0 */
  M = hmesh_new (1);
  left = hm_vof (M, c);
  long cells = (long) tree_get_nleaves () * tree_get_block_m () *
    tree_get_block_m ();
  CHECK (left == cells, "empty: %ld left, not every cell, %ld", left, cells);
  hmesh_free (M);
  scalar_free (c);
  tree_destroy ();
}
#endif

#if D == 3
/* the octahedron |x - c|_1 <= r as tetrahedra from its centre, cut to the
   box by hm_clip (): the volume of it in the box */
static double octahedron_in_box (const double * c0, double r,
                                 const double * lo, double h)
{
  Hmesh * T = hmesh_new (3);
  uint32_t o = hm_vertex_add (T, c0), v[6];
  for (int i = 0; i < 6; i++) {
    double x[3] = { c0[0], c0[1], c0[2] };
    x[i / 2] += i % 2 ? r : -r;
    v[i] = hm_vertex_add (T, x);
  }
  for (int sx = 0; sx < 2; sx++)
    for (int sy = 0; sy < 2; sy++)
      for (int sz = 0; sz < 2; sz++) {
        uint32_t t[4] = { o, v[sx], v[2 + sy], v[4 + sz] };
        uint32_t s = hm_add (T, 3, t);
        if (hm_volume_signed (T, s) < 0) {
          hm_remove (T, 3, s);
          uint32_t w = t[0];
          t[0] = t[1], t[1] = w;
          hm_add (T, 3, t);
        }
      }
  hm_clip (T, lo, h);
  double vol = 0;
  hm_foreach (T, 3, s)
    vol += hm_volume_signed (T, s);
  hmesh_free (T);
  return vol;
}

/* its surface: eight triangles, outward */
static Hmesh * octahedron (const double * c0, double r)
{
  Hmesh * M = hmesh_new (2);
  uint32_t v[6];
  for (int i = 0; i < 6; i++) {
    double x[3] = { c0[0], c0[1], c0[2] };
    x[i / 2] += i % 2 ? r : -r;
    v[i] = hm_vertex_add (M, x);
  }
  for (int sx = 0; sx < 2; sx++)
    for (int sy = 0; sy < 2; sy++)
      for (int sz = 0; sz < 2; sz++) {
        uint32_t t[3] = { v[sx], v[2 + sy], v[4 + sz] };
        /* outward: (b - a) x (c - a) along (sx, sy, sz) signs */
        if ((sx + sy + sz) % 2 == 0) {
          uint32_t w = t[0];
          t[0] = t[1], t[1] = w;
        }
        hm_add (M, 2, t);
      }
  return M;
}

/* a sphere of triangles, latitude and longitude, outward */
static Hmesh * sphere (const double * c0, double r, int nt, int np,
                       double * volume)
{
  Hmesh * M = hmesh_new (2);
  uint32_t top = hm_vertex_add (M, (double []) { c0[0], c0[1], c0[2] + r });
  uint32_t bottom = hm_vertex_add (M, (double []) { c0[0], c0[1],
                                                    c0[2] - r });
  uint32_t ring[64][64];
  for (int i = 1; i < nt; i++)
    for (int j = 0; j < np; j++) {
      double th = M_PI*i/nt, ph = 2*M_PI*j/np + 0.1*i;
      ring[i][j] = hm_vertex_add (M, (double []) {
          c0[0] + r*sin (th)*cos (ph), c0[1] + r*sin (th)*sin (ph),
          c0[2] + r*cos (th) });
    }
  for (int j = 0; j < np; j++) {
    int k = (j + 1) % np;
    hm_add (M, 2, (uint32_t []) { top, ring[1][j], ring[1][k] });
    hm_add (M, 2, (uint32_t []) { bottom, ring[nt - 1][k],
                                  ring[nt - 1][j] });
    for (int i = 1; i + 1 < nt; i++) {
      hm_add (M, 2, (uint32_t []) { ring[i][j], ring[i + 1][j],
                                    ring[i + 1][k] });
      hm_add (M, 2, (uint32_t []) { ring[i][j], ring[i + 1][k],
                                    ring[i][k] });
    }
  }
  /* its volume: the divergence theorem, from c0 */
  double v = 0;
  hm_foreach (M, 2, t) {
    double a[3][3];
    for (int i = 0; i < 3; i++)
      for (int d = 0; d < 3; d++)
        a[i][d] = hm_x (M, d)[hm_vertex (M, 2, t, i)] - c0[d];
    v += (a[0][0]*(a[1][1]*a[2][2] - a[1][2]*a[2][1])
          - a[0][1]*(a[1][0]*a[2][2] - a[1][2]*a[2][0])
          + a[0][2]*(a[1][0]*a[2][1] - a[1][1]*a[2][0])) / 6;
  }
  *volume = v;
  return M;
}

static void test_tree3 (void)
{
  grid ();
  scalar c = scalar_new ();

  /* a tilted plane, z = 0.43 + 0.17 x - 0.11 y, the phase below, its
     triangles larger than the cells, cut to the domain: every cell by the
     plane formula */
  Hmesh * M = hmesh_new (2);
  int n = 9;
  for (int j = 0; j <= n; j++)
    for (int i = 0; i <= n; i++) {
      double x = -0.3 + 1.6*i/n, y = -0.3 + 1.6*j/n;
      hm_vertex_add (M, (double []) { x, y, 0.43 + 0.17*x - 0.11*y });
    }
  uint32_t w = (uint32_t) n + 1;
  for (uint32_t j = 0; j < (uint32_t) n; j++)
    for (uint32_t i = 0; i < (uint32_t) n; i++) {
      uint32_t a = j*w + i;
      hm_add (M, 2, (uint32_t []) { a, a + 1, a + w + 1 });
      hm_add (M, 2, (uint32_t []) { a, a + w + 1, a + w });
    }
  CHECK (hm_clip (M, (double []) { 0, 0, 0 }, 1.) == 0, "plane: clip");
  long left = hm_vof (M, c);
  double emax = 0, sum = 0;
  foreach (reduction (max:emax) reduction (+:sum)) {
    double h = Delta, m[3] = { -0.17, 0.11, 1 };
    /* z <= 0.43 + 0.17 x - 0.11 y in the cell's unit coordinates */
    double alpha = (0.43 + 0.17*(x - h/2) - 0.11*(y - h/2) - (z - h/2)) / h;
    double want = plane_fraction (m, alpha);
    emax = fmax (emax, fabs (c[] - want));
    sum += c[]*h*h*h;
  }
  double want = 0.43 + 0.17*0.5 - 0.11*0.5;
  CHECK (left == 0 && emax < tol_c && fabs (sum - want) < 1e-6,
         "plane: %ld left, error %g, total %.15g, not %.15g", left, emax,
         sum, want);
  printf ("3D plane, clipped: error %.2g a cell, total off by %.2g\n", emax,
          fabs (sum - want));
  hmesh_free (M);

  /* an octahedron: every cell against its tetrahedra cut to the cell */
  double c0[3] = { 0.47, 0.52, 0.49 }, r = 0.33;
  M = octahedron (c0, r);
  left = hm_vof (M, c);
  emax = 0, sum = 0;
  foreach (reduction (max:emax) reduction (+:sum)) {
    double h = Delta, lo[3] = { x - h/2, y - h/2, z - h/2 };
    double v = octahedron_in_box (c0, r, lo, h) / (h*h*h);
    emax = fmax (emax, fabs (c[] - v));
    sum += c[]*h*h*h;
  }
  want = 4*r*r*r/3;
  CHECK (left == 0 && emax < tol_c && fabs (sum - want) < 1e-6,
         "octahedron: %ld left, error %g, total %.15g, not %.15g", left,
         emax, sum, want);
  printf ("3D octahedron: error %.2g a cell, total off by %.2g\n", emax,
          fabs (sum - want));
  hmesh_free (M);

  /* a sphere of triangles: the total */
  double vol;
  M = sphere ((double []) { 0.51, 0.48, 0.5 }, 0.3, 24, 48, &vol);
  left = hm_vof (M, c);
  sum = 0;
  double cmin = 1, cmax = 0;
  foreach (reduction (max:cmax) reduction (min:cmin) reduction (+:sum)) {
    sum += c[]*Delta*Delta*Delta;
    cmin = fmin (cmin, c[]), cmax = fmax (cmax, c[]);
  }
  CHECK (left == 0 && fabs (sum - vol) < 1e-6 && cmin >= 0 && cmax <= 1,
         "sphere: %ld left, total %.15g, not %.15g", left, sum, vol);
  printf ("3D sphere (%u triangles): total off by %.2g\n", hm_count (M, 2),
          fabs (sum - vol));
  hmesh_free (M);
  scalar_free (c);
  tree_destroy ();
}
#endif

int main (void)
{
  test_box ();
#if D == 2
  test_tree2 ();
#else
  test_tree3 ();
#endif
  if (fails)
    printf ("%d failed\n", fails);
  else
    printf ("passed, %dD\n", D);
  return fails != 0;
}
