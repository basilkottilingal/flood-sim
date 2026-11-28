/*
.. The quality of the mesh: see hmesh-quality.h. Points are padded to three
.. coordinates (zeros past D), so the cross product gives areas, volumes and
.. the sines of angles in every dimension, without the cancellation of a
.. Gram determinant for small angles.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "hmesh-quality.h"

#define D HM_DIMENSION
#define DEG (180. / 3.14159265358979323846)

static void coords (const Hmesh * M, uint32_t v, double * p)
{
  p[0] = p[1] = p[2] = 0.;
  for (int d = 0; d < D; d++)
    p[d] = M->x[d][v];
}

static void sub (const double * a, const double * b, double * c)
{
  for (int d = 0; d < 3; d++)
    c[d] = a[d] - b[d];
}

static double dot (const double * a, const double * b)
{
  return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

static void cross (const double * a, const double * b, double * c)
{
  c[0] = a[1]*b[2] - a[2]*b[1];
  c[1] = a[2]*b[0] - a[0]*b[2];
  c[2] = a[0]*b[1] - a[1]*b[0];
}

static double norm (const double * a)
{
  return sqrt (dot (a, a));
}

/* the angle between u and w, in degrees, by atan2 (|u x w|, u . w) */
static double angle (const double * u, const double * w)
{
  double c[3];
  cross (u, w, c);
  return atan2 (norm (c), dot (u, w)) * DEG;
}

/* the angle at p between q and r */
static double angle_at (const double * p, const double * q, const double * r)
{
  double u[3], w[3];
  sub (q, p, u);
  sub (r, p, w);
  return angle (u, w);
}

static double tri_area (const double * a, const double * b, const double * c)
{
  double u[3], w[3], n[3];
  sub (b, a, u);
  sub (c, a, w);
  cross (u, w, n);
  return 0.5 * norm (n);
}

/* the dihedral angle along the edge ij of a tetrahedron, its other points
   k, l: the angle between e x (pk - pi) and e x (pl - pi) */
static double dihedral (const double * pi, const double * pj,
                        const double * pk, const double * pl)
{
  double e[3], u[3], w[3], nu[3], nw[3];
  sub (pj, pi, e);
  sub (pk, pi, u);
  sub (pl, pi, w);
  cross (e, u, nu);
  cross (e, w, nw);
  return angle (nu, nw);
}

int hm_shape (const Hmesh * M, int k, uint32_t s, HmShape * q)
{
  if (!M || !q || k < 1 || k > M->dim || k > 3 || !hm_alive (M, k, s))
    return -1;
  double p[4][3], l[4][4];
  for (int i = 0; i <= k; i++)
    coords (M, hm_vertex (M, k, s, i), p[i]);
  q->lmin = INFINITY, q->lmax = 0.;
  for (int i = 0; i <= k; i++)
    for (int j = i + 1; j <= k; j++) {
      double e[3];
      sub (p[j], p[i], e);
      l[i][j] = l[j][i] = norm (e);
      q->lmin = fmin (q->lmin, l[i][j]);
      q->lmax = fmax (q->lmax, l[i][j]);
    }
  q->edge_ratio = q->lmin > 0. ? q->lmax / q->lmin : INFINITY;
  if (k == 1) {
    q->measure = q->lmax;
    q->edge_ratio = q->aspect = q->radius_ratio = 1.;
    q->angle_min = q->angle_max = NAN;
    return 0;
  }
  if (k == 2) {
    double A = tri_area (p[0], p[1], p[2]);
    double a = l[1][2], b = l[0][2], c = l[0][1], P = a + b + c;
    q->measure = A;
    double t0 = angle_at (p[0], p[1], p[2]), t1 = angle_at (p[1], p[2], p[0]);
    double t2 = angle_at (p[2], p[0], p[1]);
    q->angle_min = fmin (t0, fmin (t1, t2));
    q->angle_max = fmax (t0, fmax (t1, t2));
    if (A > HM_QUALITY_EPS * q->lmax * q->lmax) {
      q->aspect = q->lmax * P / (4. * sqrt (3.) * A);
      q->radius_ratio = a * b * c * P / (16. * A * A);
    }
    else
      q->aspect = q->radius_ratio = INFINITY;
    return 0;
  }
  /* k == 3 */
  double u[3], v[3], w[3], n[3];
  sub (p[1], p[0], u);
  sub (p[2], p[0], v);
  sub (p[3], p[0], w);
  cross (v, w, n);
  double V = fabs (dot (u, n)) / 6.;
  double S = tri_area (p[1], p[2], p[3]) + tri_area (p[0], p[2], p[3]) +
    tri_area (p[0], p[1], p[3]) + tri_area (p[0], p[1], p[2]);
  q->measure = V;
  static const int pair[6][4] = { {0, 1, 2, 3}, {0, 2, 1, 3}, {0, 3, 1, 2},
                                  {1, 2, 0, 3}, {1, 3, 0, 2}, {2, 3, 0, 1} };
  q->angle_min = INFINITY, q->angle_max = -INFINITY;
  for (int e = 0; e < 6; e++) {
    const int * c = pair[e];
    double t = dihedral (p[c[0]], p[c[1]], p[c[2]], p[c[3]]);
    q->angle_min = fmin (q->angle_min, t);
    q->angle_max = fmax (q->angle_max, t);
  }
  double lmax = q->lmax;
  if (V > HM_QUALITY_EPS * lmax * lmax * lmax) {
    double pp = l[0][1] * l[2][3], qq = l[0][2] * l[1][3];
    double rr = l[0][3] * l[1][2];
    double prod = (pp + qq + rr) * (pp + qq - rr) * (pp - qq + rr) *
      (-pp + qq + rr);
    double R = sqrt (fmax (prod, 0.)) / (24. * V);
    q->aspect = lmax * S / (6. * sqrt (6.) * V);
    q->radius_ratio = R * S / (9. * V);
  }
  else
    q->aspect = q->radius_ratio = INFINITY;
  return 0;
}

/* -------- statistics -------- */

typedef struct {
  HmStat * st;
  double   sum;
} Acc;

static void acc_add (Acc * a, double x, uint32_t id)
{
  HmStat * st = a->st;
  if (st->n == 0 || x < st->min)
    st->min = x, st->argmin = id;
  if (st->n == 0 || x > st->max)
    st->max = x, st->argmax = id;
  st->n++;
  a->sum += x;
}

static void acc_end (Acc * a)
{
  HmStat * st = a->st;
  if (st->n)
    st->mean = a->sum / st->n;
  else {
    st->min = st->mean = st->max = 0.;
    st->argmin = st->argmax = HM_NONE;
  }
}

enum { MEASURE, EDGE, EDGE_RATIO, ASPECT, RADIUS, AMIN, AMAX, TURNING,
       VALENCE, EDGE_CELL, NACC };

int hm_quality (const Hmesh * M, HmQuality * Q)
{
  if (!M || !Q)
    return -1;
  memset (Q, 0, sizeof (*Q));
  int m = M->dim;
  Q->dim = m;
  Q->points = hm_count (M, 0);
  Q->simplices = hm_count (M, m);
  Acc a[NACC];
  HmStat * st[NACC] = { &Q->measure, &Q->edge, &Q->edge_ratio, &Q->aspect,
                        &Q->radius_ratio, &Q->angle_min, &Q->angle_max,
                        &Q->turning, &Q->valence, &Q->edge_cell };
  for (int i = 0; i < NACC; i++)
    a[i].st = st[i], a[i].sum = 0.;

  if (m >= 1) {
    uint32_t end = hm_end (M, 0);
    uint32_t * val = calloc (end ? end : 1, sizeof (uint32_t));
    uint32_t * nb = m == 1 ? malloc (2 * (size_t) (end ? end : 1) *
                                     sizeof (uint32_t)) : NULL;
    if (!val || (m == 1 && !nb)) {
      free (val), free (nb);
      return -1;
    }
    /* the edges: lengths, valences, and a curve's neighbours */
    hm_foreach (M, 1, e) {
      uint32_t v[2] = { hm_vertex (M, 1, e, 0), hm_vertex (M, 1, e, 1) };
      double pa[3], pb[3], d[3];
      coords (M, v[0], pa);
      coords (M, v[1], pb);
      sub (pb, pa, d);
      double len = norm (d);
      acc_add (&a[EDGE], len, e);
#ifdef EULERIAN_LAGRANGIAN
      if (M->cell_size) {
        double h = fmin (M->cell_size (pa, M->cell_ctx),
                         M->cell_size (pb, M->cell_ctx));
        acc_add (&a[EDGE_CELL], len / h, e);
      }
#endif
      for (int i = 0; i < 2; i++) {
        if (nb && val[v[i]] < 2)
          nb[2 * (size_t) v[i] + val[v[i]]] = v[1 - i];
        val[v[i]]++;
      }
    }
    hm_foreach (M, 0, v) {
      acc_add (&a[VALENCE], (double) val[v], v);
      if (nb && val[v] == 2) {
        double p[3], q0[3], q1[3], u[3], w[3];
        coords (M, v, p);
        coords (M, nb[2 * (size_t) v], q0);
        coords (M, nb[2 * (size_t) v + 1], q1);
        sub (q0, p, u);
        sub (q1, p, w);
        double lu = norm (u), lw = norm (w);
        if (lu > 0. && lw > 0.) {
          acc_add (&a[EDGE_RATIO], fmax (lu, lw) / fmin (lu, lw), v);
          acc_add (&a[TURNING], 180. - angle (u, w), v);
        }
      }
    }
    free (val), free (nb);

    /* the top simplices */
    hm_foreach (M, m, s) {
      HmShape q;
      hm_shape (M, m, s, &q);
      acc_add (&a[MEASURE], q.measure, s);
      if (m < 2)
        continue;
      if (isinf (q.aspect)) {
        Q->degenerate++;
        continue;
      }
      acc_add (&a[EDGE_RATIO], q.edge_ratio, s);
      acc_add (&a[ASPECT], q.aspect, s);
      acc_add (&a[RADIUS], q.radius_ratio, s);
      acc_add (&a[AMIN], q.angle_min, s);
      acc_add (&a[AMAX], q.angle_max, s);
    }
  }
  for (int i = 0; i < NACC; i++)
    acc_end (&a[i]);
  return 0;
}

static void row (FILE * fp, const char * name, const HmStat * st)
{
  if (st->n)
    fprintf (fp, "  %-13s %10u %13.6g %13.6g %13.6g\n", name, st->n,
             st->min, st->mean, st->max);
}

void hm_quality_print (FILE * fp, const HmQuality * Q)
{
  static const char * measure[4] = { "-", "length", "area", "volume" };
  fprintf (fp, "hmesh quality: m = %d, %u points, %u simplices, "
           "%u degenerate\n", Q->dim, Q->points, Q->simplices,
           Q->degenerate);
  fprintf (fp, "  %-13s %10s %13s %13s %13s\n", "", "n", "min", "mean",
           "max");
  row (fp, measure[Q->dim < 0 || Q->dim > 3 ? 0 : Q->dim], &Q->measure);
  row (fp, "edge", &Q->edge);
  row (fp, "edge ratio", &Q->edge_ratio);
  row (fp, "aspect ratio", &Q->aspect);
  row (fp, "radius ratio", &Q->radius_ratio);
  row (fp, Q->dim == 3 ? "min dihedral" : "min angle", &Q->angle_min);
  row (fp, Q->dim == 3 ? "max dihedral" : "max angle", &Q->angle_max);
  row (fp, "turning", &Q->turning);
  row (fp, "valence", &Q->valence);
  row (fp, "edge / cell", &Q->edge_cell);
}
