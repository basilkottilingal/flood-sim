/*
.. The manifold from a function, the interface f = 0: see hmesh-implicit.h.
..
.. Marching simplices: each lattice cell split into D! simplices along its
.. main diagonal (Freudenthal-Kuhn), the interface's piece in each from the
.. signs at its corners. As every corner path is monotone, any two corners
.. of a simplex are a lattice node a and a + delta, delta in {0, 1}^D: the
.. lattice edge (a, delta) names the point on it, found once.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "hmesh-implicit.h"

#define D HM_DIMENSION

static double eval (HmFunction f, const double * x)
{
#if D == 1
  return f (x[0]);
#elif D == 2
  return f (x[0], x[1]);
#else
  return f (x[0], x[1], x[2]);
#endif
}

/* -------- the points already made, by lattice edge or node -------- */

typedef struct {
  uint64_t * key;             /* key + 1, 0 empty */
  uint32_t * id;
  size_t     cap, n;
} Map;

static uint32_t * map_find (Map * m, uint64_t key, int insert)
{
  if (insert && 2 * (m->n + 1) > m->cap) {
    size_t cap = m->cap ? 2 * m->cap : 1024;
    uint64_t * k = calloc (cap, sizeof (uint64_t));
    uint32_t * v = malloc (cap * sizeof (uint32_t));
    if (!k || !v) {
      fprintf (stderr, "hmesh: out of memory\n");
      abort ();
    }
    for (size_t i = 0; i < m->cap; i++)
      if (m->key[i]) {
        size_t j = (m->key[i] * 0x9E3779B97F4A7C15ull) & (cap - 1);
        while (k[j])
          j = (j + 1) & (cap - 1);
        k[j] = m->key[i], v[j] = m->id[i];
      }
    free (m->key), free (m->id);
    m->key = k, m->id = v, m->cap = cap;
  }
  if (!m->cap)
    return NULL;
  uint64_t k1 = key + 1;
  size_t j = (k1 * 0x9E3779B97F4A7C15ull) & (m->cap - 1);
  while (m->key[j] && m->key[j] != k1)
    j = (j + 1) & (m->cap - 1);
  if (m->key[j])
    return &m->id[j];
  if (!insert)
    return NULL;
  m->key[j] = k1;
  m->id[j] = 0;                        /* not made yet: ids are stored + 1 */
  m->n++;
  return &m->id[j];
}

typedef struct {
  Hmesh    * M;
  HmFunction f;
  double     lo[3], h;
  int        n;
  Map        map;
  int        fail;
} March;

static uint64_t node_index (const March * C, const int * a)
{
  uint64_t i = 0;
  for (int d = D - 1; d >= 0; d--)
    i = i * (uint64_t) (C->n + 1) + (uint64_t) a[d];
  return i;
}

static void node_position (const March * C, const int * a, double * x)
{
  for (int d = 0; d < D; d++)
    x[d] = C->lo[d] + a[d] * C->h;
}

/*
.. The root of f on the segment from xa to xb, fa and fb of opposite signs
.. (neither 0): regula falsi with the Illinois halving, to round-off
*/
static double root (HmFunction f, const double * xa, const double * xb,
                    double fa, double fb)
{
  double t0 = 0, t1 = 1, g0 = fa, g1 = fb, t = 0.5, x[3];
  int side = 0;
  for (int it = 0; it < 100; it++) {
    t = (t0 * g1 - t1 * g0) / (g1 - g0);
    if (!(t > t0 && t < t1))
      t = 0.5 * (t0 + t1);
    for (int d = 0; d < D; d++)
      x[d] = xa[d] + t * (xb[d] - xa[d]);
    double g = eval (f, x);
    if (g == 0)
      break;
    if ((g > 0) == (g1 > 0)) {
      t1 = t, g1 = g;
      if (side == -1)
        g0 *= 0.5;
      side = -1;
    }
    else {
      t0 = t, g0 = g;
      if (side == 1)
        g1 *= 0.5;
      side = 1;
    }
    if (t1 - t0 <= 1e-15)
      break;
  }
  return t;
}

/* the point on the lattice edge from node a to node b = a + delta, whose
   ends differ in sign: made once */
static inline uint32_t point_on (March * C, const int * a, double fa,
                                 const int * b, double fb)
{
  uint64_t key;
  int at = -1;                         /* a root at a node: 0 a, 1 b */
  if (fa == 0)
    at = 0, key = node_index (C, a) << 3;
  else if (fb == 0)
    at = 1, key = node_index (C, b) << 3;
  else {
    uint64_t delta = 0;
    for (int d = 0; d < D; d++)
      delta |= (uint64_t) (b[d] - a[d]) << d;
    key = node_index (C, a) << 3 | delta;
  }
  uint32_t * slot = map_find (&C->map, key, 1);
  if (*slot)                           /* made before (its id + 1) */
    return *slot - 1;
  double xa[3], xb[3], x[3];
  node_position (C, a, xa);
  node_position (C, b, xb);
  if (at == 0)
    memcpy (x, xa, sizeof (x));
  else if (at == 1)
    memcpy (x, xb, sizeof (x));
  else {
    double t = root (C->f, xa, xb, fa, fb);
    for (int d = 0; d < D; d++)
      x[d] = xa[d] + t * (xb[d] - xa[d]);
  }
  uint32_t v = hm_vertex_add (C->M, x);
  if (v == HM_NONE) {
    C->fail = 1;
    v = 0;
  }
  *slot = v + 1;
  return v;
}

/* the simplex of points p [0 ... m], its normal towards out (the outside
   corners' centre minus the inside corners'); none if it is flat */
static inline long add_oriented (March * C, uint32_t * p, const double * out)
{
  Hmesh * M = C->M;
  int m = D - 1;
  for (int i = 0; i <= m; i++)
    for (int j = 0; j < i; j++)
      if (p[i] == p[j])
        return 0;
  double x[3][3] = { { 0 } }, nrm[3] = { 0, 0, 0 };
  for (int i = 0; i <= m; i++)
    for (int d = 0; d < D; d++)
      x[i][d] = M->x[d][p[i]];
#if D == 2
  nrm[0] = x[1][1] - x[0][1], nrm[1] = - (x[1][0] - x[0][0]);
#else
  double e1[3], e2[3];
  for (int d = 0; d < 3; d++)
    e1[d] = x[1][d] - x[0][d], e2[d] = x[2][d] - x[0][d];
  nrm[0] = e1[1]*e2[2] - e1[2]*e2[1];
  nrm[1] = e1[2]*e2[0] - e1[0]*e2[2];
  nrm[2] = e1[0]*e2[1] - e1[1]*e2[0];
#endif
  double s = 0;
  for (int d = 0; d < D; d++)
    s += nrm[d] * out[d];
  if (s == 0)
    return 0;
  if (s < 0) {
    uint32_t t = p[0];
    p[0] = p[1], p[1] = t;
  }
  if (hm_add (M, m, p) == HM_NONE) {
    C->fail = 1;
    return 0;
  }
  return 1;
}

/* the interface's piece in the lattice simplex of corners c [0 ... D] */
static long march_simplex (March * C, int (* c)[3], const double * fc)
{
#if D == 1
  (void) C, (void) c, (void) fc;
  return 0;
#else
  int in[4], out[4], ni = 0, no = 0;
  for (int i = 0; i <= D; i++) {
    if (fc[i] <= 0)
      in[ni++] = i;
    else
      out[no++] = i;
  }
  if (!ni || !no)
    return 0;
  /* the outside corners' centre minus the inside corners' */
  double dir[3] = { 0, 0, 0 };
  for (int d = 0; d < D; d++) {
    for (int i = 0; i < no; i++)
      dir[d] += c[out[i]][d] / (double) no;
    for (int i = 0; i < ni; i++)
      dir[d] -= c[in[i]][d] / (double) ni;
  }
# define P(i, o) point_on (C, c[(i) < (o) ? (i) : (o)],                 \
                           fc[(i) < (o) ? (i) : (o)],                   \
                           c[(i) < (o) ? (o) : (i)],                    \
                           fc[(i) < (o) ? (o) : (i)])
  /* corners in path order: the lower of two is the earlier */
# if D == 2
  uint32_t p[2];
  int lone = ni == 1 ? in[0] : out[0];
  int k = 0;
  for (int i = 0; i <= 2; i++)
    if (i != lone)
      p[k++] = P (lone, i);
  return add_oriented (C, p, dir);
# else
  if (ni == 1 || no == 1) {
    int lone = ni == 1 ? in[0] : out[0];
    uint32_t p[3];
    int k = 0;
    for (int i = 0; i <= 3; i++)
      if (i != lone)
        p[k++] = P (lone, i);
    return add_oriented (C, p, dir);
  }
  /* two in, two out: a quadrangle, in order about it, in two triangles by
     its shorter diagonal */
  uint32_t q[4] = { P (in[0], out[0]), P (in[0], out[1]),
                    P (in[1], out[1]), P (in[1], out[0]) };
  double d02 = 0, d13 = 0;
  for (int d = 0; d < 3; d++) {
    double a = C->M->x[d][q[0]] - C->M->x[d][q[2]];
    double b = C->M->x[d][q[1]] - C->M->x[d][q[3]];
    d02 += a*a, d13 += b*b;
  }
  uint32_t t1[3], t2[3];
  if (d02 <= d13) {
    t1[0] = q[0], t1[1] = q[1], t1[2] = q[2];
    t2[0] = q[0], t2[1] = q[2], t2[2] = q[3];
  }
  else {
    t1[0] = q[0], t1[1] = q[1], t1[2] = q[3];
    t2[0] = q[1], t2[1] = q[2], t2[2] = q[3];
  }
  return add_oriented (C, t1, dir) + add_oriented (C, t2, dir);
# endif
# undef P
#endif
}

long hm_from_function (Hmesh * M, HmFunction f, const double * lo,
                       double size, int n)
{
  if (D < 2 || M->dim != D - 1 || n < 1 || !(size > 0))
    return -1;
  March C = { .M = M, .f = f, .h = size / n, .n = n };
  uint32_t before = hm_count (M, D - 1);
  for (int d = 0; d < D; d++)
    C.lo[d] = lo[d];

  /* the corner paths of the cell's simplices: one a permutation of the
     axes, each step along one */
  int perm[6][3], np = 0;
#if D == 2
  int base[2][3] = { { 0, 1, 2 }, { 1, 0, 2 } };
  for (np = 0; np < 2; np++)
    memcpy (perm[np], base[np], sizeof (perm[np]));
#else
  int base[6][3] = { { 0, 1, 2 }, { 0, 2, 1 }, { 1, 0, 2 }, { 1, 2, 0 },
                     { 2, 0, 1 }, { 2, 1, 0 } };
  for (np = 0; np < 6; np++)
    memcpy (perm[np], base[np], sizeof (perm[np]));
#endif

  /* f on two layers of nodes at a time (rows in 2D, planes in 3D) */
  size_t side = (size_t) n + 1;
  size_t layer = D == 2 ? side : side * side;
  double * f0 = malloc (layer * sizeof (double)),
    * f1 = malloc (layer * sizeof (double));
  int top = D == 2 ? 1 : 2;            /* the axis of the layers */
  for (int s = 0; s <= n; s++) {
    double * fl = s == 0 ? f0 : f1;
    for (size_t i = 0; i < layer; i++) {
      int a[3] = { 0, 0, 0 };
      a[0] = (int) (i % (size_t) (n + 1));
      if (D == 3)
        a[1] = (int) (i / (size_t) (n + 1));
      a[top] = s;
      double x[3];
      node_position (&C, a, x);
      fl[i] = eval (f, x);
    }
    if (s == 0)
      continue;
    /* the cells between layers s - 1 (f0) and s (f1) */
    size_t cells = D == 2 ? (size_t) n : (size_t) n * (size_t) n;
    for (size_t ci = 0; ci < cells && !C.fail; ci++) {
      int a0[3] = { 0, 0, 0 };
      a0[0] = (int) (ci % (size_t) n);
      if (D == 3)
        a0[1] = (int) (ci / (size_t) n);
      a0[top] = s - 1;
      for (int p = 0; p < np; p++) {
        int c[4][3];
        double fc[4];
        memcpy (c[0], a0, sizeof (c[0]));
        for (int k = 1; k <= D; k++) {
          memcpy (c[k], c[k - 1], sizeof (c[k]));
          c[k][perm[p][k - 1]]++;
        }
        for (int k = 0; k <= D; k++) {
          size_t i = (size_t) c[k][0] +
            (D == 3 ? (size_t) c[k][1] * (size_t) (n + 1) : 0);
          fc[k] = c[k][top] == s ? f1[i] : f0[i];
        }
        march_simplex (&C, c, fc);
      }
    }
    if (C.fail)
      break;
    double * t = f0;
    f0 = f1, f1 = t;
  }
  free (f0), free (f1);
  free (C.map.key), free (C.map.id);
  if (C.fail)
    return -1;
  return (long) hm_count (M, D - 1) - (long) before;
}
