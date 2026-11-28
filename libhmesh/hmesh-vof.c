/*
.. Volume fractions from the manifold, Front2VOF (Pan et al. 2026): see
.. hmesh-vof.h. Built with cwrap (hm_vof () writes a scalar of libb-amr in a
.. foreach ()).
..
.. The paper's algorithms, and where they are here:
..   Algorithm 1 (FLAG)   snap (): a point within HM_VOF_EPS of a face of the
..                        box moved onto it, its faces as bits
..   Algorithm 2 (CLIP3D) clip_polygon (), clip_segment (): an element cut
..                        to the box, axis by axis, the lower then the upper
..                        face; a new point exactly on the plane, with the
..                        faces both ends share
..   Algorithm 3 (F2V)    box_measure (), face_measure (), line_measure ():
..                        the divergence theorem in the box, then on its face
..                        x = 1, then on that face's edge, written once for
..                        2D and 3D; F_v as the wrap where the face (edge)
..                        below gave nothing
..   Algorithm 4 (F2VEF)  hm_vof_inside (): the side of the closest element,
..                        by the closest point on it (not its centroid)
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "tree.h"
#include "block-ghosts.h"
#include "hmesh-vof.h"
#include "hmesh-morton.h"
#include "hmesh-tree.h"

#define D HM_DIMENSION

static void * grow (void * p, size_t n, size_t e)
{
  void * q = realloc (p, (n ? n : 1) * e);
  if (!q) {
    fprintf (stderr, "hmesh: out of memory\n");
    abort ();
  }
  return q;
}

/* -------- the pieces of the elements in a box -------- */

/* a triangle cut by six planes has nine points at most */
#define PV 12

typedef struct {
  int    n;                   /* points: 2 (a segment, 2D), 3 and more
                                 (a polygon, 3D) */
  double u[PV][3];            /* in the box's unit coordinates */
  int    f[PV];               /* the faces of the box each is on: bit 2 a
                                 for u_a = 0, 2 a + 1 for u_a = 1 */
  double A[3], c[3];          /* area vector, out of the phase; centroid */
  int    on_face;             /* every point on one face of the box */
} Piece;

static int bit (int a, int s)
{
  return 1 << (2 * a + s);
}

/* Algorithm 1: onto the face within HM_VOF_EPS, and flagged */
static void snap (double * u, int * f)
{
  for (int a = 0; a < D; a++)
    if (fabs (u[a]) < HM_VOF_EPS)
      u[a] = 0., *f |= bit (a, 0);
    else if (fabs (u[a] - 1.) < HM_VOF_EPS)
      u[a] = 1., *f |= bit (a, 1);
}

/* how far inside the face (a, s) the point is: > 0 in, 0 on, < 0 out */
static double inside (const double * u, int a, int s)
{
  return s ? 1. - u[a] : u[a];
}

/* the point where the edge from p to q crosses the face (a, s), exactly on
   it, with the faces both ends are on */
static void crossing (const Piece * P, int p, int q, int a, int s,
                      double * u, int * f)
{
  double dp = inside (P->u[p], a, s), dq = inside (P->u[q], a, s);
  double t = dp / (dp - dq);
  for (int d = 0; d < 3; d++)
    u[d] = P->u[p][d] + t * (P->u[q][d] - P->u[p][d]);
  u[a] = s;
  *f = (P->f[p] & P->f[q]) | bit (a, s);
  snap (u, f);
}

/* Algorithm 2, one face: Sutherland-Hodgman */
static void clip_polygon (Piece * P, int a, int s)
{
  double u[PV][3];
  int f[PV], n = 0;
  for (int i = 0; i < P->n; i++) {
    int j = (i + 1) % P->n;
    double di = inside (P->u[i], a, s), dj = inside (P->u[j], a, s);
    if (di >= 0 && n < PV) {
      memcpy (u[n], P->u[i], sizeof (u[n]));
      f[n++] = P->f[i];
    }
    if (((di > 0 && dj < 0) || (di < 0 && dj > 0)) && n < PV) {
      crossing (P, i, j, a, s, u[n], &f[n]);
      n++;
    }
  }
  memcpy (P->u, u, sizeof (u));
  memcpy (P->f, f, sizeof (f));
  P->n = n;
}

/* the same for a segment: 0 if nothing of it is strictly inside */
static int clip_segment (Piece * P, int a, int s)
{
  double d0 = inside (P->u[0], a, s), d1 = inside (P->u[1], a, s);
  if (d0 >= 0 && d1 >= 0)
    return 1;
  if (d0 <= 0 && d1 <= 0)
    return 0;
  int o = d0 < 0 ? 0 : 1;
  double u[3];
  int f;
  crossing (P, 1 - o, o, a, s, u, &f);
  memcpy (P->u[o], u, sizeof (u));
  P->f[o] = f;
  return 1;
}

/* the element cut to the box: 0 if nothing of it is in it */
static int element_piece (const Hmesh * M, uint32_t e, const double * lo,
                          double h, Piece * P)
{
  int m = M->dim;
  P->n = m + 1;
  for (int i = 0; i <= m; i++) {
    uint32_t v = hm_vertex (M, m, e, i);
    P->f[i] = 0;
    for (int d = 0; d < 3; d++)
      P->u[i][d] = d < D ? (M->x[d][v] - lo[d]) / h : 0.;
    snap (P->u[i], &P->f[i]);
  }
  for (int a = 0; a < D; a++)
    for (int s = 0; s < 2; s++)
      if (m == 1) {
        if (!clip_segment (P, a, s))
          return 0;
      }
      else {
        clip_polygon (P, a, s);
        if (P->n < 3)
          return 0;
      }

  /* the area vector and the centroid */
  for (int d = 0; d < 3; d++)
    P->A[d] = P->c[d] = 0.;
  if (m == 1) {
    P->A[0] = P->u[1][1] - P->u[0][1];
    P->A[1] = - (P->u[1][0] - P->u[0][0]);
    for (int d = 0; d < 2; d++)
      P->c[d] = 0.5 * (P->u[0][d] + P->u[1][d]);
  }
  else {
    double w = 0;
    for (int i = 1; i + 1 < P->n; i++) {
      double e1[3], e2[3], a[3];
      for (int d = 0; d < 3; d++) {
        e1[d] = P->u[i][d] - P->u[0][d];
        e2[d] = P->u[i + 1][d] - P->u[0][d];
      }
      a[0] = 0.5 * (e1[1]*e2[2] - e1[2]*e2[1]);
      a[1] = 0.5 * (e1[2]*e2[0] - e1[0]*e2[2]);
      a[2] = 0.5 * (e1[0]*e2[1] - e1[1]*e2[0]);
      double s = sqrt (a[0]*a[0] + a[1]*a[1] + a[2]*a[2]);
      for (int d = 0; d < 3; d++) {
        P->A[d] += a[d];
        P->c[d] += s * (P->u[0][d] + P->u[i][d] + P->u[i + 1][d]) / 3.;
      }
      w += s;
    }
    if (w > 0)
      for (int d = 0; d < 3; d++)
        P->c[d] /= w;
  }
  if (P->A[0] == 0 && P->A[1] == 0 && P->A[2] == 0)
    return 0;
  int common = ~0;
  for (int i = 0; i < P->n; i++)
    common &= P->f[i];
  P->on_face = common != 0;
  return 1;
}

typedef struct {
  Piece  * p;
  int      n, cap;
} Pieces;

static void pieces_of (const Hmesh * M, const double * lo, double h,
                       const uint32_t * elem, uint32_t n, Pieces * S)
{
  S->n = 0;
  for (uint32_t i = 0; i < n; i++) {
    if (S->n == S->cap) {
      S->cap = S->cap ? 2 * S->cap : 16;
      S->p = grow (S->p, (size_t) S->cap, sizeof (Piece));
    }
    if (element_piece (M, elem[i], lo, h, &S->p[S->n]))
      S->n++;
  }
}

/* -------- Algorithm 3: the measures -------- */

/* a point where a trace crosses an edge of the box: its coordinate along
   the edge, the sign of its normal along it (out of the phase) */
typedef struct {
  double x;
  int    s;
} Crossing;

static int by_x (const void * a, const void * b)
{
  double x = ((const Crossing *) a)->x, y = ((const Crossing *) b)->x;
  return x < y ? -1 : x > y;
}

/*
.. The length of the phase on [0, 1], from the crossings: sum of s x, plus 1
.. if the end 1 is in the phase, which the last crossing that does not
.. cancel tells (s < 0: the phase lies after it). Unknown if every crossing
.. cancels (the curve only touches): wholly in or out.
*/
static inline double line_measure (Crossing * p, int n, int * known)
{
  qsort (p, (size_t) n, sizeof (Crossing), by_x);
  double v = 0;
  for (int i = 0; i < n; i++)
    v += p[i].s * p[i].x;
  int end = -1;
  for (int i = n - 1; i >= 0 && end < 0; ) {
    int j = i, net = 0;
    while (j >= 0 && p[j].x == p[i].x)
      net += p[j--].s;
    if (net)
      end = net < 0;
    i = j;
  }
  if (end < 0) {
    *known = 0;
    return 0.;
  }
  *known = 1;
  v += end;
  return v < 0 ? 0 : v > 1 ? 1 : v;
}

static inline int sign (double x)
{
  return (x > 0) - (x < 0);
}

/* below this a sum is 0 for the wrap: the round-off of the unit box */
#define WRAP_TOL 1e-13

#define CROSSINGS_MAX 4096

/*
.. The fraction of the face u_a = s of the box in the phase.
..
.. 2D: the face is a segment along b; its crossings, the pieces' ends on
.. it (not at its corners), each of the sign of the piece's normal along b.
..
.. 3D: along b, the divergence theorem in the face: the traces (the pieces'
.. edges on the face, not on its edges) each give nu_b mid_b, nu its normal
.. in the face (the piece's normal projected) times its length; plus the
.. length on the face's edge u_b = 1, from the traces' ends on it.
*/
static double face_measure (const Pieces * S, int a, int s, int * known)
{
  Crossing p[CROSSINGS_MAX];
  int np = 0;
#if D == 2
  int b = 1 - a;
  for (int i = 0; i < S->n; i++) {
    const Piece * P = &S->p[i];
    if (P->on_face)
      continue;
    for (int k = 0; k < 2; k++)
      if ((P->f[k] & bit (a, s)) && !(P->f[k] & (bit (b, 0) | bit (b, 1))) &&
          sign (P->A[b]) && np < CROSSINGS_MAX)
        p[np++] = (Crossing) { P->u[k][b], sign (P->A[b]) };
  }
  return line_measure (p, np, known);
#elif D == 1
  (void) S, (void) a, (void) s, (void) p, (void) np;
  *known = 0;                          /* no manifold of dimension 0 here */
  return 0.;
#else
  int b = (a + 1) % 3, c = (a + 2) % 3, traces = 0;
  double v = 0;
  for (int i = 0; i < S->n; i++) {
    const Piece * P = &S->p[i];
    if (P->on_face)
      continue;
    for (int k = 0; k < P->n; k++) {
      int l = (k + 1) % P->n, fe = P->f[k] & P->f[l];
      if (!(fe & bit (a, s)) || (fe & ~bit (a, s)))
        continue;
      /* nu l: (t_c, -t_b) in the face, turned to the piece's normal */
      double tb = P->u[l][b] - P->u[k][b], tc = P->u[l][c] - P->u[k][c];
      double nb = tc, nc = -tb, dot = nb * P->A[b] + nc * P->A[c];
      if (dot == 0)
        continue;
      if (dot < 0)
        nb = -nb, nc = -nc;
      v += nb * 0.5 * (P->u[k][b] + P->u[l][b]);
      traces++;
      /* its ends on the face's edge u_b = 1, along c */
      int ends[2] = { k, l };
      for (int e = 0; e < 2; e++) {
        int q = ends[e];
        if ((P->f[q] & bit (b, 1)) && !(P->f[q] & (bit (c, 0) | bit (c, 1)))
            && sign (nc) && np < CROSSINGS_MAX)
          p[np++] = (Crossing) { P->u[q][c], sign (nc) };
      }
    }
  }
  int kl;
  double L = line_measure (p, np, &kl);
  *known = traces > 0 || kl;
  if (kl)
    v += L;
  else if (v < -WRAP_TOL)              /* the edge was in the phase: F_v */
    v += 1;
  return v < 0 ? 0 : v > 1 ? 1 : v;
#endif
}

/* the fraction of the box: the pieces along x, plus its face x = 1 */
static double box_measure (const Pieces * S, int * known)
{
  double v = 0;
  int pieces = 0;
  for (int i = 0; i < S->n; i++)
    if (!S->p[i].on_face) {
      v += S->p[i].A[0] * S->p[i].c[0];
      pieces++;
    }
  if (!pieces) {
    *known = 0;
    return 0.;
  }
  int kf;
  double F = face_measure (S, 0, 1, &kf);
  if (kf)
    v += F;
  else if (v < -WRAP_TOL)              /* the face was in the phase: F_v */
    v += 1;
  *known = 1;
  return v < 0 ? 0 : v > 1 ? 1 : v;
}

static Pieces scratch;

double hm_vof_box (const Hmesh * M, const double * lo, double h,
                   const uint32_t * elem, uint32_t n, int * known)
{
  pieces_of (M, lo, h, elem, n, &scratch);
  return box_measure (&scratch, known);
}

double hm_vof_face (const Hmesh * M, const double * lo, double h,
                    const uint32_t * elem, uint32_t n, int a, int s,
                    int * known)
{
  pieces_of (M, lo, h, elem, n, &scratch);
  return face_measure (&scratch, a, s, known);
}

/* -------- Algorithm 4: the side of the closest element -------- */

static double dot3 (const double * a, const double * b)
{
  return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

/* the point of the segment ab closest to p */
static void closest_segment (const double * p, const double * a,
                             const double * b, double * q)
{
  double ab[3], ap[3];
  for (int d = 0; d < 3; d++)
    ab[d] = b[d] - a[d], ap[d] = p[d] - a[d];
  double t = dot3 (ap, ab) / dot3 (ab, ab);
  t = t < 0 ? 0 : t > 1 ? 1 : t;
  for (int d = 0; d < 3; d++)
    q[d] = a[d] + t * ab[d];
}

/* the point of the triangle abc closest to p (Ericson, Real-Time Collision
   Detection, 2005, 5.1.5) */
static void closest_triangle (const double * p, const double * a,
                              const double * b, const double * c, double * q)
{
  double ab[3], ac[3], ap[3], bp[3], cp[3];
  for (int d = 0; d < 3; d++) {
    ab[d] = b[d] - a[d], ac[d] = c[d] - a[d], ap[d] = p[d] - a[d];
    bp[d] = p[d] - b[d], cp[d] = p[d] - c[d];
  }
  double d1 = dot3 (ab, ap), d2 = dot3 (ac, ap);
  if (d1 <= 0 && d2 <= 0) {
    memcpy (q, a, 3 * sizeof (double));
    return;
  }
  double d3 = dot3 (ab, bp), d4 = dot3 (ac, bp);
  if (d3 >= 0 && d4 <= d3) {
    memcpy (q, b, 3 * sizeof (double));
    return;
  }
  double vc = d1*d4 - d3*d2;
  if (vc <= 0 && d1 >= 0 && d3 <= 0) {
    double v = d1 / (d1 - d3);
    for (int d = 0; d < 3; d++)
      q[d] = a[d] + v * ab[d];
    return;
  }
  double d5 = dot3 (ab, cp), d6 = dot3 (ac, cp);
  if (d6 >= 0 && d5 <= d6) {
    memcpy (q, c, 3 * sizeof (double));
    return;
  }
  double vb = d5*d2 - d1*d6;
  if (vb <= 0 && d2 >= 0 && d6 <= 0) {
    double w = d2 / (d2 - d6);
    for (int d = 0; d < 3; d++)
      q[d] = a[d] + w * ac[d];
    return;
  }
  double va = d3*d6 - d5*d4;
  if (va <= 0 && (d4 - d3) >= 0 && (d5 - d6) >= 0) {
    double w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
    for (int d = 0; d < 3; d++)
      q[d] = b[d] + w * (c[d] - b[d]);
    return;
  }
  double den = 1. / (va + vb + vc), v = vb * den, w = vc * den;
  for (int d = 0; d < 3; d++)
    q[d] = a[d] + ab[d] * v + ac[d] * w;
}

/* the element's normal, out of the phase (not unit) */
static void element_normal (const Hmesh * M, uint32_t e, double * n)
{
  int m = M->dim;
  double v[3][3] = { { 0 } };
  for (int i = 0; i <= m; i++)
    for (int d = 0; d < D; d++)
      v[i][d] = M->x[d][hm_vertex (M, m, e, i)];
  if (m == 1) {
    n[0] = v[1][1] - v[0][1], n[1] = - (v[1][0] - v[0][0]), n[2] = 0;
    return;
  }
  double e1[3], e2[3];
  for (int d = 0; d < 3; d++)
    e1[d] = v[1][d] - v[0][d], e2[d] = v[2][d] - v[0][d];
  n[0] = e1[1]*e2[2] - e1[2]*e2[1];
  n[1] = e1[2]*e2[0] - e1[0]*e2[2];
  n[2] = e1[0]*e2[1] - e1[1]*e2[0];
}

int hm_vof_inside (const Hmesh * M, const double * x, const uint32_t * elem,
                   uint32_t n)
{
  int m = M->dim, in = -1;
  double dmin = HUGE_VAL, cbest = -1, p[3] = { 0, 0, 0 };
  for (int d = 0; d < D; d++)
    p[d] = x[d];
  for (uint32_t i = 0; i < n; i++) {
    double v[3][3] = { { 0 } }, q[3], r[3], nn[3];
    for (int k = 0; k <= m; k++)
      for (int d = 0; d < D; d++)
        v[k][d] = M->x[d][hm_vertex (M, m, elem[i], k)];
    if (m == 1)
      closest_segment (p, v[0], v[1], q);
    else
      closest_triangle (p, v[0], v[1], v[2], q);
    for (int d = 0; d < 3; d++)
      r[d] = p[d] - q[d];
    double dist = sqrt (dot3 (r, r));
    element_normal (M, elem[i], nn);
    double nl = sqrt (dot3 (nn, nn));
    /* among the closest (to round-off), the one x faces most squarely: at
       an edge or a corner shared by several, the side is theirs */
    double cs = dist > 0 && nl > 0 ? fabs (dot3 (r, nn)) / (dist * nl) : 1;
    int closer = dist < dmin * (1 - 1e-12);
    if (closer || (dist <= dmin * (1 + 1e-12) && cs > cbest)) {
      dmin = closer ? dist : fmin (dmin, dist);
      cbest = cs;
      in = dot3 (r, nn) <= 0;          /* x not on the side of the normal */
    }
  }
  if (in < 0)
    return -1;
  return dmin == 0 ? 1 : in;
}

/* -------- on the tree -------- */

typedef struct {
  const Hmesh * M;
  int        m;                  /* cells a leaf side */
  uint32_t   ncell;              /* cells a leaf: m^D */
  double     origin[3], size;
  uint32_t   nl;
  double   (* lo)[3];            /* each leaf's lower corner, */
  double   * H;                  /* its side, */
  int      * level;
  uint64_t * key;                /* its key's prefix (hm_morton_cell ()) */
  double     hmin;               /* the finest leaf's side */
  uint32_t * lstart, * lelem;    /* each leaf's elements */
  double   (* emin)[3], (* emax)[3];  /* each element's box, padded */
  double   * frac;
  uint8_t  * state;              /* 0 unknown, 1 cut (frac exact),
                                    2 given by a neighbour */
} Grid;

static int in_domain (const Grid * G, const double * x)
{
  for (int d = 0; d < D; d++)
    if (!(x[d] >= G->origin[d] && x[d] <= G->origin[d] + G->size))
      return 0;
  return 1;
}

/* the leaf x is in, by its key among the leaves' (in Morton order) */
static uint32_t locate (const Grid * G, const double * x)
{
  if (!in_domain (G, x) || !G->nl)
    return HM_NONE;
  uint64_t k = hm_morton_key (x, G->origin, G->size);
  uint32_t a = 0, b = G->nl;          /* the last leaf whose key <= k */
  while (b - a > 1) {
    uint32_t c = (a + b) / 2;
    if (G->key[c] <= k)
      a = c;
    else
      b = c;
  }
  int shift = D * (HM_MORTON_BITS - G->level[a]);
  if (G->key[a] > k ||
      (shift < 64 && k >> shift != G->key[a] >> shift))
    return HM_NONE;                    /* a pruned part */
  return a;
}

/* the cell of leaf q x is in */
static uint32_t cell_of (const Grid * G, uint32_t q, const double * x,
                         int * ijk)
{
  double delta = G->H[q] / G->m;
  uint32_t c = 0, stride = 1;
  for (int d = 0; d < D; d++) {
    int i = (int) floor ((x[d] - G->lo[q][d]) / delta);
    i = i < 0 ? 0 : i >= G->m ? G->m - 1 : i;
    if (ijk)
      ijk[d] = i;
    c += (uint32_t) i * stride;
    stride *= (uint32_t) G->m;
  }
  return c;
}

/* the leaves a box meets: halved until each part's corners are in one leaf
   (or a part is below half the finest leaf, when its corners' leaves are
   all it meets) */
typedef struct {
  uint32_t * q;
  uint32_t   n, cap;
} List;

static void list_add (List * L, uint32_t q)
{
  for (uint32_t i = 0; i < L->n; i++)
    if (L->q[i] == q)
      return;
  if (L->n == L->cap) {
    L->cap = L->cap ? 2 * L->cap : 16;
    L->q = grow (L->q, L->cap, sizeof (uint32_t));
  }
  L->q[L->n++] = q;
}

static void cover (const Grid * G, const double * bmin, const double * bmax,
                   List * L)
{
  uint32_t q0 = HM_NONE, qs[8];
  int same = 1, small = 1;
  for (int c = 0; c < 1 << D; c++) {
    double x[3];
    for (int d = 0; d < D; d++)
      x[d] = c >> d & 1 ? bmax[d] : bmin[d];
    qs[c] = locate (G, x);
    if (c == 0)
      q0 = qs[c];
    else if (qs[c] != q0)
      same = 0;
  }
  for (int d = 0; d < D; d++)
    small &= bmax[d] - bmin[d] <= 0.5 * G->hmin;
  if (same || small) {
    for (int c = 0; c < 1 << D; c++)
      if (qs[c] != HM_NONE)
        list_add (L, qs[c]);
    return;
  }
  for (int c = 0; c < 1 << D; c++) {
    double lo[3], hi[3];
    for (int d = 0; d < D; d++) {
      double mid = 0.5 * (bmin[d] + bmax[d]);
      lo[d] = c >> d & 1 ? mid : bmin[d];
      hi[d] = c >> d & 1 ? bmax[d] : mid;
    }
    cover (G, lo, hi, L);
  }
}

/* the elements of leaf q whose box meets the box [lo, lo + h]^D */
static uint32_t gather (const Grid * G, uint32_t q, const double * lo,
                        double h, uint32_t ** out, uint32_t * cap)
{
  uint32_t n = 0;
  for (uint32_t i = G->lstart[q]; i < G->lstart[q + 1]; i++) {
    uint32_t e = G->lelem[i];
    int meets = 1;
    for (int d = 0; d < D; d++)
      meets &= G->emax[e][d] >= lo[d] && G->emin[e][d] <= lo[d] + h;
    if (!meets)
      continue;
    if (n == *cap) {
      *cap = *cap ? 2 * *cap : 64;
      *out = grow (*out, *cap, sizeof (uint32_t));
    }
    (*out)[n++] = e;
  }
  return n;
}

static void cell_box (const Grid * G, uint32_t q, uint32_t c, double * lo,
                      double * h)
{
  *h = G->H[q] / G->m;
  for (int d = 0; d < D; d++) {
    lo[d] = G->lo[q][d] + (double) (c % (uint32_t) G->m) * *h;
    c /= (uint32_t) G->m;
  }
}

/* the cell across the face (a, s) of cell c of leaf q, or HM_NONE */
static int neighbour (const Grid * G, uint32_t q, uint32_t c, int a, int s,
                      uint32_t * qn, uint32_t * cn)
{
  double lo[3], h, x[3];
  cell_box (G, q, c, lo, &h);
  for (int d = 0; d < D; d++)
    x[d] = lo[d] + 0.5 * h;
  x[a] += (s ? 0.75 : -0.75) * h;
  *qn = locate (G, x);
  if (*qn == HM_NONE)
    return 0;
  *cn = cell_of (G, *qn, x, NULL);
  return 1;
}

/*
.. The state of the unknown cell (qu, cu) from the cut cell (qk, ck) across
.. its face (a, s) (s the side of the unknown cell's face): the fraction of
.. that face, on the cut cell's side, which no element crosses; at a level
.. jump the part of it the unknown cell has (a box of its size in the cut
.. cell). -1 if neither that nor the closest element tells.
*/
static int state_across (const Grid * G, uint32_t qk, uint32_t ck,
                         uint32_t qu, uint32_t cu, int a, int s,
                         uint32_t ** buf, uint32_t * cap)
{
  double lok[3], hk, lou[3], hu, lo[3], h;
  cell_box (G, qk, ck, lok, &hk);
  cell_box (G, qu, cu, lou, &hu);
  int side;                            /* of the box's face on the cell u */
  if (hk <= hu) {
    memcpy (lo, lok, sizeof (lo));
    h = hk;
    side = 1 - s;
  }
  else {
    memcpy (lo, lou, sizeof (lo));
    h = hu;
    lo[a] += s ? h : -h;
    side = 1 - s;
  }
  uint32_t n = gather (G, qk, lo, h, buf, cap);
  int known;
  double f = hm_vof_face (G->M, lo, h, *buf, n, a, side, &known);
  if (known)
    return f > 0.5;
  double x[3];
  for (int d = 0; d < D; d++)
    x[d] = lo[d] + 0.5 * h;
  x[a] = lo[a] + side * h;
  return hm_vof_inside (G->M, x, *buf, n);
}

/* the fractions, read back in a foreach (): the leaf and cell of each
   centre */
static const Grid * writing;

static double fraction_at (double x, double y, double z)
{
  double p[3] = { x, y, z };
  uint32_t q = locate (writing, p);
  if (q == HM_NONE)
    return 0.;
  return writing->frac[(size_t) q * writing->ncell +
                       cell_of (writing, q, p, NULL)];
}

static void write_scalar (scalar c)
{
#if D == 1
  foreach ()
    c[] = (Real) fraction_at (x, 0., 0.);
#elif D == 2
  foreach ()
    c[] = (Real) fraction_at (x, y, 0.);
#else
  foreach ()
    c[] = (Real) fraction_at (x, y, z);
#endif
}

/* -------- the grid's leaves, the elements binned to them -------- */

static void grid_leaves (Grid * G, Hmesh * M)
{
  memset (G, 0, sizeof (*G));
  G->M = M;
  G->m = tree_get_block_m ();
  G->ncell = 1;
  for (int d = 0; d < D; d++)
    G->ncell *= (uint32_t) G->m;
  for (int d = 0; d < 3; d++)
    G->origin[d] = d < D ? tree_get_origin (d) : 0.;
  G->size = tree_get_size ();
  const TreeCaches * tc = tree_caches ();
  G->nl = (uint32_t) tc->nleaves;
  G->lo = grow (NULL, G->nl, sizeof (*G->lo));
  G->H = grow (NULL, G->nl, sizeof (double));
  G->level = grow (NULL, G->nl, sizeof (int));
  G->key = grow (NULL, G->nl, sizeof (uint64_t));
  G->hmin = G->size;
  for (uint32_t q = 0; q < G->nl; q++) {
    Index p = tc->leaves[q].p;
    int ijk[3] = { p.i, 0, 0 };
#if D > 1
    ijk[1] = p.j;
#endif
#if D > 2
    ijk[2] = p.k;
#endif
    G->level[q] = p.level;
    G->H[q] = G->size / (double) (1L << p.level);
    for (int d = 0; d < D; d++)
      G->lo[q][d] = G->origin[d] + ijk[d] * G->H[q];
    G->key[q] = hm_morton_cell (ijk, p.level);
    G->hmin = fmin (G->hmin, G->H[q]);
  }
}

/* each element to the leaves its box meets: lstart, lelem */
static void grid_bin (Grid * G)
{
  const Hmesh * M = G->M;
  int m = M->dim;
  uint32_t ne = hm_end (M, m), npair = 0, cappair = 0;
  uint32_t (* pair)[2] = NULL;
  G->emin = grow (NULL, ne, sizeof (*G->emin));
  G->emax = grow (NULL, ne, sizeof (*G->emax));
  double pad = HM_VOF_EPS * G->hmin / G->m;
  List L = { NULL, 0, 0 };
  hm_foreach (M, m, e) {
    int out = 0;
    for (int d = 0; d < D; d++) {
      double lo = HUGE_VAL, hi = -HUGE_VAL;
      for (int i = 0; i <= m; i++) {
        double x = M->x[d][hm_vertex (M, m, e, i)];
        lo = fmin (lo, x), hi = fmax (hi, x);
      }
      G->emin[e][d] = lo - pad, G->emax[e][d] = hi + pad;
      out |= G->emax[e][d] < G->origin[d] ||
        G->emin[e][d] > G->origin[d] + G->size;
    }
    if (out)
      continue;
    double bmin[3], bmax[3];
    for (int d = 0; d < D; d++) {
      bmin[d] = fmax (G->emin[e][d], G->origin[d]);
      bmax[d] = fmin (G->emax[e][d], G->origin[d] + G->size);
    }
    L.n = 0;
    cover (G, bmin, bmax, &L);
    for (uint32_t i = 0; i < L.n; i++) {
      if (npair == cappair) {
        cappair = cappair ? 2 * cappair : 1024;
        pair = grow (pair, cappair, sizeof (*pair));
      }
      pair[npair][0] = L.q[i];
      pair[npair++][1] = e;
    }
  }
  free (L.q);
  G->lstart = calloc ((size_t) G->nl + 1, sizeof (uint32_t));
  G->lelem = grow (NULL, npair, sizeof (uint32_t));
  for (uint32_t i = 0; i < npair; i++)
    G->lstart[pair[i][0] + 1]++;
  for (uint32_t q = 0; q < G->nl; q++)
    G->lstart[q + 1] += G->lstart[q];
  uint32_t * fill = grow (NULL, G->nl, sizeof (uint32_t));
  memcpy (fill, G->lstart, G->nl * sizeof (uint32_t));
  for (uint32_t i = 0; i < npair; i++)
    G->lelem[fill[pair[i][0]]++] = pair[i][1];
  free (fill);
  free (pair);
}

static void grid_unbin (Grid * G)
{
  free (G->lstart), free (G->lelem), free (G->emin), free (G->emax);
  G->lstart = G->lelem = NULL;
  G->emin = G->emax = NULL;
}

static void grid_free (Grid * G)
{
  grid_unbin (G);
  free (G->lo), free (G->H), free (G->level), free (G->key);
  free (G->frac), free (G->state);
}

/* every cell any element's box meets: visit (G, q, cell, its elements) */
typedef void (* CellVisit) (Grid * G, uint32_t q, uint32_t cc,
                            const uint32_t * elem, uint32_t n, void * ctx);

static void grid_cells (Grid * G, CellVisit visit, void * ctx)
{
  uint32_t * cstart = grow (NULL, G->ncell + 1, sizeof (uint32_t));
  uint32_t * celem = NULL, ccap = 0;
  for (uint32_t q = 0; q < G->nl; q++) {
    if (G->lstart[q] == G->lstart[q + 1])
      continue;
    double delta = G->H[q] / G->m;
    memset (cstart, 0, (G->ncell + 1) * sizeof (uint32_t));
    for (int pass = 0; pass < 2; pass++) {
      for (uint32_t i = G->lstart[q]; i < G->lstart[q + 1]; i++) {
        uint32_t e = G->lelem[i];
        int r0[3] = { 0, 0, 0 }, r1[3] = { 0, 0, 0 };
        for (int d = 0; d < D; d++) {
          int a = (int) floor ((G->emin[e][d] - G->lo[q][d]) / delta);
          int b = (int) floor ((G->emax[e][d] - G->lo[q][d]) / delta);
          r0[d] = a < 0 ? 0 : a >= G->m ? G->m - 1 : a;
          r1[d] = b < 0 ? 0 : b >= G->m ? G->m - 1 : b;
        }
        for (int k = r0[2]; k <= r1[2]; k++)
          for (int j = r0[1]; j <= r1[1]; j++)
            for (int i2 = r0[0]; i2 <= r1[0]; i2++) {
              uint32_t cc = (uint32_t) ((k * G->m + j) * G->m + i2);
              if (pass == 0)
                cstart[cc + 1]++;
              else
                celem[cstart[cc]++] = e;
            }
      }
      if (pass == 0) {
        for (uint32_t cc = 0; cc < G->ncell; cc++)
          cstart[cc + 1] += cstart[cc];
        if (cstart[G->ncell] > ccap) {
          ccap = cstart[G->ncell];
          celem = grow (celem, ccap, sizeof (uint32_t));
        }
      }
    }
    /* after the second pass cstart[cc] is the end of cell cc's list */
    for (uint32_t cc = 0, begin = 0; cc < G->ncell; cc++) {
      uint32_t end = cstart[cc];
      if (end > begin)
        visit (G, q, cc, celem + begin, end - begin, ctx);
      begin = end;
    }
  }
  free (cstart);
  free (celem);
}

/* -------- hm_vof (): every cell -------- */

static void visit_exact (Grid * G, uint32_t q, uint32_t cc,
                         const uint32_t * elem, uint32_t n, void * ctx)
{
  (void) ctx;
  double lo[3], h;
  int known;
  cell_box (G, q, cc, lo, &h);
  double f = hm_vof_box (G->M, lo, h, elem, n, &known);
  if (known) {
    G->frac[(size_t) q * G->ncell + cc] = f;
    G->state[(size_t) q * G->ncell + cc] = 1;
  }
}

long hm_vof (Hmesh * M, scalar c)
{
  if (M->dim != D - 1 || D < 2)
    return -1;
  Grid G;
  grid_leaves (&G, M);
  grid_bin (&G);
  size_t ncells = (size_t) G.nl * G.ncell;
  G.frac = calloc (ncells ? ncells : 1, sizeof (double));
  G.state = calloc (ncells ? ncells : 1, sizeof (uint8_t));
  grid_cells (&G, visit_exact, NULL);

  /*
  .. The unknown cells: from a cut neighbour, across a face no element
  .. crosses, then from cell to cell. Both ways round, so that a level jump
  .. is crossed whichever side finds the other.
  */
  uint32_t * buf = NULL, bcap = 0;
  uint64_t * queue = grow (NULL, ncells, sizeof (uint64_t));
  size_t head = 0, tail = 0;
  for (size_t i = 0; i < ncells; i++) {
    uint32_t q = (uint32_t) (i / G.ncell), cc = (uint32_t) (i % G.ncell);
    for (int a = 0; a < D; a++)
      for (int s = 0; s < 2; s++) {
        uint32_t qn, cn;
        if (!neighbour (&G, q, cc, a, s, &qn, &cn))
          continue;
        size_t j = (size_t) qn * G.ncell + cn;
        int st = -1;
        if (G.state[i] == 1 && G.state[j] == 0)        /* this cut */
          st = state_across (&G, q, cc, qn, cn, a, 1 - s, &buf, &bcap);
        if (G.state[i] == 0 && G.state[j] == 1) {      /* that cut */
          st = state_across (&G, qn, cn, q, cc, a, s, &buf, &bcap);
          if (st >= 0) {
            G.frac[i] = st, G.state[i] = 2;
            queue[tail++] = i;
          }
          continue;
        }
        if (st >= 0) {
          G.frac[j] = st, G.state[j] = 2;
          queue[tail++] = j;
        }
      }
  }
  while (head < tail) {
    size_t i = queue[head++];
    uint32_t q = (uint32_t) (i / G.ncell), cc = (uint32_t) (i % G.ncell);
    for (int a = 0; a < D; a++)
      for (int s = 0; s < 2; s++) {
        uint32_t qn, cn;
        if (!neighbour (&G, q, cc, a, s, &qn, &cn))
          continue;
        size_t j = (size_t) qn * G.ncell + cn;
        if (G.state[j] == 0) {
          G.frac[j] = G.frac[i], G.state[j] = 2;
          queue[tail++] = j;
        }
      }
  }
  long left = 0;
  for (size_t i = 0; i < ncells; i++)
    left += G.state[i] == 0;
  free (queue);
  free (buf);

  writing = &G;
  write_scalar (c);
  writing = NULL;
  grid_free (&G);
  return left;
}

/* -------- manifold_advect () -------- */

/* the value of cell cc (interior, i fastest) of leaf q in a scalar's block:
   the leaf's node from the tree's list, its block, the cell's index in it;
   no lookup in the tree's level index */
/* leaf q is this process's: its blocks are its own (under MPI another
   leaf may have none, or ghosts' copies) */
static int own_leaf (uint32_t q)
{
  return tree_caches ()->leaves[q].node->pid == tree_rank ();
}

static Real * cell_value (scalar s, uint32_t q, uint32_t cc, int m,
                          BlockShape sh)
{
  int x[3] = { 0, 0, 0 };
  for (int d = 0; d < D; d++) {
    x[d] = (int) (cc % (uint32_t) m);
    cc /= (uint32_t) m;
  }
  return block_data (tree_caches ()->leaves[q].node->s, s) +
    block_at (sh, x);
}

/*
.. Component a of a face vector at the point s (in cells of leaf q, from its
.. lower corner): linear along a between the faces i and i + 1 (face m in
.. the ghost layer), and along each other axis between the cell centres j
.. and j + 1 (j from -1, a ghost, to m - 1). Exact for a linear field.
*/
static double face_interp (scalar f, uint32_t q, const double * s, int a,
                           int m, BlockShape sh)
{
  const Real * b = block_data (tree_caches ()->leaves[q].node->s, f);
  int i0[3] = { 0, 0, 0 };
  double w[3] = { 0, 0, 0 };
  for (int d = 0; d < D; d++) {
    double t = d == a ? s[d] : s[d] - 0.5;
    int lo = d == a ? 0 : -1;
    int i = (int) floor (t);
    i = i < lo ? lo : i > m - 1 ? m - 1 : i;
    i0[d] = i;
    w[d] = t - i;
  }
  double v = 0;
  for (int corner = 0; corner < 1 << D; corner++) {
    int x[3] = { 0, 0, 0 };
    double wt = 1;
    for (int d = 0; d < D; d++) {
      int up = corner >> d & 1;
      x[d] = i0[d] + up;
      wt *= up ? w[d] : 1 - w[d];
    }
    if (wt != 0)
      v += wt * b[block_at (sh, x)];
  }
  return v;
}

/* the velocity of the face vector (components comp) at x in leaf q */
static void velocity (const Grid * G, const int * comp, uint32_t q,
                      const double * x, BlockShape sh, double * u)
{
  double s[3] = { 0, 0, 0 }, delta = G->H[q] / G->m;
  for (int d = 0; d < D; d++)
    s[d] = (x[d] - G->lo[q][d]) / delta;
  for (int a = 0; a < D; a++)
    u[a] = face_interp (comp[a], q, s, a, G->m, sh);
}

/* the velocity to move by at x (leaf q): uf's, or extrapolated in time
   from uf and uf_prev */
static void velocity_at (const Grid * G, const int * cu, const int * cp,
                         int second, double c1, double c0, uint32_t q,
                         const double * x, BlockShape sh, double * u)
{
  velocity (G, cu, q, x, sh, u);
  if (second) {
    double up[3] = { 0, 0, 0 };
    velocity (G, cp, q, x, sh, up);
    for (int d = 0; d < D; d++)
      u[d] = c1 * u[d] + c0 * up[d];
  }
}

/* x moved by s u: along a face of the box it is on, stopped on the box */
static void step (const Grid * G, const double * x, double s,
                  const double * u, double * xn)
{
  for (int d = 0; d < D; d++) {
    double lo = G->origin[d], hi = G->origin[d] + G->size;
    xn[d] = x[d] == lo || x[d] == hi ? x[d] : x[d] + s * u[d];
    xn[d] = xn[d] < lo ? lo : xn[d] > hi ? hi : xn[d];
  }
}

/* a set of cells by key q ncell + cc, each with a fraction and a state */
typedef struct {
  uint64_t * key;            /* key + 1, 0 empty */
  double   * f;
  uint8_t  * st;             /* 0 a candidate, 1 cut (f exact) */
  size_t     cap, n;
} CellMap;

static size_t cellmap_slot (CellMap * C, uint64_t key, int insert)
{
  if (insert && 2 * (C->n + 1) > C->cap) {
    CellMap N = { NULL, NULL, NULL, C->cap ? 2 * C->cap : 1024, 0 };
    N.key = calloc (N.cap, sizeof (uint64_t));
    N.f = grow (NULL, N.cap, sizeof (double));
    N.st = grow (NULL, N.cap, sizeof (uint8_t));
    if (!N.key) {
      fprintf (stderr, "hmesh: out of memory\n");
      abort ();
    }
    for (size_t i = 0; i < C->cap; i++)
      if (C->key[i]) {
        size_t j = (C->key[i] * 0x9E3779B97F4A7C15ull) & (N.cap - 1);
        while (N.key[j])
          j = (j + 1) & (N.cap - 1);
        N.key[j] = C->key[i], N.f[j] = C->f[i], N.st[j] = C->st[i];
        N.n++;
      }
    free (C->key), free (C->f), free (C->st);
    *C = N;
  }
  if (!C->cap)
    return SIZE_MAX;
  uint64_t k1 = key + 1;
  size_t j = (k1 * 0x9E3779B97F4A7C15ull) & (C->cap - 1);
  while (C->key[j] && C->key[j] != k1)
    j = (j + 1) & (C->cap - 1);
  if (C->key[j])
    return j;
  if (!insert)
    return SIZE_MAX;
  C->key[j] = k1, C->f[j] = 0, C->st[j] = 0;
  C->n++;
  return j;
}

static void cellmap_free (CellMap * C)
{
  free (C->key), free (C->f), free (C->st);
}

static void visit_candidate (Grid * G, uint32_t q, uint32_t cc,
                             const uint32_t * elem, uint32_t n, void * ctx)
{
  (void) elem, (void) n;
  cellmap_slot (ctx, (uint64_t) q * G->ncell + cc, 1);
}

static void visit_cut (Grid * G, uint32_t q, uint32_t cc,
                       const uint32_t * elem, uint32_t n, void * ctx)
{
  CellMap * C = ctx;
  double lo[3], h;
  int known;
  cell_box (G, q, cc, lo, &h);
  double f = hm_vof_box (G->M, lo, h, elem, n, &known);
  size_t j = cellmap_slot (C, (uint64_t) q * G->ncell + cc, 1);
  if (known)
    C->f[j] = f, C->st[j] = 1;
}

/* the component ids of a face vector */
static void components (face_vector f, int * comp)
{
  comp[0] = f.x;
#if D > 1
  comp[1] = f.y;
#endif
#if D > 2
  comp[2] = f.z;
#endif
}

struct HmAdvect {
  Hmesh    * M;
  Grid       G;
  CellMap    old;
  double  (* dx)[3];
  int        status;
};

HmAdvect * hm_advect_begin (Hmesh * M, face_vector uf, double dt,
                            face_vector uf_prev, double dt_prev,
                            const char * which)
{
  HmAdvect * A = calloc (1, sizeof (HmAdvect));
  if (!A) {
    fprintf (stderr, "hmesh: out of memory\n");
    abort ();
  }
  A->M = M;
  const HmTreeCache * TC;
  if (M->dim != D - 1 || D < 2 || !(dt > 0) || !(TC = hm_tree_cache (M)) ||
      TC->outside) {
    A->status = -2;
    return A;
  }
  int second = dt_prev > 0, cu[3] = { 0, 0, 0 }, cp[3] = { 0, 0, 0 };
  components (uf, cu);
  ghosts_update_list (D, cu, 0);
  if (second) {
    components (uf_prev, cp);
    ghosts_update_list (D, cp, 0);
  }
  BlockShape sh = tree_block_shape ();
  Grid * G = &A->G;
  grid_leaves (G, M);

  /*
  .. Each point's displacement, from the velocities at its leaf (the point
  .. cache's merge walk: no search). First order: x + dt u (x). Second
  .. order: the midpoint rule with the velocity extrapolated to t + dt / 2,
  .. u' = (1 + dt / 2 dt_prev) u - dt / 2 dt_prev u_prev, x' = x + dt / 2
  .. u' (x), x + dt u' (x'). (Adams-Bashforth would need each point's own
  .. velocity at its previous position: u_prev read where the point is now
  .. is u for a steady flow, and that is first order.) A point on the box
  .. moves along it; one that would leave it stops on it. None moves if one
  .. would go HM_ADVECT_CFL of the cell it leaves or enters, or further.
  */
  double c1 = 1, c0 = 0;
  if (second)
    c1 = 1 + dt / (2 * dt_prev), c0 = - dt / (2 * dt_prev);
  uint32_t end = hm_end (M, 0);
  double (* dx)[3] = grow (NULL, end, sizeof (*dx));
  int refused = 0;
  hm_foreach (M, 0, v) {
    for (int d = 0; d < 3; d++)
      dx[v][d] = 0.;
    if (which && !which[v])
      continue;
    uint32_t q = TC->leaf_of[v];
    double x[3] = { 0, 0, 0 }, u[3] = { 0, 0, 0 };
    for (int d = 0; d < D; d++)
      x[d] = M->x[d][v];
    velocity_at (G, cu, cp, second, c1, c0, q, x, sh, u);
    if (second) {
      /* the midpoint, and its leaf: this one, or found */
      double xm[3] = { 0, 0, 0 };
      step (G, x, 0.5 * dt, u, xm);
      uint32_t qm = q;
      for (int d = 0; d < D; d++)
        if (xm[d] < G->lo[q][d] || xm[d] >= G->lo[q][d] + G->H[q])
          qm = HM_NONE;
      if (qm == HM_NONE && (qm = hm_tree_locate (M, xm)) == HM_NONE) {
        refused = -2;
        break;
      }
      velocity_at (G, cu, cp, second, c1, c0, qm, xm, sh, u);
    }
    double xn[3] = { 0, 0, 0 }, l = 0;
    step (G, x, dt, u, xn);
    for (int d = 0; d < D; d++) {
      dx[v][d] = xn[d] - x[d];
      l += dx[v][d] * dx[v][d];
    }
    double h = fmin (G->H[q] / G->m, hm_tree_cell_size (M, xn));
    if (!(sqrt (l) < HM_ADVECT_CFL * h)) {
      refused = -1;
      break;
    }
  }
  A->dx = dx;
  if (refused) {
    A->status = refused;
    return A;
  }

  /* the cells the interface may cut before the move */
  grid_bin (G);
  grid_cells (G, visit_candidate, &A->old);
  grid_unbin (G);
  return A;
}

int hm_advect_status (const HmAdvect * A)
{
  return A->status;
}

void hm_advect_free (HmAdvect * A)
{
  if (!A)
    return;
  free (A->dx);
  cellmap_free (&A->old);
  grid_free (&A->G);                  /* zeros if never made: nothing */
  free (A);
}

void hm_advect_move (HmAdvect * A)
{
  Hmesh * M = A->M;
  if (A->status)
    return;
  hm_foreach (M, 0, v)
    for (int d = 0; d < D; d++)
      M->x[d][v] += A->dx[v][d];
  hm_touch (M);
}

int hm_advect_end (HmAdvect * A, scalar c)
{
  if (A->status) {
    int st = A->status;
    hm_advect_free (A);
    return st;
  }
  Hmesh * M = A->M;
  BlockShape sh = tree_block_shape ();
  Grid G = A->G;
  CellMap old = A->old, now = { NULL, NULL, NULL, 0, 0 };

  /* the cells cut after it: exact */
  grid_bin (&G);
  grid_cells (&G, visit_cut, &now);
  int m = G.m;
  for (size_t j = 0; j < now.cap; j++)
    if (now.key[j] && now.st[j] == 1 &&
        own_leaf ((uint32_t) ((now.key[j] - 1) / G.ncell))) {
      uint64_t k = now.key[j] - 1;
      *cell_value (c, (uint32_t) (k / G.ncell), (uint32_t) (k % G.ncell), m,
                   sh) = (Real) now.f[j];
    }

  /*
  .. The cells not cut now that may have been, or are touched without a cut:
  .. an old c of 0 or 1 is kept; else the exact fraction of a face shared
  .. with a cut cell, which no element crosses; else the old c rounded.
  .. Every other cell is left as it was.
  */
  uint32_t * buf = NULL, bcap = 0;
  for (int pass = 0; pass < 2; pass++) {
    CellMap * S = pass ? &now : &old;
    for (size_t j = 0; j < S->cap; j++) {
      if (!S->key[j] || (pass && S->st[j] == 1))
        continue;
      uint64_t k = S->key[j] - 1;
      if (!pass) {
        size_t i = cellmap_slot (&now, k, 0);
        if (i != SIZE_MAX)
          continue;                    /* in now: the second pass */
      }
      uint32_t q = (uint32_t) (k / G.ncell), cc = (uint32_t) (k % G.ncell);
      if (!own_leaf (q))
        continue;
      Real * cv = cell_value (c, q, cc, m, sh);
      double o = *cv;
      if (o <= 0 || o >= 1)
        continue;
      int st = -1;
      for (int a = 0; a < D && st < 0; a++)
        for (int s = 0; s < 2 && st < 0; s++) {
          uint32_t qn, cn;
          if (!neighbour (&G, q, cc, a, s, &qn, &cn))
            continue;
          size_t jn = cellmap_slot (&now, (uint64_t) qn * G.ncell + cn, 0);
          if (jn != SIZE_MAX && now.st[jn] == 1)
            st = state_across (&G, qn, cn, q, cc, a, s, &buf, &bcap);
        }
      *cv = (Real) (st >= 0 ? st : o > 0.5);
    }
  }
  free (buf);
  cellmap_free (&old);
  cellmap_free (&now);
  ghosts_mark_dirty (c);

  /* the edges now longer than their cells */
  int longer = 0;
  hm_foreach (M, 1, e) {
    double a[3] = { 0, 0, 0 }, b[3] = { 0, 0, 0 }, l = 0;
    for (int d = 0; d < D; d++) {
      a[d] = M->x[d][hm_vertex (M, 1, e, 0)];
      b[d] = M->x[d][hm_vertex (M, 1, e, 1)];
      l += (b[d] - a[d]) * (b[d] - a[d]);
    }
    longer += !(sqrt (l) < fmin (hm_tree_cell_size (M, a),
                                 hm_tree_cell_size (M, b)));
  }
  free (A->dx);
  free (A);
  grid_free (&G);
  return longer;
}

int manifold_advect (Hmesh * M, face_vector uf, scalar c, double dt,
                     face_vector uf_prev, double dt_prev)
{
  HmAdvect * A = hm_advect_begin (M, uf, dt, uf_prev, dt_prev, NULL);
  int st = hm_advect_status (A);
  if (st) {
    hm_advect_free (A);
    return st;
  }
  hm_advect_move (A);
  return hm_advect_end (A, c);
}
