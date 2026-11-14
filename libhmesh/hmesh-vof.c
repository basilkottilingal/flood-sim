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

long hm_vof (Hmesh * M, scalar c)
{
  if (M->dim != D - 1 || D < 2)
    return -1;
  Grid G = { .M = M, .m = tree_get_block_m () };
  G.ncell = 1;
  for (int d = 0; d < D; d++)
    G.ncell *= (uint32_t) G.m;
  for (int d = 0; d < 3; d++)
    G.origin[d] = d < D ? tree_get_origin (d) : 0.;
  G.size = tree_get_size ();

  /* the leaves */
  const TreeCaches * tc = tree_caches ();
  G.nl = (uint32_t) tc->nleaves;
  G.lo = grow (NULL, G.nl, sizeof (*G.lo));
  G.H = grow (NULL, G.nl, sizeof (double));
  G.level = grow (NULL, G.nl, sizeof (int));
  G.key = grow (NULL, G.nl, sizeof (uint64_t));
  G.hmin = G.size;
  for (uint32_t q = 0; q < G.nl; q++) {
    Index p = tc->leaves[q].p;
    int ijk[3] = { p.i, 0, 0 };
#if D > 1
    ijk[1] = p.j;
#endif
#if D > 2
    ijk[2] = p.k;
#endif
    G.level[q] = p.level;
    G.H[q] = G.size / (double) (1L << p.level);
    for (int d = 0; d < D; d++)
      G.lo[q][d] = G.origin[d] + ijk[d] * G.H[q];
    G.key[q] = hm_morton_cell (ijk, p.level);
    G.hmin = fmin (G.hmin, G.H[q]);
  }

  /* each element to the leaves its box meets, as (leaf, element) pairs */
  int m = M->dim;
  uint32_t ne = hm_end (M, m), npair = 0, cappair = 0;
  uint32_t (* pair)[2] = NULL;
  G.emin = grow (NULL, ne, sizeof (*G.emin));
  G.emax = grow (NULL, ne, sizeof (*G.emax));
  double pad = HM_VOF_EPS * G.hmin / G.m;
  List L = { NULL, 0, 0 };
  hm_foreach (M, m, e) {
    int out = 0;
    for (int d = 0; d < D; d++) {
      double lo = HUGE_VAL, hi = -HUGE_VAL;
      for (int i = 0; i <= m; i++) {
        double x = M->x[d][hm_vertex (M, m, e, i)];
        lo = fmin (lo, x), hi = fmax (hi, x);
      }
      G.emin[e][d] = lo - pad, G.emax[e][d] = hi + pad;
      out |= G.emax[e][d] < G.origin[d] ||
        G.emin[e][d] > G.origin[d] + G.size;
    }
    if (out)
      continue;
    double bmin[3], bmax[3];
    for (int d = 0; d < D; d++) {
      bmin[d] = fmax (G.emin[e][d], G.origin[d]);
      bmax[d] = fmin (G.emax[e][d], G.origin[d] + G.size);
    }
    L.n = 0;
    cover (&G, bmin, bmax, &L);
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
  G.lstart = calloc ((size_t) G.nl + 1, sizeof (uint32_t));
  G.lelem = grow (NULL, npair, sizeof (uint32_t));
  for (uint32_t i = 0; i < npair; i++)
    G.lstart[pair[i][0] + 1]++;
  for (uint32_t q = 0; q < G.nl; q++)
    G.lstart[q + 1] += G.lstart[q];
  uint32_t * fill = grow (NULL, G.nl, sizeof (uint32_t));
  memcpy (fill, G.lstart, G.nl * sizeof (uint32_t));
  for (uint32_t i = 0; i < npair; i++)
    G.lelem[fill[pair[i][0]]++] = pair[i][1];
  free (fill);
  free (pair);

  /* each leaf: its elements to its cells, and the cells with any cut */
  size_t ncells = (size_t) G.nl * G.ncell;
  G.frac = calloc (ncells ? ncells : 1, sizeof (double));
  G.state = calloc (ncells ? ncells : 1, sizeof (uint8_t));
  uint32_t * cstart = grow (NULL, G.ncell + 1, sizeof (uint32_t));
  uint32_t * celem = NULL, ccap = 0, * buf = NULL, bcap = 0;
  for (uint32_t q = 0; q < G.nl; q++) {
    if (G.lstart[q] == G.lstart[q + 1])
      continue;
    double delta = G.H[q] / G.m;
    memset (cstart, 0, (G.ncell + 1) * sizeof (uint32_t));
    for (int pass = 0; pass < 2; pass++) {
      for (uint32_t i = G.lstart[q]; i < G.lstart[q + 1]; i++) {
        uint32_t e = G.lelem[i];
        int r0[3] = { 0, 0, 0 }, r1[3] = { 0, 0, 0 };
        for (int d = 0; d < D; d++) {
          int a = (int) floor ((G.emin[e][d] - G.lo[q][d]) / delta);
          int b = (int) floor ((G.emax[e][d] - G.lo[q][d]) / delta);
          r0[d] = a < 0 ? 0 : a >= G.m ? G.m - 1 : a;
          r1[d] = b < 0 ? 0 : b >= G.m ? G.m - 1 : b;
        }
        for (int k = r0[2]; k <= r1[2]; k++)
          for (int j = r0[1]; j <= r1[1]; j++)
            for (int i2 = r0[0]; i2 <= r1[0]; i2++) {
              uint32_t cc = (uint32_t) ((k * G.m + j) * G.m + i2);
              if (pass == 0)
                cstart[cc + 1]++;
              else
                celem[cstart[cc]++] = e;
            }
      }
      if (pass == 0) {
        for (uint32_t cc = 0; cc < G.ncell; cc++)
          cstart[cc + 1] += cstart[cc];
        if (cstart[G.ncell] > ccap) {
          ccap = cstart[G.ncell];
          celem = grow (celem, ccap, sizeof (uint32_t));
        }
      }
    }
    /* after the second pass cstart[cc] is the end of cell cc's list */
    for (uint32_t cc = 0, begin = 0; cc < G.ncell; cc++) {
      uint32_t end = cstart[cc];
      if (end > begin) {
        double lo[3], h;
        int known;
        cell_box (&G, q, cc, lo, &h);
        double f = hm_vof_box (M, lo, h, celem + begin, end - begin, &known);
        if (known) {
          G.frac[(size_t) q * G.ncell + cc] = f;
          G.state[(size_t) q * G.ncell + cc] = 1;
        }
      }
      begin = end;
    }
  }
  free (cstart);
  free (celem);

  /*
  .. The unknown cells: from a cut neighbour, across a face no element
  .. crosses, then from cell to cell. Both ways round, so that a level jump
  .. is crossed whichever side finds the other.
  */
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

  free (G.lo), free (G.H), free (G.level), free (G.key);
  free (G.lstart), free (G.lelem), free (G.emin), free (G.emax);
  free (G.frac), free (G.state);
  return left;
}
