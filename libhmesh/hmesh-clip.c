/*
.. A manifold cut to a box: see hm_clip () in hmesh.h.
..
.. References
..   I. E. Sutherland, G. W. Hodgman, Reentrant polygon clipping, Commun.
..     ACM 17 (1), 1974: one plane at a time, a polygon's vertices in order,
..     each kept if in, a new one where an edge crosses.
..   J. Dompierre, P. Labbe, M.-G. Vallet, R. Camarero, How to subdivide
..     pyramids, prisms and hexahedra into tetrahedra, 8th International
..     Meshing Roundtable, 1999: the faces' diagonals chosen from the
..     lowest vertex of a global order, so that neighbours agree, and a
..     prism always splits into three tetrahedra without a new point.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "hmesh.h"

static void * grow (void * p, size_t n, size_t e)
{
  void * q = realloc (p, (n ? n : 1) * e);
  if (!q) {
    fprintf (stderr, "hmesh: out of memory\n");
    abort ();
  }
  return q;
}

/* a point of the part inside of a simplex: its id, and its barycentric
   coordinates in that simplex */
typedef struct {
  uint32_t id;
  double   b[4];
} Cv;

/* a plane x[axis] = p, the inside x > p (lower face) or x < p (upper) */
typedef struct {
  int        axis, upper;
  double     p;
  int8_t   * side;          /* by point: 1 in, 0 on the plane, -1 out */
  uint32_t * cut;           /* by edge: its new point, or HM_NONE */
} Plane;

/* a new simplex: level k, its points, the simplex it comes from */
typedef struct {
  int      k;
  uint32_t v[4], parent;
} Piece;

typedef struct {
  Piece  * piece;
  uint32_t npiece, cappiece;
  uint32_t (* rem)[2];      /* level, id: the simplices that go */
  uint32_t nrem, caprem;
} Work;

static void piece_push (Work * W, int k, const Cv * const * v,
                        uint32_t parent)
{
  if (W->npiece == W->cappiece) {
    W->cappiece = W->cappiece ? 2 * W->cappiece : 64;
    W->piece = grow (W->piece, W->cappiece, sizeof (Piece));
  }
  Piece * P = &W->piece[W->npiece++];
  P->k = k;
  P->parent = parent;
  for (int i = 0; i <= k; i++)
    P->v[i] = v[i]->id;
}

static void rem_push (Work * W, int k, uint32_t s)
{
  if (W->nrem == W->caprem) {
    W->caprem = W->caprem ? 2 * W->caprem : 64;
    W->rem = grow (W->rem, W->caprem, sizeof (*W->rem));
  }
  W->rem[W->nrem][0] = (uint32_t) k;
  W->rem[W->nrem++][1] = s;
}

static int side (const Plane * P, const Cv * c)
{
  return P->side[c->id];
}

/* the point where the edge from a to b, one in and one out, crosses */
static Cv cut_point (const Hmesh * M, const Plane * P, const Cv * a,
                     const Cv * b)
{
  uint32_t ab[2] = { a->id, b->id };
  Cv c = { P->cut[hm_find (M, 1, ab)], { 0, 0, 0, 0 } };
  double xa = M->x[P->axis][a->id], xb = M->x[P->axis][b->id];
  double t = (P->p - xa) / (xb - xa);
  for (int i = 0; i < 4; i++)
    c.b[i] = (1 - t)*a->b[i] + t*b->b[i];
  return c;
}

/* the polygon v [0 ... n - 1] cut by the plane, in order: each point in or
   on kept, a new one where an edge crosses (Sutherland and Hodgman) */
static int clip_polygon (const Hmesh * M, const Plane * P, const Cv * v,
                         int n, Cv * out)
{
  int m = 0;
  for (int i = 0; i < n; i++) {
    const Cv * a = &v[i], * b = &v[(i + 1) % n];
    if (side (P, a) >= 0)
      out[m++] = *a;
    if (side (P, a) * side (P, b) < 0)
      out[m++] = cut_point (M, P, a, b);
  }
  return m;
}

/* the triangles of the fan from the lowest point of the polygon */
static int fan (const Cv * v, int n, int (* tri)[3])
{
  int r = 0;
  for (int i = 1; i < n; i++)
    if (v[i].id < v[r].id)
      r = i;
  for (int i = 1; i + 1 < n; i++) {
    tri[i - 1][0] = r;
    tri[i - 1][1] = (r + i) % n;
    tri[i - 1][2] = (r + i + 1) % n;
  }
  return n > 2 ? n - 2 : 0;
}

/* the determinant of the barycentric coordinates of the k + 1 points:
   the piece's volume over its parent's, of the sign of its orientation */
static double det_bary (int k, const Cv * const * v)
{
  int n = k + 1;
  double a[4][4], d = 1;
  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++)
      a[i][j] = v[i]->b[j];
  for (int c = 0; c < n; c++) {
    int p = c;
    for (int r = c + 1; r < n; r++)
      if (fabs (a[r][c]) > fabs (a[p][c]))
        p = r;
    if (a[p][c] == 0.)
      return 0.;
    if (p != c) {
      for (int j = 0; j < n; j++) {
        double t = a[c][j];
        a[c][j] = a[p][j];
        a[p][j] = t;
      }
      d = -d;
    }
    d *= a[c][c];
    for (int r = c + 1; r < n; r++) {
      double f = a[r][c] / a[c][c];
      for (int j = c; j < n; j++)
        a[r][j] -= f*a[c][j];
    }
  }
  return d;
}

/* a piece, in its parent's orientation; none if it is flat */
static void piece_oriented (Work * W, int k, const Cv ** v, uint32_t parent)
{
  double d = det_bary (k, v);
  if (d == 0.)
    return;
  if (d < 0) {
    const Cv * t = v[0];
    v[0] = v[1];
    v[1] = t;
  }
  piece_push (W, k, v, parent);
}

/* the part inside of the k-simplex s, with points in and out, as pieces */
static void pieces (const Hmesh * M, const Plane * P, Work * W, int k,
                    uint32_t s)
{
  Cv v[4];
  for (int i = 0; i <= k; i++) {
    v[i].id = hm_vertex (M, k, s, i);
    for (int j = 0; j < 4; j++)
      v[i].b[j] = i == j;
  }

  if (k == 1) {
    /* the end out replaced by the cut */
    Cv c[2];
    for (int i = 0; i < 2; i++)
      c[i] = side (P, &v[i]) < 0 ? cut_point (M, P, &v[1 - i], &v[i]) : v[i];
    const Cv * p[2] = { &c[0], &c[1] };
    piece_oriented (W, 1, p, s);
    return;
  }

  if (k == 2) {
    /* a triangle or a quadrangle: the fan from its lowest point */
    Cv poly[6];
    int tri[4][3], n = clip_polygon (M, P, v, 3, poly), nt = fan (poly, n, tri);
    for (int t = 0; t < nt; t++) {
      const Cv * p[3] = { &poly[tri[t][0]], &poly[tri[t][1]],
                          &poly[tri[t][2]] };
      piece_oriented (W, 2, p, s);
    }
    return;
  }

#if HM_DIMENSION == 3
  /* k = 3: the polyhedron's faces, the tetrahedron's four cut and the cap
     on the plane, each by its fan; then the cone from the lowest point to
     the triangles of the faces it is not on */
  Cv face[5][8], all[16];
  int nf[5], nface = 0, nall = 0;
  for (int i = 0; i < 4; i++) {
    Cv tri3[3];
    for (int j = 0, a = 0; j < 4; j++)
      if (j != i)
        tri3[a++] = v[j];
    nf[nface] = clip_polygon (M, P, tri3, 3, face[nface]);
    if (nf[nface] >= 3)
      nface++;
  }
  /* the cap: every point of the faces on the plane, once, in order about
     their centre in the plane's two axes */
  int n = 0;
  Cv * cap = face[nface];
  for (int f = 0; f < nface; f++)
    for (int a = 0; a < nf[f]; a++) {
      const Cv * c = &face[f][a];
      int seen = 0;
      for (int b = 0; b < nall; b++)
        seen |= all[b].id == c->id;
      if (!seen)
        all[nall++] = *c;
      if (side (P, c) == 0) {
        seen = 0;
        for (int b = 0; b < n; b++)
          seen |= cap[b].id == c->id;
        if (!seen && n < 8)
          cap[n++] = *c;
      }
    }
  if (n >= 3) {
    int u = (P->axis + 1) % 3, w = (P->axis + 2) % 3;
    double cu = 0, cw = 0, ang[8];
    for (int a = 0; a < n; a++)
      cu += M->x[u][cap[a].id] / n, cw += M->x[w][cap[a].id] / n;
    for (int a = 0; a < n; a++)
      ang[a] = atan2 (M->x[w][cap[a].id] - cw, M->x[u][cap[a].id] - cu);
    for (int a = 1; a < n; a++)        /* by angle: four points at most */
      for (int b = a; b > 0 && ang[b] < ang[b - 1]; b--) {
        double t = ang[b];
        ang[b] = ang[b - 1], ang[b - 1] = t;
        Cv c = cap[b];
        cap[b] = cap[b - 1], cap[b - 1] = c;
      }
    nf[nface++] = n;
  }
  int apex = 0;
  for (int a = 1; a < nall; a++)
    if (all[a].id < all[apex].id)
      apex = a;
  for (int f = 0; f < nface; f++) {
    int on = 0, tri[6][3];
    for (int a = 0; a < nf[f]; a++)
      on |= face[f][a].id == all[apex].id;
    if (on)                            /* flat with the apex */
      continue;
    int nt = fan (face[f], nf[f], tri);
    for (int t = 0; t < nt; t++) {
      const Cv * p[4] = { &all[apex], &face[f][tri[t][0]],
                          &face[f][tri[t][1]], &face[f][tri[t][2]] };
      piece_oriented (W, 3, p, s);
    }
  }
#endif
}

/* the fields of level k: s takes the values of parent */
static void fields_copy (Hmesh * M, int k, uint32_t s, uint32_t parent)
{
  for (int f = 0; f < M->nfield; f++)
    if (M->field[f].level == k) {
      size_t e = M->field[f].esize;
      for (int c = 0; c < M->field[f].ncomp; c++) {
        char * d = M->field[f].data[c];
        memmove (d + (size_t) s * e, d + (size_t) parent * e, e);
      }
    }
}

/* the points' fields at the new point c, t of the way from a to b */
static void fields_cut (Hmesh * M, uint32_t c, uint32_t a, uint32_t b,
                        double t)
{
  for (int f = 0; f < M->nfield; f++)
    if (M->field[f].level == 0)
      for (int i = 0; i < M->field[f].ncomp; i++) {
        if (M->field[f].real) {
          HmReal * d = M->field[f].data[i];
          d[c] = (HmReal) ((1 - t)*d[a] + t*d[b]);
        }
        else {
          size_t e = M->field[f].esize;
          char * d = M->field[f].data[i];
          memmove (d + (size_t) c * e, d + (size_t) (t < 0.5 ? a : b) * e, e);
        }
      }
}

static void clip_plane (Hmesh * M, int axis, int upper, double p, double snap)
{
  int m = M->dim;
  Plane P = { axis, upper, p, NULL, NULL };
  uint32_t end0 = hm_end (M, 0);
  P.side = grow (NULL, end0, sizeof (int8_t));
  hm_foreach (M, 0, v) {
    double x = M->x[axis][v];
    if (fabs (x - p) <= snap) {
      M->x[axis][v] = p;
      P.side[v] = 0;
    }
    else
      P.side[v] = (upper ? x < p : x > p) ? 1 : -1;
  }

  /* a point on each edge that crosses, from its lower id end */
  if (m >= 1) {
    uint32_t end1 = hm_end (M, 1), ncut = 0;
    P.cut = grow (NULL, end1, sizeof (uint32_t));
    uint32_t * cuts = grow (NULL, hm_count (M, 1), sizeof (uint32_t));
    for (uint32_t e = 0; e < end1; e++)
      P.cut[e] = HM_NONE;
    hm_foreach (M, 1, e) {
      uint32_t a = hm_vertex (M, 1, e, 0), b = hm_vertex (M, 1, e, 1);
      if (P.side[a] * P.side[b] < 0)
        cuts[ncut++] = e;
    }
    for (uint32_t i = 0; i < ncut; i++) {
      uint32_t e = cuts[i], a = hm_vertex (M, 1, e, 0),
        b = hm_vertex (M, 1, e, 1);
      if (b < a) {
        uint32_t t = a;
        a = b, b = t;
      }
      double xa = M->x[axis][a], xb = M->x[axis][b];
      double t = (p - xa) / (xb - xa), x[HM_DIMENSION];
      for (int d = 0; d < HM_DIMENSION; d++)
        x[d] = d == axis ? p : M->x[d][a] + t*(M->x[d][b] - M->x[d][a]);
      uint32_t c = hm_vertex_add (M, x);
      fields_cut (M, c, a, b, t);
      P.cut[e] = c;
      cuts[i] = c;
    }
    P.side = grow (P.side, hm_end (M, 0), sizeof (int8_t));
    for (uint32_t i = 0; i < ncut; i++)
      P.side[cuts[i]] = 0;
    free (cuts);
  }

  /* the maximal simplices, from the top: those across, their pieces, and
     those with nothing in, out */
  Work W = { NULL, 0, 0, NULL, 0, 0 };
  for (int k = m; k >= 1; k--)
    hm_foreach (M, k, s) {
      if (k < m && hm_use (M, k, s) != HM_NONE)
        continue;
      int in = 0, out = 0;
      for (int i = 0; i <= k; i++) {
        int sd = P.side[hm_vertex (M, k, s, i)];
        in += sd > 0, out += sd < 0;
      }
      if (!out)
        continue;
      if (in)
        pieces (M, &P, &W, k, s);
      rem_push (&W, k, s);
    }
  hm_foreach (M, 0, v)
    if ((m == 0 || hm_use (M, 0, v) == HM_NONE) && P.side[v] < 0)
      rem_push (&W, 0, v);

  /* the pieces in first, so that the faces they share with what goes stay */
  for (uint32_t i = 0; i < W.npiece; i++) {
    Piece * Q = &W.piece[i];
    uint32_t s = hm_add (M, Q->k, Q->v);
    if (s != HM_NONE)
      fields_copy (M, Q->k, s, Q->parent);
  }
  for (uint32_t i = 0; i < W.nrem; i++)
    if (hm_alive (M, (int) W.rem[i][0], W.rem[i][1]))
      hm_remove_closure (M, (int) W.rem[i][0], W.rem[i][1]);

  free (W.piece);
  free (W.rem);
  free (P.side);
  free (P.cut);
}

int hm_clip (Hmesh * M, const double * origin, double size)
{
  if (!(size > 0.))
    return -1;
  M->domain = 0;                       /* the new points are on it */
  double snap = HM_CLIP_SNAP * size;
  for (int d = 0; d < HM_DIMENSION; d++) {
    clip_plane (M, d, 0, origin[d], snap);
    clip_plane (M, d, 1, origin[d] + size, snap);
  }
  M->generation++;
  return hm_set_domain (M, origin, size);
}

int hm_domain_faces (const Hmesh * M, int k, uint32_t s)
{
  if (!M->domain)
    return 0;
  int faces = (1 << (2 * HM_DIMENSION)) - 1;
  for (int i = 0; i <= k; i++) {
    uint32_t v = hm_vertex (M, k, s, i);
    int on = 0;
    for (int d = 0; d < HM_DIMENSION; d++) {
      on |= (M->x[d][v] == M->origin[d]) << (2 * d);
      on |= (M->x[d][v] == M->origin[d] + M->size) << (2 * d + 1);
    }
    faces &= on;
  }
  return faces;
}

int hm_half_on_domain (const Hmesh * M, int k, uint32_t h)
{
  if (k < 1 || !hm_half_boundary (M, k, h))
    return 0;
  return hm_domain_faces (M, k - 1, hm_half_facet (M, k, h));
}
