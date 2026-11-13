/*
.. libhmesh's tests, in the euclidean dimension it is built with:
..   the pool (ids reused, holes skipped, compaction's map)
..   a curve with a junction, a triangulated square, a closed surface (3D),
..   a cube of tetrahedra (3D): counts, Euler characteristic, cycles, signs,
..   d d = 0, measures, the one-ring by next and twin
..   removals and additions at random, then compaction: every invariant
..   (hm_check ()), and the fields moved with their simplices
..   what is refused
..   the time of a removal and an addition
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

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

static uint64_t seed = 88172645463325252ull;

static uint32_t rnd (uint32_t n)
{
  seed ^= seed << 13;
  seed ^= seed >> 7;
  seed ^= seed << 17;
  return (uint32_t) (seed % n);
}

static uint32_t point (Hmesh * M, double x, double y, double z)
{
  double p[3] = { x, y, z };
  return hm_vertex_add (M, p);
}

/* the length of the twin cycle from h */
static uint32_t cycle (const Hmesh * M, int k, uint32_t h)
{
  uint32_t n = 0, g = h;
  do {
    g = hm_half_twin (M, k, g);
    n++;
  } while (g != h);
  return n;
}

static long euler (const Hmesh * M)
{
  long chi = 0;
  for (int k = 0; k <= M->dim; k++)
    chi += (k & 1 ? -1 : 1) * (long) hm_count (M, k);
  return chi;
}

/* d d = 0: for each simplex, the signs summed per face of codimension 2 */
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

/* -------- the pool -------- */

static void test_pool (void)
{
  Hmesh * M = hmesh_new (0);
  enum { N = 1000 };
  uint32_t id[N];
  int alive[N];
  for (int i = 0; i < N; i++) {
    id[i] = point (M, i, 0, 0);
    alive[i] = 1;
    CHECK (id[i] == (uint32_t) i, "a new point's id %u, not %d", id[i], i);
  }
  int n = N;
  for (int r = 0; r < N / 2; r++) {
    uint32_t i = rnd (N);
    if (alive[i]) {
      CHECK (hm_remove (M, 0, i) == 0, "point %u not removed", i);
      alive[i] = 0;
      n--;
    }
    else
      CHECK (hm_remove (M, 0, i) == -1, "point %u removed twice", i);
  }
  CHECK (hm_count (M, 0) == (uint32_t) n, "count %u, not %d",
         hm_count (M, 0), n);
  CHECK (!hm_dense (M, 0), "dense with holes");
  int seen = 0, inorder = 1;
  uint32_t last = 0;
  hm_foreach (M, 0, s) {
    CHECK (alive[s], "visited %u, removed", s);
    if (seen && s <= last)
      inorder = 0;
    last = s;
    seen++;
  }
  CHECK (seen == n && inorder, "visited %d of %d, in order %d", seen, n,
         inorder);
  CHECK (hm_check (M) == 0, "invariants");

  /* reuse: the last freed first */
  uint32_t top = M->level[0].freed[M->level[0].nfreed - 1];
  uint32_t a = point (M, top, 0, 0);
  CHECK (a == top, "reused %u, not the last freed %u", a, top);
  alive[a] = 1;
  n++;

  /* compaction: dense, in the order of ids, coordinates moved along */
  uint32_t * map[HM_LEVELS];
  hm_compact (M, map);
  CHECK (hm_dense (M, 0) && hm_count (M, 0) == (uint32_t) n, "not dense");
  uint32_t next = 0;
  for (int i = 0; i < N; i++)
    if (alive[i]) {
      CHECK (map[0][i] == next, "map %d -> %u, not %u", i, map[0][i], next);
      CHECK (hm_x (M, 0)[map[0][i]] == i, "point %d not moved along", i);
      next++;
    }
    else
      CHECK (map[0][i] == HM_NONE, "a removed point mapped");
  free (map[0]);
  CHECK (hm_check (M) == 0, "invariants after compaction");
  CHECK (point (M, 0, 0, 0) == (uint32_t) n, "after compaction, not the end");
  hmesh_free (M);
}

/* -------- meshes -------- */

/* a chain of n edges along x, and a branch at its middle point */
static Hmesh * curve (int n, uint32_t * junction)
{
  Hmesh * M = hmesh_new (1);
  for (int i = 0; i <= n; i++)
    point (M, i, 0, 0);
  for (uint32_t i = 0; i < (uint32_t) n; i++) {
    uint32_t v[2] = { i, i + 1 };
    hm_add (M, 1, v);
  }
  *junction = (uint32_t) n / 2;
  uint32_t p = point (M, *junction, 1, 1);
  uint32_t v[2] = { *junction, p };
  hm_add (M, 1, v);
  return M;
}

/* the surface z = (x + y)/2 over the unit square (z only in 3D), n x n
   squares in two triangles each, counterclockwise seen from above */
#if HM_DIMENSION >= 2
static Hmesh * square (int n)
{
  Hmesh * M = hmesh_new (2);
  for (int j = 0; j <= n; j++)
    for (int i = 0; i <= n; i++)
      point (M, (double) i / n, (double) j / n, 0.5*(i + j)/n);
  for (uint32_t j = 0; j < (uint32_t) n; j++)
    for (uint32_t i = 0; i < (uint32_t) n; i++) {
      uint32_t w = (uint32_t) n + 1, a = j*w + i, b = a + 1, c = a + w,
        d = c + 1;
      uint32_t t0[3] = { a, b, d }, t1[3] = { a, d, c };
      hm_add (M, 2, t0);
      hm_add (M, 2, t1);
    }
  return M;
}
#endif

#if HM_DIMENSION == 3
static Hmesh * octahedron (void)
{
  Hmesh * M = hmesh_new (2);
  double p[6][3] = {{1,0,0},{-1,0,0},{0,1,0},{0,-1,0},{0,0,1},{0,0,-1}};
  for (int i = 0; i < 6; i++)
    point (M, p[i][0], p[i][1], p[i][2]);
  /* outward: (a, b, c) with (b - a) x (c - a) . a > 0 */
  uint32_t t[8][3] = {{0,2,4},{2,1,4},{1,3,4},{3,0,4},
                      {2,0,5},{1,2,5},{3,1,5},{0,3,5}};
  for (int i = 0; i < 8; i++)
    hm_add (M, 2, t[i]);
  return M;
}

/* the unit cube, n^3 cubes in six tetrahedra each (Kuhn), each of positive
   volume */
static Hmesh * cube (int n)
{
  Hmesh * M = hmesh_new (3);
  uint32_t w = (uint32_t) n + 1;
  for (int k = 0; k <= n; k++)
    for (int j = 0; j <= n; j++)
      for (int i = 0; i <= n; i++)
        point (M, (double) i / n, (double) j / n, (double) k / n);
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
            v[0] = v[1];
            v[1] = x;
            t = hm_add (M, 3, v);
          }
        }
  return M;
}
#endif

/* -------- the meshes' structure -------- */

static void test_curve (void)
{
  uint32_t jn;
  Hmesh * M = curve (10, &jn);
  CHECK (hm_check (M) == 0, "curve: invariants");
  CHECK (hm_count (M, 0) == 12 && hm_count (M, 1) == 11, "curve: %u, %u",
         hm_count (M, 0), hm_count (M, 1));
  CHECK (euler (M) == 1, "curve: Euler %ld", euler (M));
  hm_foreach (M, 0, v) {
    uint32_t n = cycle (M, 1, hm_use (M, 0, v));
    uint32_t want = v == jn ? 3 : v == 0 || v == 10 || v == 11 ? 1 : 2;
    CHECK (n == want, "curve: %u edges at point %u, not %u", n, v, want);
    /* the boundary of the whole chain, at each point: +1 at its end */
    int sum = 0;
    uint32_t h = hm_use (M, 0, v), g = h;
    do {
      sum += hm_half_sign (M, 1, g);
      g = hm_half_twin (M, 1, g);
    } while (g != h);
    int want_sum = v == 0 ? -1 : v == 10 || v == 11 ? 1 : v == jn ? -1 : 0;
    CHECK (sum == want_sum, "curve: d at point %u is %d, not %d", v, sum,
           want_sum);
  }
  /* an edge's two ends: next, and the sign of each */
  uint32_t e = hm_half_simplex (1, hm_use (M, 0, 3));
  uint32_t h0 = hm_half (1, e, 0), h1 = hm_half_next (1, h0);
  CHECK (hm_half_next (1, h1) == h0 && hm_half_facet (M, 1, h0) ==
         hm_vertex (M, 1, e, 1) && hm_half_sign (M, 1, h0) == 1 &&
         hm_half_sign (M, 1, h1) == -1, "curve: an edge's ends");
  double len = 0;
  hm_foreach (M, 1, s)
    len += hm_measure (M, 1, s);
  double want = HM_DIMENSION == 1 ? 10 : 10 + (HM_DIMENSION == 2 ? 1 :
                                                 sqrt (2.));
  CHECK (fabs (len - want) < 1e-12, "curve: length %g, not %g", len, want);
  hmesh_free (M);
}

static void test_square (void)
{
#if HM_DIMENSION >= 2
  int n = 8;
  Hmesh * M = square (n);
  uint32_t V = (uint32_t) ((n + 1)*(n + 1)), E = (uint32_t) (3*n*n + 2*n),
    F = (uint32_t) (2*n*n);
  CHECK (hm_check (M) == 0, "square: invariants");
  CHECK (hm_count (M, 0) == V && hm_count (M, 1) == E && hm_count (M, 2) == F,
         "square: %u %u %u", hm_count (M, 0), hm_count (M, 1),
         hm_count (M, 2));
  CHECK (euler (M) == 1, "square: Euler %ld", euler (M));
  CHECK (boundary_of_boundary (M) == 0, "square: d d is not 0");

  /* each edge: two triangles of opposite signs inside, one on the boundary */
  uint32_t boundary = 0;
  hm_foreach (M, 1, e) {
    uint32_t h = hm_use (M, 1, e), n2 = cycle (M, 2, h);
    if (n2 == 1) {
      boundary++;
      CHECK (hm_half_boundary (M, 2, h), "square: boundary");
    }
    else {
      uint32_t t = hm_half_twin (M, 2, h);
      CHECK (n2 == 2 && hm_half_twin (M, 2, t) == h, "square: twin");
      CHECK (hm_half_sign (M, 2, h) == -hm_half_sign (M, 2, t),
             "square: an edge's two halves of the same sign");
      /* opposite directions: h from a to b, t from b to a */
      CHECK (hm_half_origin (M, h) ==
             hm_half_vertex (M, 2, hm_half_next (2, hm_half_next (2, t))),
             "square: twins not opposite");
    }
  }
  CHECK (boundary == (uint32_t) (4*n), "square: %u boundary edges",
         boundary);

  /* the one-ring of a point inside, by twin (prev (h)) */
  uint32_t c = (uint32_t) (n/2*(n + 1) + n/2), h = HM_NONE;
  int brute = 0;
  hm_foreach (M, 2, t)
    for (int i = 0; i < 3; i++)
      if (hm_vertex (M, 2, t, i) == c) {
        brute++;
        h = hm_half (2, t, (i + 2) % 3);    /* from c */
      }
  int ring = 0;
  uint32_t g = h;
  do {
    CHECK (hm_half_origin (M, g) == c, "square: the ring left its point");
    g = hm_half_twin (M, 2, hm_half_next (2, hm_half_next (2, g)));
    ring++;
  } while (g != h && ring < 100);
  CHECK (ring == brute && ring == 6, "square: a ring of %d, %d at the point",
         ring, brute);

  double area = 0, signed_area = 0;
  hm_foreach (M, 2, t) {
    area += hm_measure (M, 2, t);
    signed_area += hm_volume_signed (M, t);
  }
  double want = HM_DIMENSION == 2 ? 1 : sqrt (1.5);
  CHECK (fabs (area - want) < 1e-12, "square: area %.15g, not %.15g", area,
         want);
  CHECK (fabs (signed_area - want) < 1e-12, "square: signed area %.15g",
         signed_area);
  hmesh_free (M);
#endif
}

static void test_closed (void)
{
#if HM_DIMENSION == 3
  Hmesh * M = octahedron ();
  CHECK (hm_check (M) == 0, "octahedron: invariants");
  CHECK (hm_count (M, 1) == 12 && euler (M) == 2, "octahedron: %u edges, "
         "Euler %ld", hm_count (M, 1), euler (M));
  hm_foreach (M, 1, e) {
    uint32_t h = hm_use (M, 1, e), t = hm_half_twin (M, 2, h);
    CHECK (cycle (M, 2, h) == 2 && hm_half_sign (M, 2, h) ==
           -hm_half_sign (M, 2, t), "octahedron: edge %u", e);
  }
  CHECK (boundary_of_boundary (M) == 0, "octahedron: d d is not 0");
  hmesh_free (M);

  int n = 3;
  M = cube (n);
  uint32_t T = (uint32_t) (6*n*n*n);
  CHECK (hm_check (M) == 0, "cube: invariants");
  CHECK (hm_count (M, 3) == T && hm_count (M, 0) == (uint32_t)
         ((n + 1)*(n + 1)*(n + 1)), "cube: %u tetrahedra", hm_count (M, 3));
  CHECK (euler (M) == 1, "cube: Euler %ld", euler (M));
  CHECK (boundary_of_boundary (M) == 0, "cube: d d is not 0");
  uint32_t boundary = 0;
  hm_foreach (M, 2, f) {
    uint32_t h = hm_use (M, 2, f), c = cycle (M, 3, h);
    if (c == 1)
      boundary++;
    else
      CHECK (c == 2 && hm_half_sign (M, 3, h) ==
             -hm_half_sign (M, 3, hm_half_twin (M, 3, h)),
             "cube: face %u, %u tetrahedra", f, c);
  }
  CHECK (boundary == (uint32_t) (12*n*n), "cube: %u boundary faces",
         boundary);
  double vol = 0;
  hm_foreach (M, 3, t)
    vol += hm_volume_signed (M, t);
  CHECK (fabs (vol - 1) < 1e-12, "cube: volume %.15g", vol);
  /* the tetrahedra about an inner edge, through the cycles of its faces */
  hmesh_free (M);
#endif
}

/* -------- removals and additions, then compaction -------- */

/* a curve, a square or a cube, small (for the tests) or large (the time) */
static Hmesh * mesh (int m, int large)
{
  uint32_t jn;
  if (m == 1)
    return curve (large ? 1 << 17 : 40, &jn);
#if HM_DIMENSION >= 2
  if (m == 2)
    return square (large ? 256 : 12);
#endif
#if HM_DIMENSION == 3
  return cube (large ? 24 : 4);
#endif
  return NULL;
}

/* each level's scalar holds the measure of its simplex, the points' their
   first coordinate */
static void fill (Hmesh * M, const int * f, int all)
{
  for (int k = 0; k <= M->dim; k++)
    hm_foreach (M, k, s)
      if (all || hm_val (M, f[k], 0)[s] == 0)
        hm_val (M, f[k], 0)[s] = (HmReal) (k ? hm_measure (M, k, s) :
                                           1 + hm_x (M, 0)[s]);
}

static int filled (const Hmesh * M, const int * f)
{
  for (int k = 0; k <= M->dim; k++)
    hm_foreach (M, k, s) {
      double want = k ? hm_measure (M, k, s) : 1 + hm_x (M, 0)[s];
      if (fabs (hm_val (M, f[k], 0)[s] - want) > 1e-6*want)
        return -1;
    }
  return 0;
}

static void test_churn (int m)
{
  Hmesh * M = mesh (m, 0);
  int f[HM_LEVELS];
  for (int k = 0; k <= m; k++)
    f[k] = hm_scalar_new (M, k, "measure");
  int vec = hm_vector_new (M, m, "u");
  fill (M, f, 1);
  uint32_t count[HM_LEVELS];
  for (int k = 0; k <= m; k++)
    count[k] = hm_count (M, k);
  long chi = euler (M);

  /* half the top simplices out, with what they leave unused */
  uint32_t ntop = hm_count (M, m), nout = 0;
  uint32_t (* out)[4] = malloc (ntop * sizeof (*out));
  for (uint32_t r = 0; r < ntop / 2; r++) {
    uint32_t s = hm_alive_from (M, m, rnd (hm_end (M, m)));
    if (s == HM_NONE)
      continue;
    for (int i = 0; i <= m; i++)
      out[nout][i] = hm_vertex (M, m, s, i);
    nout++;
    CHECK (hm_remove_closure (M, m, s) == 0, "churn %d: not removed", m);
  }
  CHECK (hm_check (M) == 0, "churn %d: invariants after removals", m);
  CHECK (filled (M, f) == 0, "churn %d: fields after removals", m);

  /* nothing left unused: every point on an edge, every face on a coface */
  for (int k = 0; k < m; k++)
    hm_foreach (M, k, s)
      CHECK (hm_use (M, k, s) != HM_NONE, "churn %d: level %d, %u unused",
             m, k, s);
  CHECK (hm_count (M, m) == ntop - nout, "churn %d: %u left of %u - %u", m,
         hm_count (M, m), ntop, nout);
  hmesh_free (M);
  free (out);

  /* the same, keeping the points: a simplex out alone with its unused faces
     down to the edges, then in again */
  M = mesh (m, 0);
  for (int k = 0; k <= m; k++)
    f[k] = hm_scalar_new (M, k, "measure");
  vec = hm_vector_new (M, m, "u");
  fill (M, f, 1);
  out = malloc (ntop * sizeof (*out));
  nout = 0;
  for (uint32_t r = 0; r < ntop / 2; r++) {
    uint32_t s = hm_alive_from (M, m, rnd (hm_end (M, m)));
    if (s == HM_NONE)
      continue;
    for (int i = 0; i <= m; i++)
      out[nout][i] = hm_vertex (M, m, s, i);
    nout++;
    /* the faces of s, then those of them left unused, level by level */
    uint32_t face[HM_LEVELS][64];
    int nface[HM_LEVELS] = { 0 };
    face[m][nface[m]++] = s;
    for (int k = m; k >= 1; k--)
      for (int a = 0; a < nface[k]; a++) {
        uint32_t t = face[k][a];
        for (int i = 0; i <= k; i++) {
          uint32_t x = hm_half_facet (M, k, hm_half (k, t, i)), b = 0;
          while (b < (uint32_t) nface[k - 1] && face[k - 1][b] != x)
            b++;
          if (b == (uint32_t) nface[k - 1] && k - 1 >= 1)
            face[k - 1][nface[k - 1]++] = x;
        }
      }
    for (int k = m; k >= 1; k--)
      for (int a = 0; a < nface[k]; a++)
        if (k == m || hm_use (M, k, face[k][a]) == HM_NONE)
          CHECK (hm_remove (M, k, face[k][a]) == 0, "churn %d: level %d",
                 m, k);
  }
  CHECK (hm_check (M) == 0, "churn %d: invariants after removals", m);
  CHECK (hm_count (M, 0) == count[0], "churn %d: points removed", m);
  for (uint32_t r = 0; r < nout; r++) {
    uint32_t s = hm_add (M, m, out[r]);
    CHECK (s != HM_NONE, "churn %d: not added again", m);
    hm_val (M, vec, 0)[s] = 1;
  }
  fill (M, f, 0);
  CHECK (hm_check (M) == 0, "churn %d: invariants after additions", m);
  for (int k = 0; k <= m; k++)
    CHECK (hm_count (M, k) == count[k], "churn %d: level %d, %u, not %u", m,
           k, hm_count (M, k), count[k]);
  CHECK (euler (M) == chi, "churn %d: Euler %ld, not %ld", m, euler (M), chi);
  CHECK (filled (M, f) == 0, "churn %d: fields after additions", m);
  CHECK (boundary_of_boundary (M) == 0, "churn %d: d d is not 0", m);

  /* compaction: dense, the same structure, fields along */
  uint32_t marked = 0;
  hm_foreach (M, m, s)
    marked += hm_val (M, vec, 0)[s] == 1;
  hm_compact (M, NULL);
  CHECK (hm_check (M) == 0, "churn %d: invariants after compaction", m);
  for (int k = 0; k <= m; k++)
    CHECK (hm_dense (M, k) && hm_count (M, k) == count[k],
           "churn %d: level %d not dense", m, k);
  CHECK (filled (M, f) == 0, "churn %d: fields after compaction", m);
  uint32_t marked2 = 0;
  for (uint32_t s = 0; s < hm_count (M, m); s++)
    marked2 += hm_val (M, vec, 0)[s] == 1;
  CHECK (marked == marked2 && marked == nout, "churn %d: %u marked, %u, %u",
         m, marked, marked2, nout);
  CHECK (euler (M) == chi && boundary_of_boundary (M) == 0,
         "churn %d: after compaction", m);
  free (out);
  hmesh_free (M);
}

/* -------- what is refused -------- */

static void test_refused (void)
{
  CHECK (hmesh_new (HM_DIMENSION + 1) == NULL && hmesh_new (-1) == NULL,
         "a manifold beyond the dimension");
  CHECK (manifold (HM_DIMENSION + 1) == -1 && hmesh == NULL, "manifold ()");
  CHECK (manifold (1) == 0 && hmesh && hmesh->dim == 1, "manifold (1)");
  Hmesh * M = hmesh;
  uint32_t a = point (M, 0, 0, 0), b = point (M, 1, 0, 0),
    c = point (M, 2, 0, 0);
  uint32_t aa[2] = { a, a }, ab[2] = { a, b }, ba[2] = { b, a },
    bc[2] = { b, c }, abc[3] = { a, b, c }, ax[2] = { a, 77 };
  CHECK (hm_add (M, 1, aa) == HM_NONE, "an edge from a point to itself");
  CHECK (hm_add (M, 1, ax) == HM_NONE, "an edge to a point not there");
  CHECK (hm_add (M, 2, abc) == HM_NONE, "a triangle on a curve");
  uint32_t e = hm_add (M, 1, ab);
  CHECK (hm_add (M, 1, ba) == e && hm_find (M, 1, ba) == e &&
         hm_vertex (M, 1, e, 0) == a, "the same edge, its orientation kept");
  CHECK (hm_find (M, 1, bc) == HM_NONE, "an edge not there");
  CHECK (hm_remove (M, 0, a) == -1, "a point removed under its edge");
  CHECK (hm_remove (M, 1, e) == 0 && hm_remove (M, 0, a) == 0,
         "removed in order");
  CHECK (hm_check (M) == 0, "invariants");
  manifold_free ();
  CHECK (hmesh == NULL, "manifold_free ()");

#if HM_DIMENSION >= 2
  M = square (2);
  uint32_t edge = hm_half_facet (M, 2, hm_half (2, 0, 0));
  CHECK (hm_remove (M, 1, edge) == -1, "an edge removed under triangles");
  CHECK (hm_remove_closure (M, 1, edge) == -1, "the same, with its faces");
  CHECK (hm_check (M) == 0, "invariants");
  hmesh_free (M);
#endif
}

/* -------- hm_check () finds what is broken -------- */

static void corrupt (Hmesh * M, uint32_t * p, uint32_t bad, const char * what)
{
  uint32_t keep = *p;
  *p = bad;
  fprintf (stderr, "     (expected) ");
  CHECK (hm_check (M) == -1, "hm_check (): %s not found", what);
  *p = keep;
}

static void test_corrupt (void)
{
  int m = HM_DIMENSION;
  Hmesh * M = mesh (m, 0);
  HmLevel * l = &M->level[m], * lf = &M->level[m - 1];
  CHECK (hm_check (M) == 0, "corrupt: invariants");
  /* a half inside, h, and its twin */
  uint32_t h = 0;
  while (hm_half_boundary (M, m, h))
    h++;
  uint32_t t = l->twin[h], f = hm_half_facet (M, m, h);
  corrupt (M, &l->twin[h], h, "a twin cut short");
  corrupt (M, &l->twin[h], hm_half_next (m, h), "a twin on another facet");
  corrupt (M, &l->facet[h], l->facet[h] ^ HM_FLIP, "a wrong sign");
  corrupt (M, &l->facet[h], hm_half_facet (M, m, hm_half_next (m, h)),
           "a wrong facet");
  corrupt (M, &l->vert[h], l->vert[t], "a wrong vertex");
  corrupt (M, &lf->use[f], hm_half_next (m, h), "a use on another facet");
  corrupt (M, &l->n, l->n + 1, "a count");
  {
    uint64_t keep = l->alive[0];
    l->alive[0] &= ~(uint64_t) 1;
    fprintf (stderr, "     (expected) ");
    CHECK (hm_check (M) == -1, "hm_check (): a bit cleared not found");
    l->alive[0] = keep;
  }
  CHECK (hm_check (M) == 0, "corrupt: invariants restored");
  hmesh_free (M);
}

/* -------- the time of a removal and an addition -------- */

static double now (void)
{
  struct timespec t;
  clock_gettime (CLOCK_MONOTONIC, &t);
  return (double) t.tv_sec + 1e-9*(double) t.tv_nsec;
}

static void bench (void)
{
  int m = HM_DIMENSION;
  Hmesh * M = mesh (m, 1);
  hm_scalar_new (M, m, "s");
  enum { R = 200000 };
  double t0 = now ();
  for (int r = 0; r < R; r++) {
    uint32_t s = hm_alive_from (M, m, rnd (hm_end (M, m))), v[4];
    if (s == HM_NONE)
      continue;
    for (int i = 0; i <= m; i++)
      v[i] = hm_vertex (M, m, s, i);
    hm_remove (M, m, s);
    hm_add (M, m, v);
  }
  double t = now () - t0;
  CHECK (hm_check (M) == 0, "bench: invariants");
  printf ("a %d-simplex out and in again: %.0f ns (of %u)\n", m, 1e9*t/R,
          hm_count (M, m));
  hmesh_free (M);

  /* the pool alone: points with two fields, out at random and in again */
  M = hmesh_new (0);
  for (uint32_t i = 0; i < 1u << 20; i++)
    point (M, i, 0, 0);
  hm_scalar_new (M, 0, "s");
  hm_vector_new (M, 0, "u");
  double p[3] = { 0, 0, 0 };
  t0 = now ();
  for (int r = 0; r < R; r++) {
    hm_remove (M, 0, rnd (1u << 20));
    hm_vertex_add (M, p);
  }
  t = now () - t0;
  printf ("a point out and in again: %.0f ns (of %u)\n", 1e9*t/R,
          hm_count (M, 0));
  hmesh_free (M);
}

int main (void)
{
  test_pool ();
  test_curve ();
  test_square ();
  test_closed ();
  for (int m = 1; m <= HM_DIMENSION; m++)
    test_churn (m);
  test_refused ();
  test_corrupt ();
  bench ();
  if (fails)
    printf ("%d failed\n", fails);
  else
    printf ("passed, %dD\n", HM_DIMENSION);
  return fails != 0;
}
