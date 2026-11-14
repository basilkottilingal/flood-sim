/*
.. libhmesh in libb-amr's tree (hmesh-tree.h), in the dimension it is built
.. with, on a tree refined over half its domain, the domain moved and
.. resized:
..   the box: points outside, NaN, on the upper faces; a tree too deep
..   the cache: each point in the one leaf whose box holds it (brute force,
..   faces shared by two leaves to the upper one), the leaves' ranges, the
..   keys sorted, duplicates; again after the points moved, after the tree
..   was refined, after the domain moved under the points
..   hm_compact_morton (): every level along the curve, fields moved, a
..   leaf's points the ids of its range
..   the time of the cache: sorted from nothing, again after small moves
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include "tree.h"
#include "block-ghosts.h"
#include "hmesh-tree.h"

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

static double urand (void)
{
  seed ^= seed << 13;
  seed ^= seed >> 7;
  seed ^= seed << 17;
  return (double) (seed >> 11) / 9007199254740992.;
}

static double now (void)
{
  struct timespec t;
  clock_gettime (CLOCK_MONOTONIC, &t);
  return (double) t.tv_sec + 1e-9*(double) t.tv_nsec;
}

/* the domain: not the unit cube */
static const double X0[3] = { -1.25, 2., 0.5 }, H0 = 3.;

static void point_at (double * x, const double * u)
{
  for (int d = 0; d < HM_DIMENSION; d++)
    x[d] = X0[d] + H0*u[d];
}

/* the leaves coarser than `to` with their centre in the lower half along x
   refined, as far as the 2:1 balance lets */
static void refine_half (int to)
{
  for (int pass = 0; pass < 4 * to; pass++) {
    const TreeCaches * tc = tree_caches ();
    int n = 0;
    Index * todo = malloc ((size_t) (tc->nleaves ? tc->nleaves : 1) *
                           sizeof (Index));
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

/* leaf q's box, [lo, hi] along each axis */
static void leaf_box (uint32_t q, double * lo, double * hi)
{
  Index p = tree_caches ()->leaves[q].p;
  int ijk[3] = { p.i, 0, 0 };
#if HM_DIMENSION > 1
  ijk[1] = p.j;
#endif
#if HM_DIMENSION > 2
  ijk[2] = p.k;
#endif
  double h = H0 / (double) (1 << p.level);
  for (int d = 0; d < HM_DIMENSION; d++) {
    lo[d] = X0[d] + ijk[d]*h;
    hi[d] = lo[d] + h;
  }
}

/* the leaf holding x, by brute force: [lo, hi) along each axis, [lo, hi]
   on the domain's upper face; HM_NONE if none, or more than one */
static uint32_t leaf_brute (const double * x)
{
  const TreeCaches * tc = tree_caches ();
  uint32_t found = HM_NONE, nfound = 0;
  for (uint32_t q = 0; q < (uint32_t) tc->nleaves; q++) {
    double lo[3], hi[3];
    leaf_box (q, lo, hi);
    int in = 1;
    for (int d = 0; d < HM_DIMENSION; d++) {
      int top = hi[d] >= X0[d] + H0;
      in = in && x[d] >= lo[d] && (x[d] < hi[d] || (top && x[d] <= hi[d]));
    }
    if (in)
      found = q, nfound++;
  }
  return nfound == 1 ? found : HM_NONE;
}

/* the cache against brute force */
static void check_cache (Hmesh * M, const char * when)
{
  const HmTreeCache * C = hm_tree_cache (M);
  CHECK (C, "%s: no cache", when);
  if (!C)
    return;
  CHECK (C->n == hm_count (M, 0), "%s: %u sorted, %u points", when, C->n,
         hm_count (M, 0));
  CHECK (C->nleaves == (uint32_t) tree_caches ()->nleaves && C->outside == 0,
         "%s: %u leaves, %u outside", when, C->nleaves, C->outside);
  /* sorted to the cells; the keys those of the points */
  int lo = HM_DIMENSION * (HM_MORTON_BITS - C->cell_bits);
  for (uint32_t a = 1; a < C->n; a++)
    if (C->key[a] >> lo < C->key[a - 1] >> lo) {
      CHECK (0, "%s: keys not sorted to the cells at %u", when, a);
      break;
    }
  for (uint32_t a = 0; a < C->n; a++) {
    double x[3];
    for (int d = 0; d < HM_DIMENSION; d++)
      x[d] = hm_x (M, d)[C->order[a]];
    if (C->key[a] != hm_morton_key (x, M->origin, M->size)) {
      CHECK (0, "%s: key %u not its point's", when, a);
      break;
    }
  }
  uint32_t in = 0, wrong = 0, last = 0;
  for (uint32_t q = 0; q < C->nleaves; q++) {
    CHECK (C->start[q] >= last && C->stop[q] >= C->start[q],
           "%s: leaf %u's range", when, q);
    last = C->stop[q];
    hm_foreach_leaf_point (C, q, v) {
      in++;
      double x[3];
      for (int d = 0; d < HM_DIMENSION; d++)
        x[d] = hm_x (M, d)[v];
      if (C->leaf_of[v] != q || leaf_brute (x) != q)
        wrong++;
    }
  }
  CHECK (in == C->n && wrong == 0, "%s: %u points in the leaves, %u in the "
         "wrong one", when, in, wrong);
  CHECK (hm_check (M) == 0, "%s: invariants", when);
}

static void test_tree (void)
{
  int depth = HM_DIMENSION == 3 ? 4 : 6;
  tree (depth, 4);
  tree_shift_resize (X0[0], X0[1], X0[2], H0);
  refine_half (depth - 1);

  Hmesh * M = hmesh_new (1);
  CHECK (hm_tree_bind (M) == 0, "bind");

  /* the box */
  double u[3], x[3];
  u[0] = u[1] = u[2] = 0.5;
  point_at (x, u);
  x[0] = X0[0] - 1e-9;
  CHECK (hm_vertex_add (M, x) == HM_NONE, "a point below the box");
  x[0] = X0[0] + H0 + 1e-9;
  CHECK (hm_vertex_add (M, x) == HM_NONE, "a point above the box");
  x[0] = NAN;
  CHECK (hm_vertex_add (M, x) == HM_NONE, "a NaN");

  /* points on the faces: the domain's corners, the middle (a level jump),
     leaf corners */
  for (int c = 0; c < 1 << HM_DIMENSION; c++) {
    for (int d = 0; d < HM_DIMENSION; d++)
      u[d] = c >> d & 1;
    point_at (x, u);
    CHECK (hm_vertex_add (M, x) != HM_NONE, "a corner of the box");
  }
  for (int d = 0; d < 3; d++)
    u[d] = 0.5;
  point_at (x, u);
  uint32_t mid = hm_vertex_add (M, x);
  for (int r = 0; r < 64; r++) {
    for (int d = 0; d < HM_DIMENSION; d++)
      u[d] = floor (urand ()*16)/16;
    point_at (x, u);
    hm_vertex_add (M, x);
  }
  /* two at the same place */
  hm_vertex_add (M, x);

  /* random points, and a curve through some */
  uint32_t first = hm_end (M, 0);
  for (int r = 0; r < 2000; r++) {
    for (int d = 0; d < HM_DIMENSION; d++)
      u[d] = urand ();
    point_at (x, u);
    hm_vertex_add (M, x);
  }
  for (uint32_t v = first; v + 1 < first + 200; v++)
    hm_add (M, 1, (uint32_t []) { v, v + 1 });

  x[0] = X0[0] + H0 + 1e-9;
  CHECK (hm_vertex_move (M, mid, x) == -1, "moved out of the box");

  check_cache (M, "built");
  const HmTreeCache * C = hm_tree_cache (M);
  CHECK (C && C->duplicates >= 1, "the duplicate not found");
  CHECK (hm_tree_cache (M) == C, "the same cache");

  /* the points moved a little, staying in the box */
  for (uint32_t v = first; v < hm_end (M, 0); v++) {
    for (int d = 0; d < HM_DIMENSION; d++) {
      x[d] = hm_x (M, d)[v] + 1e-3*H0*(urand () - 0.5);
      x[d] = fmin (fmax (x[d], X0[d]), X0[d] + H0);
    }
    CHECK (hm_vertex_move (M, v, x) == 0, "moved");
  }
  check_cache (M, "moved");

  /* written through hm_x (), then hm_touch () */
  hm_x (M, 0)[first] = X0[0] + 0.999*H0;
  hm_touch (M);
  check_cache (M, "touched");

  /* points removed, others added */
  for (uint32_t v = first + 300; v < first + 600; v++)
    hm_remove (M, 0, v);
  for (int r = 0; r < 100; r++) {
    for (int d = 0; d < HM_DIMENSION; d++)
      u[d] = urand ();
    point_at (x, u);
    hm_vertex_add (M, x);
  }
  check_cache (M, "removed and added");

  /* the tree refined further */
  uint32_t before = (uint32_t) tree_caches ()->nleaves;
  refine_half (depth);
  CHECK ((uint32_t) tree_caches ()->nleaves > before, "not refined");
  check_cache (M, "refined");

  /* the domain moved from under a point: no cache; back, a cache */
  tree_shift_resize (X0[0] + 0.5, X0[1], X0[2], H0);
  fprintf (stderr, "(expected) ");
  CHECK (hm_tree_cache (M) == NULL, "a point left outside the domain");
  tree_shift_resize (X0[0], X0[1], X0[2], H0);
  check_cache (M, "domain back");

  /* every level along the curve */
  int fx = hm_scalar_new (M, 0, "x"), fl = hm_scalar_new (M, 1, "length");
  hm_foreach (M, 0, v)
    hm_val (M, fx, 0)[v] = (HmReal) hm_x (M, 0)[v];
  hm_foreach (M, 1, e)
    hm_val (M, fl, 0)[e] = (HmReal) hm_measure (M, 1, e);
  uint32_t * remap[HM_LEVELS], end0 = hm_end (M, 0);
  double * oldx = malloc (end0 * sizeof (double));
  for (uint32_t v = 0; v < end0; v++)
    oldx[v] = hm_x (M, 0)[v];
  CHECK (hm_compact_morton (M, remap) == 0, "hm_compact_morton ()");
  C = hm_tree_cache (M);
  CHECK (hm_dense (M, 0) && hm_dense (M, 1), "not dense");
  uint32_t bad = 0;
  for (uint32_t a = 0; C && a < C->n; a++)
    bad += C->order[a] != a;
  CHECK (C && bad == 0, "%u points not in the order of their keys", bad);
  for (uint32_t v = 0; v < end0; v++)
    if (remap[0][v] != HM_NONE && hm_x (M, 0)[remap[0][v]] != oldx[v])
      bad++;
  CHECK (bad == 0, "remap: %u points not where it says", bad);
  hm_foreach (M, 0, v)
    bad += hm_val (M, fx, 0)[v] != (HmReal) hm_x (M, 0)[v];
  hm_foreach (M, 1, e)
    bad += fabs (hm_val (M, fl, 0)[e] - hm_measure (M, 1, e)) >
      1e-6*hm_measure (M, 1, e);
  CHECK (bad == 0, "%u field values not moved along", bad);
  uint64_t prev = 0;
  hm_foreach (M, 1, e) {
    double c[3];
    hm_centroid (M, 1, e, c);
    uint64_t k = hm_morton_key (c, M->origin, M->size);
    bad += e && k < prev;
    prev = k;
  }
  CHECK (bad == 0, "%u edges out of the curve's order", bad);
  for (uint32_t q = 0; C && q < C->nleaves; q++)
    hm_foreach_leaf_point (C, q, v)
      bad += v != a_v;
  CHECK (bad == 0, "a leaf's points not the ids of its range");
  check_cache (M, "compacted");
  free (oldx);
  for (int k = 0; k <= M->dim; k++)
    free (remap[k]);
  hmesh_free (M);
  tree_destroy ();
}

/* a mesh across the tree's domain, cut to it: its points on the domain's
   faces in the leaves there */
static void test_clip (void)
{
  int depth = HM_DIMENSION == 3 ? 4 : 6;
  tree (depth, 4);
  tree_shift_resize (X0[0], X0[1], X0[2], H0);
  refine_half (depth - 1);
  int m = HM_DIMENSION == 1 ? 1 : 2;
  Hmesh * M = hmesh_new (m);
  int n = 16;
  double u[3], x[3];
  if (m == 1) {
    for (int i = 0; i <= n; i++) {
      u[0] = -0.3 + 1.6*i/n;
      point_at (x, u);
      hm_vertex_add (M, x);
    }
    for (uint32_t i = 0; i < (uint32_t) n; i++)
      hm_add (M, 1, (uint32_t []) { i, i + 1 });
  }
  else {
    /* a square larger than the domain, rotated, at z = 0.4 */
    for (int j = 0; j <= n; j++)
      for (int i = 0; i <= n; i++) {
        double a = -0.4 + 1.8*i/n - 0.5, b = -0.4 + 1.8*j/n - 0.5;
        u[0] = 0.5 + 0.9*a - 0.44*b, u[1] = 0.5 + 0.44*a + 0.9*b, u[2] = 0.4;
        point_at (x, u);
        hm_vertex_add (M, x);
      }
    uint32_t w = (uint32_t) n + 1;
    for (uint32_t j = 0; j < (uint32_t) n; j++)
      for (uint32_t i = 0; i < (uint32_t) n; i++) {
        uint32_t a = j*w + i;
        hm_add (M, 2, (uint32_t []) { a, a + 1, a + w + 1 });
        hm_add (M, 2, (uint32_t []) { a, a + w + 1, a + w });
      }
  }
  CHECK (hm_tree_clip (M) == 0, "clip to the tree");
  uint32_t on = 0;
  hm_foreach (M, 0, v)
    on += hm_domain_faces (M, 0, v) != 0;
  CHECK (on > 0, "clip: no point on the domain's faces");
  check_cache (M, "clipped");
  double area = 0;
  hm_foreach (M, m, s)
    area += hm_measure (M, m, s);
  double want = m == 1 ? H0 : H0*H0;
  CHECK (fabs (area - want) < 1e-11*want, "clip: measure %.15g, not %.15g",
         area, want);
  hmesh_free (M);
  tree_destroy ();
}

static void test_too_deep (void)
{
#if HM_DIMENSION == 3
  tree (18, 4);                  /* 18 + 4 bits an axis, the keys 21 */
  Hmesh * M = hmesh_new (0);
  fprintf (stderr, "     (expected) ");
  CHECK (hm_tree_bind (M) == -1, "a tree too deep for the keys");
  hmesh_free (M);
  tree_destroy ();
  tree (17, 4);
  M = hmesh_new (0);
  CHECK (hm_tree_bind (M) == 0, "a tree just deep enough");
  hmesh_free (M);
  tree_destroy ();
#endif
}

/* every leaf refined to level `to` */
static void refine_all (int to)
{
  for (int pass = 0; pass < to; pass++) {
    const TreeCaches * tc = tree_caches ();
    int n = 0;
    Index * todo = malloc ((size_t) tc->nleaves * sizeof (Index));
    for (int q = 0; q < tc->nleaves; q++)
      if (tc->leaves[q].p.level < to)
        todo[n++] = tc->leaves[q].p;
    for (int a = 0; a < n; a++)
      tree_refine (todo[a]);
    free (todo);
  }
}

/* 2^20 points in about as many cells: a uniform tree of blocks of 16 */
static void bench (void)
{
  int level = HM_DIMENSION == 1 ? 16 : HM_DIMENSION == 2 ? 6 : 3;
  tree (level + 1, 4);
  refine_all (level);
  tree_shift_resize (0., 0., 0., 1.);
  Hmesh * M = hmesh_new (0);
  hm_tree_bind (M);
  uint32_t n = 1u << 20;
  for (uint32_t r = 0; r < n; r++) {
    double x[3];
    for (int d = 0; d < HM_DIMENSION; d++)
      x[d] = urand ();
    hm_vertex_add (M, x);
  }
  double t0 = now ();
  hm_tree_cache (M);
  double t1 = now ();
  for (uint32_t v = 0; v < n; v++) {
    double x[3];
    for (int d = 0; d < HM_DIMENSION; d++)
      x[d] = fmin (fmax (hm_x (M, d)[v] + 1e-5*(urand () - 0.5), 0.), 1.);
    for (int d = 0; d < HM_DIMENSION; d++)
      hm_x (M, d)[v] = x[d];
  }
  hm_touch (M);
  double t2 = now ();
  hm_tree_cache (M);
  double t3 = now ();
  hm_compact_morton (M, NULL);
  hm_touch (M);
  double t4 = now ();
  hm_tree_cache (M);
  double t5 = now ();
  /* in memory order, then moved: the points' ids still close to the order
     of their keys */
  for (uint32_t v = 0; v < n; v++)
    for (int d = 0; d < HM_DIMENSION; d++)
      hm_x (M, d)[v] = fmin (fmax (hm_x (M, d)[v] + 1e-5*(urand () - 0.5),
                                   0.), 1.);
  hm_touch (M);
  double t6 = now ();
  hm_tree_cache (M);
  double t7 = now ();
  CHECK (hm_count (M, 0) == n, "bench: %u points of %u", hm_count (M, 0), n);
  printf ("the cache of %u points, %d leaves, ns a point: %.1f from nothing, "
          "%.1f after moves of 1e-5; in memory order (hm_compact_morton ()) "
          "%.1f, %.1f after moves\n", n, tree_caches ()->nleaves,
          1e9*(t1 - t0)/n, 1e9*(t3 - t2)/n, 1e9*(t5 - t4)/n,
          1e9*(t7 - t6)/n);
  hmesh_free (M);
  tree_destroy ();
}

/* -------- zsimplex and nsimplex -------- */

/* a curve of short edges, inside the domain: in 1D along x, else a circle
   (z = 0.5 in 3D) */
static Hmesh * short_curve (int n)
{
  Hmesh * M = hmesh_new (1);
  hm_tree_bind (M);
  double u[3], x[3];
  for (int i = 0; i < n; i++) {
    double t = 2*M_PI*i/n;
    u[0] = HM_DIMENSION == 1 ? 0.05 + 0.9*i/n : 0.5 + 0.37*cos (t);
    u[1] = 0.5 + 0.37*sin (t), u[2] = 0.5;
    point_at (x, u);
    hm_vertex_add (M, x);
  }
  for (uint32_t i = 0; i + (HM_DIMENSION == 1) < (uint32_t) n; i++)
    hm_add (M, 1, (uint32_t []) { i, (i + 1) % (uint32_t) n });
  return M;
}

/* every leaf cell's run: its simplices' centroids in it, each simplex in
   one run */
static void check_zn (Hmesh * M, int k, scalar zs, scalar ns,
                      const char * when)
{
  long count = hm_tree_zn (M, k, zs, ns);
  const HmTreeSorted * S = hm_tree_sorted (M, k);
  CHECK (count == (long) hm_count (M, k) && S, "%s, level %d: %ld, not %u",
         when, k, count, hm_count (M, k));
  if (!S)
    return;
  uint8_t * seen = calloc (hm_end (M, k) + 1, 1);
  long total = 0, outside = 0, twice = 0;
  foreach () {
    double lo[3], h = Delta;
    lo[0] = x - h/2;
#if HM_DIMENSION > 1
    lo[1] = y - h/2;
#endif
#if HM_DIMENSION > 2
    lo[2] = z - h/2;
#endif
    for (int a = (int) zs[]; a < (int) zs[] + (int) ns[]; a++) {
      uint32_t s = S->order[a];
      double c[3];
      hm_centroid (M, k, s, c);
      for (int d = 0; d < HM_DIMENSION; d++)
        outside += c[d] < lo[d] - 1e-12 || c[d] > lo[d] + h + 1e-12;
      twice += seen[s]++ > 0;
    }
    total += (long) ns[];
  }
  CHECK (total == count && !outside && !twice, "%s, level %d: %ld in runs "
         "of %ld, %ld centroids outside their cell, %ld twice", when, k,
         total, count, outside, twice);
  free (seen);
}

static void test_zn (void)
{
  int depth = HM_DIMENSION == 3 ? 4 : 6;
  tree (depth, 4);
  tree_shift_resize (X0[0], X0[1], X0[2], H0);
  refine_half (depth - 1);
  scalar zs = scalar_new (), ns = scalar_new ();
  Hmesh * M = short_curve (HM_DIMENSION == 3 ? 400 : 1500);
  check_zn (M, 0, zs, ns, "built");
  check_zn (M, 1, zs, ns, "built");
  /* the points moved: the runs follow (M->stamp) */
  hm_foreach (M, 0, v) {
    double x[3];
    for (int d = 0; d < HM_DIMENSION; d++)
      x[d] = hm_x (M, d)[v] + 2e-3*H0*(urand () - 0.5);
    hm_vertex_move (M, v, x);
  }
  check_zn (M, 1, zs, ns, "moved");
  /* edges removed and added back: the runs follow */
  uint32_t ne = hm_count (M, 1);
  for (uint32_t e = 0; e < hm_end (M, 1); e += 7)
    if (hm_alive (M, 1, e))
      hm_remove (M, 1, e);
  check_zn (M, 1, zs, ns, "removed");
  hmesh_free (M);
  M = short_curve (HM_DIMENSION == 3 ? 400 : 1500);
  CHECK (hm_count (M, 1) == ne, "zn: the curve again");
  /* the tree refined: the runs of the finer cells */
  refine_half (depth);
  check_zn (M, 1, zs, ns, "refined");
  hmesh_free (M);
  scalar_free (zs), scalar_free (ns);
  tree_destroy ();
}

/* -------- edges shorter than the cells (-DEULERIAN_LAGRANGIAN) -------- */

#ifdef EULERIAN_LAGRANGIAN
static void test_el (void)
{
  int depth = HM_DIMENSION == 3 ? 4 : 6;
  tree (depth, 4);
  tree_shift_resize (X0[0], X0[1], X0[2], H0);
  refine_half (depth - 1);
  Hmesh * M = hmesh_new (HM_DIMENSION >= 2 ? 2 : 1);
  CHECK (hm_tree_bind (M) == 0, "el: bind");
  double u[3] = { 0.25, 0.5, 0.5 }, a[3], b[3], c[3];
  point_at (a, u);
  double hf = hm_tree_cell_size (M, a);
  u[0] = 0.75;
  point_at (b, u);
  double hc = hm_tree_cell_size (M, b);
  CHECK (hf > 0 && hc > hf, "el: cells %g and %g, not finer on the left",
         hf, hc);
  /* in the fine half: 0.6 of a cell fits, 1.2 does not */
  uint32_t va = hm_vertex_add (M, a);
  for (int d = 0; d < HM_DIMENSION; d++)
    b[d] = a[d], c[d] = a[d];
  b[0] += 0.6*hf, c[0] += 1.2*hf;
  uint32_t vb = hm_vertex_add (M, b), vc = hm_vertex_add (M, c);
  CHECK (hm_add (M, 1, (uint32_t []) { va, vc }) == HM_NONE,
         "el: an edge of 1.2 cells made");
  uint32_t e = hm_add (M, 1, (uint32_t []) { va, vb });
  CHECK (e != HM_NONE, "el: an edge of 0.6 cells refused");
#if HM_DIMENSION >= 2
  /* a triangle with one edge too long: refused, nothing made */
  uint32_t edges = hm_count (M, 1);
  CHECK (hm_add (M, 2, (uint32_t []) { va, vb, vc }) == HM_NONE &&
         hm_count (M, 1) == edges, "el: a triangle with a long edge made");
  double d3[3] = { a[0], a[1], a[2] };
  d3[1] += 0.5*hf;
  uint32_t vd = hm_vertex_add (M, d3);
  CHECK (hm_add (M, 2, (uint32_t []) { va, vb, vd }) != HM_NONE,
         "el: a small triangle refused");
#endif
  /* across the level jump at x = 0.5: 1.5 of the fine cell, under the
     coarse one; the smaller cell bounds it */
  u[0] = 0.5 - 1e-3;
  point_at (a, u);
  double ha = hm_tree_cell_size (M, a);
  for (int d = 0; d < HM_DIMENSION; d++)
    b[d] = a[d];
  b[0] += 1.5*ha;
  double hb = hm_tree_cell_size (M, b);
  uint32_t p = hm_vertex_add (M, a), q = hm_vertex_add (M, b);
  CHECK (hb == 2*ha && hm_add (M, 1, (uint32_t []) { p, q }) == HM_NONE,
         "el: an edge across the jump (cells %g, %g), longer than the fine "
         "cell, made", ha, hb);
  /* a move that stretches an edge: refused, the point left */
  double x0 = hm_x (M, 0)[vb];
  for (int d = 0; d < HM_DIMENSION; d++)
    c[d] = hm_x (M, d)[vb];
  c[0] += 0.6*hf;
  CHECK (hm_vertex_move (M, vb, c) == -1 && hm_x (M, 0)[vb] == x0,
         "el: a move to 1.2 cells made");
  c[0] = x0 + 0.2*hf;
  CHECK (hm_vertex_move (M, vb, c) == 0, "el: a move to 0.8 cells refused");
  CHECK (hm_long_edges (M) == 0 && hm_check (M) == 0, "el: invariants");
  /* the tree refined: the cells halve, the edges no longer fit */
  refine_half (depth);
  uint32_t longer = hm_long_edges (M);
  CHECK (longer > 0, "el: no long edge after refinement");
  fprintf (stderr, "     (expected) ");
  CHECK (hm_check (M) == -1, "el: hm_check () after refinement");
  printf ("el: %u edges too long after the refinement\n", longer);
  hmesh_free (M);

  /* clipped: every piece made, long or not, counted after */
  M = hmesh_new (1);
  u[0] = -0.2, u[1] = 0.5, u[2] = 0.5;
  point_at (a, u);
  u[0] = 1.2;
  point_at (b, u);
  p = hm_vertex_add (M, a), q = hm_vertex_add (M, b);
  hm_add (M, 1, (uint32_t []) { p, q });   /* no cell size yet: made */
  CHECK (hm_tree_clip (M) == 0 && hm_count (M, 1) == 1 &&
         hm_long_edges (M) == 1, "el: clip: %u edges, %u long",
         hm_count (M, 1), hm_long_edges (M));
  hmesh_free (M);
  tree_destroy ();
}
#endif

int main (void)
{
#ifdef EULERIAN_LAGRANGIAN
  /* the other tests make edges longer than the cells */
  (void) test_tree, (void) test_clip, (void) test_too_deep, (void) bench;
  test_el ();
  test_zn ();
#else
  test_tree ();
  test_clip ();
  test_too_deep ();
  test_zn ();
  bench ();
#endif
  if (fails)
    printf ("%d failed\n", fails);
  else
    printf ("passed, %dD\n", HM_DIMENSION);
  return fails != 0;
}
