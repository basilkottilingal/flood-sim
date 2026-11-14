/*
.. The manifold in libb-amr's tree: see hmesh-tree.h. Built with cwrap, as
.. tree.h is.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "tree.h"
#include "block-ghosts.h"
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

static void cache_free (void * p)
{
  HmTreeCache * C = p;
  if (!C)
    return;
  free (C->key), free (C->order), free (C->leaf_of);
  free (C->start), free (C->stop);
  for (int k = 0; k < HM_LEVELS; k++)
    free (C->sorted[k].key), free (C->sorted[k].order);
  free (C->table_key), free (C->table_level);
  free (C);
}

static HmTreeCache * cache_of (Hmesh * M)
{
  if (!M->ext) {
    M->ext = calloc (1, sizeof (HmTreeCache));
    if (!M->ext) {
      fprintf (stderr, "hmesh: out of memory\n");
      abort ();
    }
    M->ext_free = cache_free;
  }
  return M->ext;
}

#ifdef EULERIAN_LAGRANGIAN
static double tree_cell_size (const double * x, void * ctx);
#endif

/* the bits an axis the tree's finest cells need, if at most the keys', or
   0 */
static int depth_ok (void)
{
  int m = tree_get_block_m (), bits = 0;
  while ((1 << bits) < m)
    bits++;
  int need = tree_get_depth_max () + bits;
  if (need > HM_MORTON_BITS) {
    fprintf (stderr, "hmesh: the tree's cells need %d bits an axis (depth "
             "%d, blocks of %d cells), more than the %d of the keys\n",
             need, tree_get_depth_max (), m, HM_MORTON_BITS);
    return 0;
  }
  return need;
}

static int compare_keys (const void * a, const void * b)
{
  uint64_t x = *(const uint64_t *) a, y = *(const uint64_t *) b;
  return x < y ? -1 : x > y;
}

/* the points with the key of another, within each cell's run of the order
   (sorted to the cells, not within them): pairs in a short run, sorted
   keys in a long one */
static uint32_t duplicates (const HmTreeCache * C, int lo)
{
  uint32_t dup = 0, n = C->n;
  uint64_t * tmp = NULL;
  for (uint32_t a = 0, b; a < n; a = b) {
    b = a + 1;
    while (b < n && C->key[b] >> lo == C->key[a] >> lo)
      b++;
    if (b - a <= 16) {
      for (uint32_t i = a + 1; i < b; i++)
        for (uint32_t j = a; j < i; j++)
          if (C->key[i] == C->key[j]) {
            dup++;
            break;
          }
    }
    else {
      tmp = grow (tmp, b - a, sizeof (uint64_t));
      memcpy (tmp, C->key + a, (b - a) * sizeof (uint64_t));
      qsort (tmp, b - a, sizeof (uint64_t), compare_keys);
      for (uint32_t i = 1; i < b - a; i++)
        dup += tmp[i] == tmp[i - 1];
    }
  }
  free (tmp);
  return dup;
}

static void tree_box (double * origin, double * size)
{
  for (int d = 0; d < HM_DIMENSION; d++)
    origin[d] = tree_get_origin (d);
  *size = tree_get_size ();
}

int hm_tree_bind (Hmesh * M)
{
  if (!depth_ok ())
    return -1;
  double o[HM_DIMENSION], h;
  tree_box (o, &h);
  if (hm_set_domain (M, o, h)) {
    fprintf (stderr, "hmesh: a point is outside the tree's domain\n");
    return -1;
  }
#ifdef EULERIAN_LAGRANGIAN
  hm_set_cell_size (M, tree_cell_size, M);
#endif
  return 0;
}

int hm_tree_clip (Hmesh * M)
{
  if (!depth_ok ())
    return -1;
  double o[HM_DIMENSION], h;
  tree_box (o, &h);
  if (hm_clip (M, o, h))
    return -1;
#ifdef EULERIAN_LAGRANGIAN
  hm_set_cell_size (M, tree_cell_size, M);
#endif
  return 0;
}

static int same_box (const Hmesh * M)
{
  double o[HM_DIMENSION], h;
  tree_box (o, &h);
  if (!M->domain || M->size != h)
    return 0;
  for (int d = 0; d < HM_DIMENSION; d++)
    if (M->origin[d] != o[d])
      return 0;
  return 1;
}

/* the points sorted again, from the last order: those still in use in it,
   then the others */
static void sort_points (Hmesh * M, HmTreeCache * C)
{
  uint32_t n = hm_count (M, 0), end = hm_end (M, 0);
  if (n > C->cap || C->cap == 0) {
    C->cap = n > 2 * C->cap ? n : 2 * C->cap;
    C->key = grow (C->key, C->cap, sizeof (uint64_t));
    C->order = grow (C->order, C->cap, sizeof (uint32_t));
  }
  if (end > C->leaf_cap || C->leaf_cap == 0) {
    C->leaf_cap = end > 2 * C->leaf_cap ? end : 2 * C->leaf_cap;
    C->leaf_of = grow (C->leaf_of, C->leaf_cap, sizeof (uint32_t));
  }
  uint32_t * mark = C->leaf_of, a = 0;
  memset (mark, 0, (size_t) end * sizeof (uint32_t));
  for (uint32_t b = 0; C->valid && b < C->n; b++) {
    uint32_t v = C->order[b];
    if (v < end && hm_alive (M, 0, v) && !mark[v]) {
      mark[v] = 1;
      C->order[a++] = v;
    }
  }
  hm_foreach (M, 0, v)
    if (!mark[v])
      C->order[a++] = v;
  /* the keys in the order of the ids, reading the coordinates in memory
     order, then one key gathered a point (not D coordinates) */
  uint64_t * byid = grow (NULL, end, sizeof (uint64_t));
  hm_foreach (M, 0, v) {
    double x[HM_DIMENSION];
    for (int d = 0; d < HM_DIMENSION; d++)
      x[d] = M->x[d][v];
    byid[v] = hm_morton_key (x, M->origin, M->size);
  }
  for (a = 0; a < n; a++)
    C->key[a] = byid[C->order[a]];
  free (byid);
  /* to the finest cells: the bits below them left out */
  int lo = HM_DIMENSION * (HM_MORTON_BITS - C->cell_bits);
  hm_sort_keys (C->key, C->order, n, lo, HM_DIMENSION * HM_MORTON_BITS);
  C->n = n;
  C->duplicates = duplicates (C, lo);
}

/* the leaves' ranges: the leaves and the keys are both along the curve */
static void leaves_walk (Hmesh * M, HmTreeCache * C)
{
  const TreeCaches * tc = tree_caches ();
  uint32_t nl = (uint32_t) tc->nleaves, n = C->n, a = 0;
  if (nl > C->nleaves_cap || C->nleaves_cap == 0) {
    C->nleaves_cap = nl > 2 * C->nleaves_cap ? nl : 2 * C->nleaves_cap;
    C->start = grow (C->start, C->nleaves_cap, sizeof (uint32_t));
    C->stop = grow (C->stop, C->nleaves_cap, sizeof (uint32_t));
  }
  for (uint32_t v = 0; v < hm_end (M, 0); v++)
    C->leaf_of[v] = HM_NONE;
  C->outside = 0;
  for (uint32_t q = 0; q < nl; q++) {
    Index p = tc->leaves[q].p;
    int ijk[3] = { p.i, 0, 0 };
#if HM_DIMENSION > 1
    ijk[1] = p.j;
#endif
#if HM_DIMENSION > 2
    ijk[2] = p.k;
#endif
    int shift = HM_DIMENSION * (HM_MORTON_BITS - p.level);
    uint64_t lo = hm_morton_cell (ijk, p.level);
    while (a < n && C->key[a] < lo)        /* in no leaf: a pruned part */
      a++, C->outside++;
    C->start[q] = a;
    while (a < n && (shift >= 64 || C->key[a] >> shift == lo >> shift))
      C->leaf_of[C->order[a++]] = q;
    C->stop[q] = a;
  }
  C->outside += n - a;
  C->nleaves = nl;
}

const HmTreeCache * hm_tree_cache (Hmesh * M)
{
  HmTreeCache * C = cache_of (M);
  int bits = depth_ok ();
  if (!bits)
    return NULL;
  if (!same_box (M) && hm_tree_bind (M))
    return NULL;
  int points = !C->valid || C->mesh_generation != M->generation ||
    C->cell_bits != bits;
  C->cell_bits = bits;
  if (!points && C->tree_generation == tree_generation)
    return C;
  if (points)
    sort_points (M, C);
  leaves_walk (M, C);
  C->valid = 1;
  C->mesh_generation = M->generation;
  C->tree_generation = tree_generation;
  return C;
}

int hm_compact_morton (Hmesh * M, uint32_t ** remap)
{
  if (!hm_tree_cache (M))
    return -1;
  HmTreeCache * C = M->ext;
  /* the points by their whole keys, within the cells too: each cell's run
     stays where it is, as the cache's order is sorted to the cells */
  uint32_t n0 = C->n;
  uint64_t * key0 = grow (NULL, n0, sizeof (uint64_t));
  uint32_t * order[HM_LEVELS] = { grow (NULL, n0, sizeof (uint32_t)),
                                  NULL, NULL, NULL };
  memcpy (key0, C->key, (size_t) n0 * sizeof (uint64_t));
  memcpy (order[0], C->order, (size_t) n0 * sizeof (uint32_t));
  hm_sort_keys (key0, order[0], n0, 0, HM_DIMENSION * HM_MORTON_BITS);
  uint64_t * key = NULL;
  for (int k = 1; k <= M->dim; k++) {
    uint32_t n = hm_count (M, k), a = 0;
    order[k] = grow (NULL, n, sizeof (uint32_t));
    key = grow (key, n, sizeof (uint64_t));
    hm_foreach (M, k, s) {
      double c[HM_DIMENSION];
      hm_centroid (M, k, s, c);
      key[a] = hm_morton_key (c, M->origin, M->size);
      order[k][a++] = s;
    }
    hm_sort_keys (key, order[k], n, 0, HM_DIMENSION * HM_MORTON_BITS);
  }
  free (key);
  int r = hm_compact_order (M, order, remap);
  for (int k = 0; k <= M->dim; k++)
    free (order[k]);
  if (r) {
    free (key0);
    return -1;
  }

  /* the points are in the order of their keys now: point a is order[a],
     the leaves' ranges the same */
  free (C->key);
  C->key = key0;
  C->cap = n0;
  for (uint32_t a = 0; a < C->n; a++)
    C->order[a] = a;
  for (uint32_t v = 0; v < C->n; v++)
    C->leaf_of[v] = HM_NONE;
  for (uint32_t q = 0; q < C->nleaves; q++)
    for (uint32_t a = C->start[q]; a < C->stop[q]; a++)
      C->leaf_of[a] = q;
  C->mesh_generation = M->generation;
  return 0;
}

/* -------- points in the tree -------- */

/* the leaves' keys and levels, made again when the tree changed */
static void table_update (HmTreeCache * C)
{
  if (C->table_valid && C->table_generation == tree_generation)
    return;
  const TreeCaches * tc = tree_caches ();
  uint32_t nl = (uint32_t) tc->nleaves;
  if (nl > C->table_cap || !C->table_cap) {
    C->table_cap = nl > 2 * C->table_cap ? nl : 2 * C->table_cap;
    C->table_key = grow (C->table_key, C->table_cap, sizeof (uint64_t));
    C->table_level = grow (C->table_level, C->table_cap, sizeof (int));
  }
  for (uint32_t q = 0; q < nl; q++) {
    Index p = tc->leaves[q].p;
    int ijk[3] = { p.i, 0, 0 };
#if D > 1
    ijk[1] = p.j;
#endif
#if D > 2
    ijk[2] = p.k;
#endif
    C->table_key[q] = hm_morton_cell (ijk, p.level);
    C->table_level[q] = p.level;
  }
  C->table_n = nl;
  C->table_generation = tree_generation;
  C->table_valid = 1;
}

uint32_t hm_tree_locate (Hmesh * M, const double * x)
{
  double o[D], h;
  tree_box (o, &h);
  for (int d = 0; d < D; d++)
    if (!(x[d] >= o[d] && x[d] <= o[d] + h))
      return HM_NONE;
  HmTreeCache * C = cache_of (M);
  table_update (C);
  if (!C->table_n)
    return HM_NONE;
  uint64_t k = hm_morton_key (x, o, h);
  uint32_t a = 0, b = C->table_n;        /* the last leaf whose key <= k */
  while (b - a > 1) {
    uint32_t c = (a + b) / 2;
    if (C->table_key[c] <= k)
      a = c;
    else
      b = c;
  }
  int shift = D * (HM_MORTON_BITS - C->table_level[a]);
  if (C->table_key[a] > k ||
      (shift < 64 && k >> shift != C->table_key[a] >> shift))
    return HM_NONE;                      /* a pruned part */
  return a;
}

double hm_tree_cell_size (Hmesh * M, const double * x)
{
  uint32_t q = hm_tree_locate (M, x);
  if (q == HM_NONE)
    return 0.;
  HmTreeCache * C = M->ext;
  return tree_get_size () / ((double) tree_get_block_m () *
                             (double) (1L << C->table_level[q]));
}

#ifdef EULERIAN_LAGRANGIAN
static double tree_cell_size (const double * x, void * ctx)
{
  return hm_tree_cell_size (ctx, x);
}
#endif

/* -------- the simplices of a level by cell -------- */

/* the k-simplices sorted by their centroids' keys, to the finest cells,
   from the last order: those still in use in it, then the others */
static void sort_level (Hmesh * M, HmTreeSorted * S, int k, int bits)
{
  uint32_t n = hm_count (M, k), end = hm_end (M, k), a = 0;
  if (n > S->cap || !S->cap) {
    S->cap = n > 2 * S->cap ? n : 2 * S->cap;
    S->key = grow (S->key, S->cap, sizeof (uint64_t));
    S->order = grow (S->order, S->cap, sizeof (uint32_t));
  }
  uint8_t * mark = calloc (end ? end : 1, 1);
  uint64_t * byid = grow (NULL, end, sizeof (uint64_t));
  for (uint32_t b = 0; S->valid && b < S->n; b++) {
    uint32_t s = S->order[b];
    if (s < end && hm_alive (M, k, s) && !mark[s]) {
      mark[s] = 1;
      S->order[a++] = s;
    }
  }
  hm_foreach (M, k, s) {
    if (!mark[s])
      S->order[a++] = s;
    double c[D];
    hm_centroid (M, k, s, c);
    byid[s] = hm_morton_key (c, M->origin, M->size);
  }
  for (a = 0; a < n; a++)
    S->key[a] = byid[S->order[a]];
  free (byid);
  free (mark);
  hm_sort_keys (S->key, S->order, n, D * (HM_MORTON_BITS - bits),
                D * HM_MORTON_BITS);
  S->n = n;
  S->cell_bits = bits;
  S->stamp = M->stamp;
  S->valid = 1;
}

const HmTreeSorted * hm_tree_sorted (Hmesh * M, int k)
{
  if (k < 0 || k > M->dim)
    return NULL;
  HmTreeCache * C = cache_of (M);
  int bits = depth_ok ();
  if (!bits)
    return NULL;
  if (!same_box (M) && hm_tree_bind (M))
    return NULL;
  HmTreeSorted * S = &C->sorted[k];
  if (!S->valid || S->stamp != M->stamp || S->cell_bits != bits)
    sort_level (M, S, k, bits);
  return S;
}

/* the run of the cell centred at p, of side delta, in S */
static const HmTreeSorted * zn_sorted;
static const Hmesh * zn_mesh;

static void zn_run (double x, double y, double z, double delta, Real * first,
                    Real * count)
{
  const HmTreeSorted * S = zn_sorted;
  const Hmesh * M = zn_mesh;
  double p[3] = { x, y, z };
  uint64_t k = hm_morton_key (p, M->origin, M->size);
  int cbits = (int) lround (log2 (M->size / delta));
  int shift = D * (HM_MORTON_BITS - cbits), ss = D * (HM_MORTON_BITS -
                                                      S->cell_bits);
  /* the keys of the cell, at the sort's resolution: [lo, hi) */
  uint64_t lo = shift >= 64 ? 0 : (k >> shift) << shift >> ss;
  uint64_t hi = shift >= 64 ? UINT64_MAX : lo + ((uint64_t) 1 << (shift - ss));
  uint32_t b[2];
  for (int e = 0; e < 2; e++) {
    uint64_t t = e ? hi : lo;
    uint32_t l = 0, r = S->n;          /* the first key >> ss >= t */
    while (l < r) {
      uint32_t c = (l + r) / 2;
      if (S->key[c] >> ss < t)
        l = c + 1;
      else
        r = c;
    }
    b[e] = l;
  }
  if (shift >= 64)
    b[1] = S->n;
  *first = (Real) b[0];
  *count = (Real) (b[1] - b[0]);
}

long hm_tree_zn (Hmesh * M, int k, scalar zsimplex, scalar nsimplex)
{
  const HmTreeSorted * S = hm_tree_sorted (M, k);
  if (!S)
    return -1;
  if (sizeof (Real) < sizeof (double) && S->n > (1u << 24))
    return -1;                         /* not exact in a float */
  zn_sorted = S, zn_mesh = M;
  foreach () {
    Real first, count;
#if D == 1
    zn_run (x, 0., 0., Delta, &first, &count);
#elif D == 2
    zn_run (x, y, 0., Delta, &first, &count);
#else
    zn_run (x, y, z, Delta, &first, &count);
#endif
    zsimplex[] = first;
    nsimplex[] = count;
  }
  zn_sorted = NULL, zn_mesh = NULL;
  return (long) S->n;
}
