/*
.. A simplicial manifold, as a half-edge mesh on every level: see hmesh.h.
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "hmesh.h"

const int HM_TAG (HM_DIMENSION) = HM_DIMENSION;

Hmesh * hmesh = NULL;

static void * xrealloc (void * p, size_t n)
{
  void * q = realloc (p, n ? n : 1);
  if (!q) {
    fprintf (stderr, "hmesh: out of memory (%zu bytes)\n", n);
    abort ();
  }
  return q;
}

static size_t words (uint32_t n)
{
  return ((size_t) n + 63) / 64;
}

/* -------- the mesh -------- */

Hmesh * hmesh_new (int dim)
{
  if (dim < 0 || dim > HM_DIMENSION)
    return NULL;
  Hmesh * M = calloc (1, sizeof (Hmesh));
  if (!M)
    return NULL;
  M->dim = dim;
  return M;
}

void hmesh_free (Hmesh * M)
{
  if (!M)
    return;
  for (int k = 0; k <= M->dim; k++) {
    HmLevel * l = &M->level[k];
    free (l->alive);
    free (l->freed);
    free (l->vert);
    free (l->twin);
    free (l->facet);
    free (l->use);
  }
  for (int c = 0; c < HM_DIMENSION; c++)
    free (M->x[c]);
  for (int f = 0; f < M->nfield; f++)
    hm_field_free (M, f);
  free (M->field);
  free (M);
}

int manifold (int dim)
{
  Hmesh * M = hmesh_new (dim);
  if (!M)
    return -1;
  hmesh_free (hmesh);
  hmesh = M;
  return 0;
}

void manifold_free (void)
{
  hmesh_free (hmesh);
  hmesh = NULL;
}

/* -------- the pools: ids of a level -------- */

static void level_grow (Hmesh * M, int k, uint32_t need)
{
  HmLevel * l = &M->level[k];
  if (need <= l->cap)
    return;
  uint64_t c = l->cap ? l->cap : 16;
  while (c < need)
    c *= 2;
  if (c > (uint64_t) HM_ID_MAX)
    c = HM_ID_MAX;
  uint32_t cap = (uint32_t) c;
  size_t ow = words (l->cap), nw = words (cap), a = (size_t) k + 1;
  l->alive = xrealloc (l->alive, nw * sizeof (uint64_t));
  memset (l->alive + ow, 0, (nw - ow) * sizeof (uint64_t));
  if (k >= 1) {
    l->vert  = xrealloc (l->vert,  a * cap * sizeof (uint32_t));
    l->twin  = xrealloc (l->twin,  a * cap * sizeof (uint32_t));
    l->facet = xrealloc (l->facet, a * cap * sizeof (uint32_t));
  }
  if (k < M->dim)
    l->use = xrealloc (l->use, (size_t) cap * sizeof (uint32_t));
  if (k == 0)
    for (int d = 0; d < HM_DIMENSION; d++)
      M->x[d] = xrealloc (M->x[d], (size_t) cap * sizeof (double));
  for (int f = 0; f < M->nfield; f++)
    if (M->field[f].level == k)
      for (int i = 0; i < M->field[f].ncomp; i++)
        M->field[f].data[i] = xrealloc (M->field[f].data[i],
                                        (size_t) cap * M->field[f].esize);
  l->cap = cap;
}

static uint32_t level_alloc (Hmesh * M, int k)
{
  HmLevel * l = &M->level[k];
  uint32_t s;
  if (l->nfreed)
    s = l->freed[--l->nfreed];
  else {
    if (l->end >= HM_ID_MAX)
      return HM_NONE;
    level_grow (M, k, l->end + 1);
    s = l->end++;
  }
  l->alive[s >> 6] |= (uint64_t) 1 << (s & 63);
  l->n++;
  for (int f = 0; f < M->nfield; f++)
    if (M->field[f].level == k) {
      size_t e = M->field[f].esize;
      for (int i = 0; i < M->field[f].ncomp; i++)
        memset ((char *) M->field[f].data[i] + (size_t) s * e, 0, e);
    }
  return s;
}

static void level_release (Hmesh * M, int k, uint32_t s)
{
  HmLevel * l = &M->level[k];
  l->alive[s >> 6] &= ~((uint64_t) 1 << (s & 63));
  l->n--;
  if (l->nfreed == l->capfreed) {
    l->capfreed = l->capfreed ? 2 * l->capfreed : 16;
    l->freed = xrealloc (l->freed, (size_t) l->capfreed * sizeof (uint32_t));
  }
  l->freed[l->nfreed++] = s;
}

/* -------- orientation -------- */

/*
.. The facet of half i of (v0 ... vk), k >= 2, in the orientation the
.. simplex induces, (-1)^i (v0 ... vi ... vk without vi): for odd i, its first
.. two vertices swapped, so that the order itself has the sign. For k = 1 the
.. facet is the point v(1-i), and the sign (-1)^i stays apart.
*/
static void facet_of (int k, const uint32_t * v, int i, uint32_t * u)
{
  int n = 0;
  for (int j = 0; j <= k; j++)
    if (j != i)
      u[n++] = v[j];
  if ((i & 1) && k >= 2) {
    uint32_t t = u[0];
    u[0] = u[1];
    u[1] = t;
  }
}

/* the parity of u as a permutation of w (n <= 3): 1, -1, or 0 if not one */
static int parity (int n, const uint32_t * w, const uint32_t * u)
{
  int p[3], inv = 0;
  for (int a = 0; a < n; a++) {
    p[a] = -1;
    for (int b = 0; b < n; b++)
      if (w[b] == u[a])
        p[a] = b;
    if (p[a] < 0)
      return 0;
  }
  for (int a = 0; a < n; a++)
    for (int b = a + 1; b < n; b++) {
      if (p[a] == p[b])
        return 0;
      inv += p[a] > p[b];
    }
  return inv & 1 ? -1 : 1;
}

/* the sign bit of half i on facet f, u the facet as facet_of () gives it */
static uint32_t half_sign (const Hmesh * M, int k, int i, uint32_t f,
                           const uint32_t * u)
{
  if (k == 1)
    return i ? HM_FLIP : 0;
  const uint32_t * w = &M->level[k - 1].vert[(size_t) k * f];
  return parity (k, w, u) < 0 ? HM_FLIP : 0;
}

/* -------- simplices -------- */

uint32_t hm_vertex_add (Hmesh * M, const double * x)
{
  uint32_t s = level_alloc (M, 0);
  if (s == HM_NONE)
    return HM_NONE;
  for (int d = 0; d < HM_DIMENSION; d++)
    M->x[d][s] = x[d];
  if (M->dim > 0)
    M->level[0].use[s] = HM_NONE;
  return s;
}

/*
.. The k-simplex with the vertices v, in any order: the k-simplices on the
.. facet v[1] ... v[k] are the cycle of its use, and the one with v[0] is the
.. one whose half there is opposite v[0].
*/
uint32_t hm_find (const Hmesh * M, int k, const uint32_t * v)
{
  if (k < 0 || k > M->dim)
    return HM_NONE;
  if (k == 0)
    return hm_alive (M, 0, v[0]) ? v[0] : HM_NONE;
  uint32_t f = hm_find (M, k - 1, v + 1);
  if (f == HM_NONE)
    return HM_NONE;
  const HmLevel * l = &M->level[k];
  uint32_t h0 = M->level[k - 1].use[f], h = h0;
  if (h0 == HM_NONE)
    return HM_NONE;
  do {
    if (l->vert[h] == v[0])
      return hm_half_simplex (k, h);
    h = l->twin[h];
  } while (h != h0);
  return HM_NONE;
}

uint32_t hm_add (Hmesh * M, int k, const uint32_t * v)
{
  if (k < 0 || k > M->dim)
    return HM_NONE;
  for (int i = 0; i <= k; i++) {
    if (!hm_alive (M, 0, v[i]))
      return HM_NONE;
    for (int j = 0; j < i; j++)
      if (v[i] == v[j])
        return HM_NONE;
  }
  if (k == 0)
    return v[0];
  uint32_t s = hm_find (M, k, v);
  if (s != HM_NONE)
    return s;

  /* the facets, found or made, before s, so that level k is not grown
     while s is written */
  uint32_t f[HM_LEVELS], sign[HM_LEVELS];
  for (int i = 0; i <= k; i++) {
    uint32_t u[HM_LEVELS - 1];
    facet_of (k, v, i, u);
    f[i] = k == 1 ? u[0] : hm_add (M, k - 1, u);
    if (f[i] == HM_NONE)
      return HM_NONE;
    sign[i] = half_sign (M, k, i, f[i], u);
  }

  s = level_alloc (M, k);
  if (s == HM_NONE)
    return HM_NONE;
  HmLevel * l = &M->level[k], * lf = &M->level[k - 1];
  for (int i = 0; i <= k; i++) {
    uint32_t h = hm_half (k, s, i), u0 = lf->use[f[i]];
    l->vert[h] = v[i];
    l->facet[h] = f[i] | sign[i];
    if (u0 == HM_NONE) {
      lf->use[f[i]] = h;
      l->twin[h] = h;
    }
    else {
      l->twin[h] = l->twin[u0];
      l->twin[u0] = h;
    }
  }
  if (k < M->dim)
    l->use[s] = HM_NONE;
  return s;
}

/* h out of the cycle on its facet: O(the cycle's length) */
static void unlink_half (Hmesh * M, int k, uint32_t h)
{
  HmLevel * l = &M->level[k], * lf = &M->level[k - 1];
  uint32_t f = l->facet[h] & ~HM_FLIP;
  if (l->twin[h] == h) {
    lf->use[f] = HM_NONE;
    return;
  }
  uint32_t p = h;
  while (l->twin[p] != h)
    p = l->twin[p];
  l->twin[p] = l->twin[h];
  if (lf->use[f] == h)
    lf->use[f] = l->twin[h];
}

int hm_remove (Hmesh * M, int k, uint32_t s)
{
  if (k < 0 || k > M->dim || !hm_alive (M, k, s))
    return -1;
  if (k < M->dim && M->level[k].use[s] != HM_NONE)
    return -1;
  for (int i = 0; k >= 1 && i <= k; i++)
    unlink_half (M, k, hm_half (k, s, i));
  level_release (M, k, s);
  return 0;
}

int hm_remove_closure (Hmesh * M, int k, uint32_t s)
{
  uint32_t f[HM_LEVELS];
  for (int i = 0; k >= 1 && i <= k && hm_alive (M, k, s); i++)
    f[i] = hm_half_facet (M, k, hm_half (k, s, i));
  if (hm_remove (M, k, s))
    return -1;
  for (int i = 0; k >= 1 && i <= k; i++)
    if (hm_alive (M, k - 1, f[i]) && M->level[k - 1].use[f[i]] == HM_NONE)
      hm_remove_closure (M, k - 1, f[i]);
  return 0;
}

/* -------- compaction -------- */

void hm_compact (Hmesh * M, uint32_t ** remap)
{
  int m = M->dim;
  uint32_t * map[HM_LEVELS] = { NULL };
  for (int k = 0; k <= m; k++) {
    HmLevel * l = &M->level[k];
    map[k] = xrealloc (NULL, (size_t) l->end * sizeof (uint32_t));
    uint32_t n = 0;
    for (uint32_t s = 0; s < l->end; s++)
      map[k][s] = hm_alive (M, k, s) ? n++ : HM_NONE;
  }

  /* each level in place: a simplex moves down (t <= s), into what was
     read already */
  for (int k = 0; k <= m; k++) {
    HmLevel * l = &M->level[k];
    uint32_t a = (uint32_t) k + 1, b = (uint32_t) k + 2;
    for (uint32_t s = 0; s < l->end; s++) {
      uint32_t t = map[k][s];
      if (t == HM_NONE)
        continue;
      for (uint32_t i = 0; k >= 1 && i < a; i++) {
        uint32_t src = a * s + i, dst = a * t + i;
        uint32_t tw = l->twin[src], fc = l->facet[src];
        l->vert[dst] = map[0][l->vert[src]];
        l->twin[dst] = a * map[k][tw / a] + tw % a;
        l->facet[dst] = map[k - 1][fc & ~HM_FLIP] | (fc & HM_FLIP);
      }
      if (k < m) {
        uint32_t u = l->use[s];
        l->use[t] = u == HM_NONE ? HM_NONE : b * map[k + 1][u / b] + u % b;
      }
      if (k == 0)
        for (int d = 0; d < HM_DIMENSION; d++)
          M->x[d][t] = M->x[d][s];
      for (int f = 0; f < M->nfield; f++)
        if (M->field[f].level == k && t != s) {
          size_t e = M->field[f].esize;
          for (int c = 0; c < M->field[f].ncomp; c++) {
            char * p = M->field[f].data[c];
            memcpy (p + (size_t) t * e, p + (size_t) s * e, e);
          }
        }
    }
    /* the bits: n ones, then zeros */
    uint32_t n = l->n;
    memset (l->alive, 0, words (l->end) * sizeof (uint64_t));
    for (uint32_t w = 0; w < n / 64; w++)
      l->alive[w] = ~(uint64_t) 0;
    if (n % 64)
      l->alive[n / 64] = ((uint64_t) 1 << (n % 64)) - 1;
    l->end = n;
    l->nfreed = 0;
  }

  for (int k = 0; k <= m; k++)
    if (remap)
      remap[k] = map[k];
    else
      free (map[k]);
}

/* -------- fields -------- */

int hm_field_new (Hmesh * M, int k, const char * name, int ncomp,
                  size_t esize)
{
  if (k < 0 || k > M->dim || ncomp < 1 || esize < 1)
    return -1;
  int f = 0;
  while (f < M->nfield && M->field[f].level >= 0)
    f++;
  if (f == M->nfield) {
    if (M->nfield == M->capfield) {
      M->capfield = M->capfield ? 2 * M->capfield : 8;
      M->field = xrealloc (M->field, (size_t) M->capfield * sizeof (HmField));
    }
    M->nfield++;
  }
  HmField * F = &M->field[f];
  F->level = k;
  F->ncomp = ncomp;
  F->esize = esize;
  F->name = name ? strdup (name) : NULL;
  F->data = xrealloc (NULL, (size_t) ncomp * sizeof (void *));
  size_t cap = M->level[k].cap;
  for (int c = 0; c < ncomp; c++) {
    F->data[c] = calloc (cap ? cap : 1, esize);
    if (!F->data[c]) {
      fprintf (stderr, "hmesh: out of memory\n");
      abort ();
    }
  }
  return f;
}

void hm_field_free (Hmesh * M, int f)
{
  if (f < 0 || f >= M->nfield || M->field[f].level < 0)
    return;
  HmField * F = &M->field[f];
  for (int c = 0; c < F->ncomp; c++)
    free (F->data[c]);
  free (F->data);
  free (F->name);
  memset (F, 0, sizeof (HmField));
  F->level = -1;
}

int hm_scalar_new (Hmesh * M, int k, const char * name)
{
  return hm_field_new (M, k, name, 1, sizeof (HmReal));
}

int hm_vector_new (Hmesh * M, int k, const char * name)
{
  return hm_field_new (M, k, name, HM_DIMENSION, sizeof (HmReal));
}

/* -------- geometry -------- */

static double det (int n, double a[3][3])
{
  switch (n) {
  case 0: return 1.;
  case 1: return a[0][0];
  case 2: return a[0][0]*a[1][1] - a[0][1]*a[1][0];
  default:
    return a[0][0]*(a[1][1]*a[2][2] - a[1][2]*a[2][1])
         - a[0][1]*(a[1][0]*a[2][2] - a[1][2]*a[2][0])
         + a[0][2]*(a[1][0]*a[2][1] - a[1][1]*a[2][0]);
  }
}

/* the edges from v0: e[j][d] = x (vj+1) - x (v0) */
static void edges (const Hmesh * M, int k, uint32_t s, double e[3][3])
{
  uint32_t v0 = hm_vertex (M, k, s, 0);
  for (int j = 0; j < k; j++) {
    uint32_t v = hm_vertex (M, k, s, j + 1);
    for (int d = 0; d < 3; d++)
      e[j][d] = d < HM_DIMENSION ? M->x[d][v] - M->x[d][v0] : 0.;
  }
}

static double factorial (int k)
{
  return k == 3 ? 6. : k == 2 ? 2. : 1.;
}

double hm_measure (const Hmesh * M, int k, uint32_t s)
{
  double e[3][3], g[3][3];
  edges (M, k, s, e);
  for (int a = 0; a < k; a++)
    for (int b = 0; b < k; b++) {
      g[a][b] = 0.;
      for (int d = 0; d < HM_DIMENSION; d++)
        g[a][b] += e[a][d]*e[b][d];
    }
  double G = det (k, g);
  return sqrt (G > 0. ? G : 0.)/factorial (k);
}

double hm_volume_signed (const Hmesh * M, uint32_t s)
{
  double e[3][3];
  int k = M->dim;
  if (k != HM_DIMENSION)
    return hm_measure (M, k, s);
  edges (M, k, s, e);
  return det (k, e)/factorial (k);
}

void hm_centroid (const Hmesh * M, int k, uint32_t s, double * c)
{
  for (int d = 0; d < HM_DIMENSION; d++) {
    c[d] = 0.;
    for (int i = 0; i <= k; i++)
      c[d] += M->x[d][hm_vertex (M, k, s, i)];
    c[d] /= k + 1;
  }
}

/* -------- invariants -------- */

#define FAIL(...) do {                                  \
    fprintf (stderr, "hm_check: " __VA_ARGS__);         \
    fputc ('\n', stderr);                               \
    free (count);                                       \
    return -1;                                          \
  } while (0)

int hm_check (const Hmesh * M)
{
  int m = M->dim;
  uint32_t * count = NULL;
  for (int k = 0; k <= m; k++) {
    const HmLevel * l = &M->level[k];
    uint32_t a = (uint32_t) k + 1;

    /* the pool */
    if (l->end > l->cap)
      FAIL ("level %d: end %u beyond cap %u", k, l->end, l->cap);
    uint32_t n = 0;
    for (size_t w = 0; w < words (l->cap); w++) {
      /* the ids below end, of this word */
      uint64_t mask = 64 * (w + 1) <= l->end ? ~(uint64_t) 0 :
        64 * w >= l->end ? 0 : ((uint64_t) 1 << (l->end - 64 * w)) - 1;
      uint64_t b = l->alive[w] & mask;
      if (b != l->alive[w])
        FAIL ("level %d: alive bits beyond end %u", k, l->end);
      n += (uint32_t) bits_popcount_u64 (b);
    }
    if (n != l->n)
      FAIL ("level %d: %u alive bits, count %u", k, n, l->n);
    if (l->nfreed + l->n != l->end)
      FAIL ("level %d: %u freed + %u in use, end %u", k, l->nfreed, l->n,
            l->end);
    for (uint32_t j = 0; j < l->nfreed; j++)
      if (l->freed[j] >= l->end || hm_alive (M, k, l->freed[j]))
        FAIL ("level %d: freed id %u in use", k, l->freed[j]);
    if (k == 0)
      continue;

    /* each simplex: its vertices, facets, signs, twins */
    const HmLevel * lf = &M->level[k - 1];
    hm_foreach (M, k, s) {
      uint32_t * v = &l->vert[a * s];
      for (uint32_t i = 0; i < a; i++) {
        if (!hm_alive (M, 0, v[i]))
          FAIL ("level %d: simplex %u, vertex %u not in use", k, s, v[i]);
        for (uint32_t j = 0; j < i; j++)
          if (v[i] == v[j])
            FAIL ("level %d: simplex %u, vertex %u twice", k, s, v[i]);
      }
      for (uint32_t i = 0; i < a; i++) {
        uint32_t h = a * s + i, f = l->facet[h] & ~HM_FLIP, u[3];
        if (!hm_alive (M, k - 1, f))
          FAIL ("level %d: half %u, facet %u not in use", k, h, f);
        facet_of (k, v, (int) i, u);
        if (k == 1 ? f != u[0] : !parity (k, &lf->vert[(size_t) k * f], u))
          FAIL ("level %d: half %u, facet %u not its vertices", k, h, f);
        if ((l->facet[h] & HM_FLIP) != half_sign (M, k, (int) i, f, u))
          FAIL ("level %d: half %u, wrong sign", k, h);
        uint32_t t = l->twin[h];
        if (t >= a * l->end || !hm_alive (M, k, t / a) ||
            (l->facet[t] & ~HM_FLIP) != f)
          FAIL ("level %d: half %u, twin %u not on facet %u", k, h, t, f);
      }
      if (hm_find (M, k, v) != s)
        FAIL ("level %d: simplex %u not found by its vertices", k, s);
    }

    /* each facet's cycle: every half on it, once */
    count = xrealloc (count, (size_t) (lf->end ? lf->end : 1) *
                      sizeof (uint32_t));
    memset (count, 0, (size_t) lf->end * sizeof (uint32_t));
    hm_foreach (M, k, s)
      for (uint32_t i = 0; i < a; i++)
        count[l->facet[a * s + i] & ~HM_FLIP]++;
    hm_foreach (M, k - 1, f) {
      uint32_t h0 = lf->use[f];
      if (!count[f]) {
        if (h0 != HM_NONE)
          FAIL ("level %d: simplex %u unused, use %u", k - 1, f, h0);
        continue;
      }
      if (h0 == HM_NONE || h0 >= a * l->end || !hm_alive (M, k, h0 / a) ||
          (l->facet[h0] & ~HM_FLIP) != f)
        FAIL ("level %d: simplex %u, use %u not a half on it", k - 1, f, h0);
      uint32_t h = h0, len = 0;
      do {
        h = l->twin[h];
        len++;
      } while (h != h0 && len <= count[f]);
      if (len != count[f])
        FAIL ("level %d: simplex %u, a cycle of %u halves, %u on it",
              k - 1, f, len, count[f]);
    }
  }
  free (count);
  return 0;
}
