/*
.. nano-vdb: a sparse array of pointers as a VDB tree (root hash, upper,
.. lower, leaf). See nano-vdb.h.
*/

#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "nano-vdb.h"

NanoVDB nvdb_new (int len)
{
  (void) len;
  return calloc (1, sizeof (struct _NanoVDB));
}

/* ---- nodes ---------------------------------------------------------- */

static size_t node_bytes (int l)
{
  size_t w = NVDB_WORDS (l);
#ifdef NVDB_COMPACT
  return sizeof (NvNode) + w * sizeof (uint64_t)
    + ((w * sizeof (uint32_t) + 7) & ~(size_t) 7);
#else
  return sizeof (NvNode) + w * sizeof (uint64_t)
    + (size_t) NVDB_SLOTS (l) * sizeof (void *);
#endif
}

static NvNode * node_new (NanoVDB m, int l)
{
  NvNode * n = pool_alloc (m->pool [l]);          /* zeroed */
#ifndef NVDB_COMPACT
  n->c = (void **) (n->mask + NVDB_WORDS (l));
#endif
  return n;
}

static void node_release (NanoVDB m, int l, NvNode * n)
{
#ifdef NVDB_COMPACT
  free (n->c);
#endif
  pool_free (m->pool [l], n);
}

/* slot s of n (level l) holds p, not NULL */
static void node_put (NvNode * n, int l, int s, void * p)
{
  int w = s >> 6;
  uint64_t bit = 1ull << (s & 63);
#ifdef NVDB_COMPACT
  uint32_t * prefix = (uint32_t *) (n->mask + NVDB_WORDS (l));
  uint32_t r = prefix [w] + __builtin_popcountll (n->mask [w] & (bit - 1));
  if (n->mask [w] & bit) {
    n->c [r] = p;
    return;
  }
  if (n->n == n->cap) {
    n->cap = n->cap ? 2 * n->cap : 4;
    if (n->cap > (uint32_t) NVDB_SLOTS (l))
      n->cap = NVDB_SLOTS (l);
    n->c = realloc (n->c, n->cap * sizeof (void *));
    if (!n->c)
      abort ();
  }
  memmove (n->c + r + 1, n->c + r, (n->n - r) * sizeof (void *));
  n->c [r] = p;
  for (int v = w + 1; v < NVDB_WORDS (l); v++)
    prefix [v]++;
#else
  (void) l;
  n->c [s] = p;
  if (n->mask [w] & bit)
    return;
#endif
  n->mask [w] |= bit;
  n->n++;
}

/* slot s of n (level l) empty */
static void node_clear (NvNode * n, int l, int s)
{
  int w = s >> 6;
  uint64_t bit = 1ull << (s & 63);
  if (!(n->mask [w] & bit))
    return;
#ifdef NVDB_COMPACT
  uint32_t * prefix = (uint32_t *) (n->mask + NVDB_WORDS (l));
  uint32_t r = prefix [w] + __builtin_popcountll (n->mask [w] & (bit - 1));
  memmove (n->c + r, n->c + r + 1, (n->n - r - 1) * sizeof (void *));
  for (int v = w + 1; v < NVDB_WORDS (l); v++)
    prefix [v]--;
#else
  (void) l;
  n->c [s] = NULL;
#endif
  n->mask [w] &= ~bit;
  n->n--;
}

/* ---- the root: upper nodes by key ----------------------------------- */

static void root_insert (NanoVDB m, uint64_t key, NvNode * u);

static void root_grow (NanoVDB m)
{
  NvTile * old = m->root;
  uint32_t cap = m->root_cap;
  m->root_cap = cap ? 2 * cap : 16;
  m->root = calloc (m->root_cap, sizeof (NvTile));
  if (!m->root)
    abort ();
  m->root_n = 0;
  for (uint32_t h = 0; h < cap; h++)
    if (old [h].node)
      root_insert (m, old [h].key, old [h].node);
  free (old);
}

static void root_insert (NanoVDB m, uint64_t key, NvNode * u)
{
  if (2 * (m->root_n + 1) > m->root_cap)
    root_grow (m);
  uint32_t h = nvdb_hash (key, m->root_cap);
  while (m->root [h].node)
    h = (h + 1) & (m->root_cap - 1);
  m->root [h] = (NvTile) { key, u };
  m->root_n++;
}

/* linear probing: the entries after a deleted one shift back */
static void root_delete (NanoVDB m, uint64_t key)
{
  uint32_t mask = m->root_cap - 1, i = nvdb_hash (key, m->root_cap);
  while (m->root [i].key != key || !m->root [i].node)
    i = (i + 1) & mask;
  for (uint32_t j = i; ; ) {
    j = (j + 1) & mask;
    if (!m->root [j].node)
      break;
    uint32_t k = nvdb_hash (m->root [j].key, m->root_cap);
    /* j stays if its home k lies cyclically in (i, j] */
    if (i < j ? (i < k && k <= j) : (i < k || k <= j))
      continue;
    m->root [i] = m->root [j];
    i = j;
  }
  m->root [i].node = NULL;
  m->root_n--;
}

/* ---- the interface -------------------------------------------------- */

static void start (NanoVDB m, int len)
{
  assert (len > 0 && !(len & (len - 1)));
  m->len = len;
  for (int l = 0; l < 3; l++)
    m->pool [l] = pool_new (node_bytes (l));
}

void nvdb_assign (NanoVDB m, NVDB_ARGS, int len, void * b)
{
  assert (b != NULL);
  if (!m->len)
    start (m, len);
  assert (NVDB_IN (m));
  int x [dimension] = { NVDB_IDX };
  NvNode * u;
  if (m->len <= 1 << NVDB_ROOT_SPAN) {
    if (!m->single)
      m->single = node_new (m, 2);
    u = m->single;
  }
  else {
    uint64_t key = nvdb_key (x);
    if (!(u = nvdb_root_find (m, key))) {
      u = node_new (m, 2);
      root_insert (m, key, u);
    }
  }
  int s2 = nvdb_slot (x, 2), s1 = nvdb_slot (x, 1);
  NvNode * l = nvdb_child (u, 2, s2);
  if (!l) {
    l = node_new (m, 1);
    node_put (u, 2, s2, l);
  }
  NvNode * f = nvdb_child (l, 1, s1);
  if (!f) {
    f = node_new (m, 0);
    node_put (l, 1, s1, f);
  }
  node_put (f, 0, nvdb_slot (x, 0), b);
}

void nvdb_free (NanoVDB m, NVDB_ARGS, int len, void * b)
{
  (void) len;
  (void) b;
  if (!NVDB_IN (m))
    return;
  int x [dimension] = { NVDB_IDX };
  NvNode * u = nvdb_upper (m, x);
  if (!u)
    return;
  int s2 = nvdb_slot (x, 2), s1 = nvdb_slot (x, 1);
  NvNode * l = nvdb_child (u, 2, s2);
  if (!l)
    return;
  NvNode * f = nvdb_child (l, 1, s1);
  if (!f)
    return;
  node_clear (f, 0, nvdb_slot (x, 0));
  /* empty nodes go back to their pools, from the leaf up */
  if (f->n)
    return;
  node_release (m, 0, f);
  node_clear (l, 1, s1);
  if (l->n)
    return;
  node_release (m, 1, l);
  node_clear (u, 2, s2);
  if (u->n)
    return;
  node_release (m, 2, u);
  if (m->len <= 1 << NVDB_ROOT_SPAN)
    m->single = NULL;
  else
    root_delete (m, nvdb_key (x));
}

#ifdef NVDB_COMPACT
/* the children arrays of an upper node and of every node below it */
static void free_arrays (NvNode * u)
{
  for (uint32_t a = 0; a < u->n; a++) {
    NvNode * l = u->c [a];
    for (uint32_t b = 0; b < l->n; b++)
      free (((NvNode *) l->c [b])->c);
    free (l->c);
  }
  free (u->c);
}
#endif

/* every node goes with the pools; the struct itself stays, as in memindex */
void nvdb_destroy (NanoVDB m, int len)
{
  (void) len;
#ifdef NVDB_COMPACT
  if (m->single)
    free_arrays (m->single);
  for (uint32_t h = 0; h < m->root_cap; h++)
    if (m->root [h].node)
      free_arrays (m->root [h].node);
#endif
  for (int l = 0; l < 3; l++)
    pool_destroy (m->pool [l]);
  free (m->root);
  *m = (struct _NanoVDB) {0};
}
