/* The level indices: memindex, voxel-memindex and nano-vdb, side by side
 * against a reference hash map, through random assigns, overwrites and
 * frees, clustered and scattered; lookups by get, by cursor around every
 * entry, out of range; emptying, destroying and reusing. Build with
 * -DNVDB_COMPACT for nano-vdb's compact nodes. */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>

#include "memindex.h"
#include "voxel-memindex.h"
#include "nano-vdb.h"

static long checks, fails;
#define CHECK(c) do { checks++; if (!(c)) { if (fails++ < 20) \
  printf ("FAIL line %d: %s\n", __LINE__, #c); } } while (0)

#define D dimension
typedef struct { int x [D]; } Pt;
#if D == 1
# define P(p)  (p).x [0]
#elif D == 2
# define P(p)  (p).x [0], (p).x [1]
#else
# define P(p)  (p).x [0], (p).x [1], (p).x [2]
#endif

/* ---- the reference: open addressing on the packed indices ------------- */

typedef struct { uint64_t key; void * v; } Slot;
static Slot * ref;
static size_t refcap, refn;

static uint64_t key_of (Pt p)
{
  uint64_t k = 1;                       /* never 0: 0 marks an empty slot */
  for (int d = 0; d < D; d++)
    k = k << 21 | (uint64_t) p.x [d];
  return k;
}

static size_t ref_find (uint64_t k)
{
  size_t h = (k * 0x9E3779B97F4A7C15ull) >> 40;
  for (h &= refcap - 1; ref [h].key && ref [h].key != k; h = (h + 1) & (refcap - 1))
    ;
  return h;
}

static void * ref_get (Pt p)
{
  size_t h = ref_find (key_of (p));
  return ref [h].key ? ref [h].v : NULL;
}

static void ref_set (Pt p, void * v)     /* v NULL: remove (a tombstone) */
{
  size_t h = ref_find (key_of (p));
  if (!ref [h].key) {
    if (!v)
      return;
    ref [h].key = key_of (p);
    refn++;
  }
  ref [h].v = v;
}

/* ---- the three indices ------------------------------------------------ */

static struct _Memindex mi;
static struct _VMemindex vi;
static struct _NanoVDB ni;

#define CALL(f, ...)  f (__VA_ARGS__)

static void assign (Pt p, int len, void * v)
{
  CALL (mem_assign, &mi, P (p), len, v);
  CALL (vmem_assign, &vi, P (p), len, v);
  CALL (nvdb_assign, &ni, P (p), len, v);
  ref_set (p, v);
}

static void release (Pt p, int len)
{
  CALL (mem_free, &mi, P (p), len, NULL);
  CALL (vmem_free, &vi, P (p), len, NULL);
  CALL (nvdb_free, &ni, P (p), len, NULL);
  ref_set (p, NULL);
}

/* every index answers as the reference, by get and by allocated */
static void agree (Pt p)
{
  void * r = ref_get (p);
  CHECK (CALL (mem_get, &mi, P (p)) == r);
  CHECK (CALL (vmem_get, &vi, P (p)) == r);
  CHECK (CALL (nvdb_get, &ni, P (p)) == r);
  CHECK (CALL (nvdb_allocated, &ni, P (p)) == (r != NULL));
  CHECK (CALL (vmem_allocated, &vi, P (p)) == (r != NULL));
}

static uint64_t rng = 0x2545F4914F6CDD1Dull;
static int rnd (int n)
{
  rng ^= rng << 13; rng ^= rng >> 7; rng ^= rng << 17;
  return (int) (rng % (uint64_t) n);
}

static Pt random_pt (int len, Pt centre, int clustered)
{
  Pt p;
  for (int d = 0; d < D; d++) {
    int w = len < 24 ? len : 24;
    p.x [d] = clustered ? centre.x [d] - w / 2 + rnd (w) : rnd (len);
    if (p.x [d] < 0) p.x [d] = 0;
    if (p.x [d] >= len) p.x [d] = len - 1;
  }
  return p;
}

/* an index far outside the level: past the end, or negative */
static int q_far (int d, int len) { return d ? -1 - len : len + 5; }

static void run_level (int level, int steps)
{
  int len = 1 << level;
  refcap = 1 << 16;
  ref = calloc (refcap, sizeof (Slot));
  refn = 0;
  Pt * live = malloc ((size_t) steps * sizeof (Pt));
  int nlive = 0;
  char * tag = malloc ((size_t) steps + 1);
  Pt centre;
  for (int d = 0; d < D; d++)
    centre.x [d] = rnd (len);

  for (int s = 0; s < steps; s++) {
    int op = rnd (10);
    if (op < 6 || !nlive) {                     /* assign, or overwrite */
      Pt p = random_pt (len, centre, op < 4);
      if (!ref_get (p))
        live [nlive++] = p;
      assign (p, len, tag + s + 1);             /* any non-NULL address */
    }
    else if (op < 8) {                          /* overwrite a live one */
      Pt p = live [rnd (nlive)];
      assign (p, len, tag + s + 1);
    }
    else {                                      /* free a live one */
      int q = rnd (nlive);
      release (live [q], len);
      live [q] = live [--nlive];
    }
    if (rnd (50) == 0)
      for (int d = 0; d < D; d++)
        centre.x [d] = rnd (len);
  }

  /* every live entry, the indices around it, by get and by cursor */
  MemCursor mc = {0};
  VMemCursor vc = {0};
  NvdbCursor nc = {0};
  for (int q = 0; q < nlive; q++) {
    agree (live [q]);
    for (int d = 0; d < D; d++)
      for (int o = -2; o <= 2; o++) {
        Pt p = live [q];
        p.x [d] += o;                           /* out of range too */
        agree (p);
        void * r = ref_get (p);
        CHECK (CALL (mem_get_cursor, &mi, &mc, P (p)) == r);
        CHECK (CALL (vmem_get_cursor, &vi, &vc, P (p)) == r);
        CHECK (CALL (nvdb_get_cursor, &ni, &nc, P (p)) == r);
      }
  }
  /* random indices over the level, and far outside it */
  for (int q = 0; q < 2000; q++) {
    Pt p = random_pt (len, centre, q & 1);
    agree (p);
    CHECK (CALL (nvdb_get_cursor, &ni, &nc, P (p)) == ref_get (p));
  }
  Pt far;
  for (int d = 0; d < D; d++)
    far.x [d] = q_far (d, len);
  agree (far);

  /* empty again: nano-vdb gives back every node, the upper ones too */
  while (nlive)
    release (live [--nlive], len);
  for (int q = 0; q < 500; q++)
    agree (random_pt (len, centre, q & 1));
  CHECK (ni.single == NULL && ni.root_n == 0);

  /* destroyed, and reused at another level */
  mem_destroy (&mi, len);
  mi = (struct _Memindex) {0};
  vmem_destroy (&vi, len);
  nvdb_destroy (&ni, len);
  CHECK (ni.len == 0 && vi.len == 0);
  for (int q = 0; q < 50; q++)
    CHECK (CALL (nvdb_get, &ni, P (random_pt (len, centre, 0))) == NULL);
  free (ref);
  free (live);
  free (tag);
}


int main (void)
{
  /* empty indices: destroy is harmless, every lookup NULL */
  mem_destroy (&mi, 0);
  vmem_destroy (&vi, 0);
  nvdb_destroy (&ni, 0);
  Pt o = {{0}};
  CHECK (CALL (mem_get, &mi, P (o)) == NULL);
  CHECK (CALL (vmem_get, &vi, P (o)) == NULL);
  CHECK (CALL (nvdb_get, &ni, P (o)) == NULL);

  /* small and large levels; above 12, nano-vdb's root is a hash table */
#if D == 1
  int levels [] = { 0, 1, 3, 6, 7, 12, 13, 16, 20 };
#elif D == 2
  int levels [] = { 0, 1, 3, 5, 7, 12, 13, 15 };
#else
  int levels [] = { 0, 1, 2, 4, 7, 12, 13 };
#endif
  for (size_t q = 0; q < sizeof levels / sizeof levels [0]; q++)
    run_level (levels [q], 20000);

  printf ("dimension %d%s: %ld checks, %ld failed\n", dimension,
#ifdef NVDB_COMPACT
          " (nano-vdb compact)",
#else
          "",
#endif
          checks, fails);
  return fails != 0;
}
