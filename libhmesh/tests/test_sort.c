/*
.. The sort of Morton keys (hmesh-morton.h), against a reference: qsort by
.. bits lo to hi - 1, then by the position the key had, which is the one
.. stable order. For each case, the keys and ids hm_sort_keys () gives
.. (sorted check, insertion or radix) must be the reference's bit for bit.
.. Inputs: random 64 bit keys (bits above hi set too), all equal, sorted,
.. reversed, nearly sorted (insertion), a few keys moved far (insertion
.. gives up), few distinct fields (stability), only the lowest digit
.. varying (passes skipped); n about the chunk (256) and the digits;
.. ranges of bits: whole keys, the tree's cells, odd ranges inside. The
.. scratch is exactly hm_sort_scratch (n) bytes (ASan sees overruns) and
.. full of garbage.
..
.. HM_SORT_BENCH=1: then times the sort of 2^20 and 2^22 random keys (the
.. radix), ns a key, best of 5.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "hmesh-morton.h"

static int fails = 0, checks = 0;

#define CHECK(c, ...) do {                                        \
    checks++;                                                     \
    if (!(c)) {                                                   \
      fprintf (stderr, "%s:%d: ", __FILE__, __LINE__);            \
      fprintf (stderr, __VA_ARGS__);                              \
      fputc ('\n', stderr);                                       \
      fails++;                                                    \
    }                                                             \
  } while (0)

static uint64_t seed = 88172645463325252ull;

static uint64_t rnd64 (void)
{
  seed ^= seed << 13;
  seed ^= seed >> 7;
  seed ^= seed << 17;
  return seed;
}

static uint32_t rnd (uint32_t n)
{
  return (uint32_t) (rnd64 () % n);
}

/* the reference: by the field, then by the position */
static int ref_lo;
static uint64_t ref_m;
static const uint64_t * ref_key;

static int cmp (const void * a, const void * b)
{
  uint32_t i = *(const uint32_t *) a, j = *(const uint32_t *) b;
  uint64_t fi = ref_key[i] >> ref_lo & ref_m;
  uint64_t fj = ref_key[j] >> ref_lo & ref_m;
  return fi < fj ? -1 : fi > fj ? 1 : i < j ? -1 : i > j;
}

static uint64_t mask_of (int lo, int hi)
{
  return hi - lo >= 64 ? ~(uint64_t) 0 : ((uint64_t) 1 << (hi - lo)) - 1;
}

static void * xmalloc (size_t n)
{
  void * p = malloc (n ? n : 1);
  if (!p) {
    fprintf (stderr, "out of memory\n");
    exit (1);
  }
  return p;
}

/* the sort of key[0:n] (ids id0) against the reference */
static void check (const char * name, const uint64_t * key,
                   const uint32_t * id0, uint32_t n, int lo, int hi)
{
  uint32_t * pos = xmalloc (n * sizeof (uint32_t));
  uint64_t * rk = xmalloc (n * sizeof (uint64_t));
  uint32_t * ri = xmalloc (n * sizeof (uint32_t));
  uint64_t * k = xmalloc (n * sizeof (uint64_t));
  uint32_t * v = xmalloc (n * sizeof (uint32_t));
  size_t bytes = hm_sort_scratch (n);
  void * scratch = xmalloc (bytes);
  for (uint32_t i = 0; i < n; i++)
    pos[i] = i;
  ref_lo = lo, ref_m = mask_of (lo, hi), ref_key = key;
  qsort (pos, n, sizeof (uint32_t), cmp);
  for (uint32_t i = 0; i < n; i++)
    rk[i] = key[pos[i]], ri[i] = id0[pos[i]];
  memcpy (k, key, n * sizeof (uint64_t));
  memcpy (v, id0, n * sizeof (uint32_t));
  memset (scratch, 0xa5, bytes);
  hm_sort_keys (k, v, n, lo, hi, scratch);
  CHECK (!memcmp (k, rk, n * sizeof (uint64_t)) &&
         !memcmp (v, ri, n * sizeof (uint32_t)),
         "%s, n %u, bits %d to %d: not the stable order", name, n, lo,
         hi - 1);
  free (pos), free (rk), free (ri), free (k), free (v), free (scratch);
}

enum { RANDOM, EQUAL, SORTED, REVERSED, NEARLY, FAR, FEW, LOWEST, NKINDS };
static const char * kinds[NKINDS] = { "random", "all equal", "sorted",
  "reversed", "nearly sorted", "a few moved far", "few distinct",
  "lowest digit only" };

static void fill (int kind, uint64_t * key, uint32_t n, int lo, int hi)
{
  uint64_t m = mask_of (lo, hi);
  for (uint32_t i = 0; i < n; i++)
    key[i] = rnd64 ();
  switch (kind) {
  case EQUAL:
    for (uint32_t i = 0; i < n; i++)    /* the field equal, the rest not */
      key[i] = (key[i] & ~(m << lo)) | (key[0] & m << lo);
    break;
  case FEW:
    for (uint32_t i = 0; i < n; i++)
      key[i] = (key[i] & ~(m << lo)) | ((uint64_t) rnd (5) & m) << lo;
    break;
  case LOWEST:
    for (uint32_t i = 0; i < n; i++)
      key[i] = (key[i] & ~(m << lo)) | (key[i] & m & 0x7f) << lo;
    break;
  case SORTED: case REVERSED: case NEARLY: case FAR: {
    ref_lo = lo, ref_m = m, ref_key = key;
    uint32_t * p = xmalloc (n * sizeof (uint32_t));
    uint64_t * t = xmalloc (n * sizeof (uint64_t));
    for (uint32_t i = 0; i < n; i++)
      p[i] = i;
    qsort (p, n, sizeof (uint32_t), cmp);
    for (uint32_t i = 0; i < n; i++)
      t[i] = key[p[kind == REVERSED ? n - 1 - i : i]];
    memcpy (key, t, n * sizeof (uint64_t));
    free (p), free (t);
    if (kind == NEARLY)                 /* neighbours swapped */
      for (uint32_t s = 0; n > 1 && s < n / 512 + 1; s++) {
        uint32_t i = rnd (n - 1);
        uint64_t x = key[i];
        key[i] = key[i + 1], key[i + 1] = x;
      }
    if (kind == FAR)                    /* a few moved from end to end */
      for (uint32_t s = 0; n > 1 && s < n / 512 + 1; s++) {
        uint32_t i = rnd (n / 8 + 1), j = n - 1 - rnd (n / 8 + 1);
        uint64_t x = key[i];
        key[i] = key[j], key[j] = x;
      }
    break;
  }
  default:
    break;
  }
}

static double now (void)
{
  struct timespec t;
  clock_gettime (CLOCK_MONOTONIC, &t);
  return (double) t.tv_sec + 1e-9 * (double) t.tv_nsec;
}

/* ns a key of the sort of random keys (the radix), best of 5 */
static double bench (uint32_t n, int lo, int hi)
{
  uint64_t * key0 = xmalloc (n * sizeof (uint64_t));
  uint64_t * k = xmalloc (n * sizeof (uint64_t));
  uint32_t * v = xmalloc (n * sizeof (uint32_t));
  void * scratch = xmalloc (hm_sort_scratch (n));
  for (uint32_t i = 0; i < n; i++)
    key0[i] = rnd64 ();
  double best = 1e30;
  for (int r = 0; r < 6; r++) {        /* the first a warm up */
    memcpy (k, key0, n * sizeof (uint64_t));
    for (uint32_t i = 0; i < n; i++)
      v[i] = i;
    double t = now ();
    hm_sort_keys (k, v, n, lo, hi, scratch);
    t = now () - t;
    if (r > 0 && t < best)
      best = t;
  }
  free (key0), free (k), free (v), free (scratch);
  return 1e9 * best / n;
}

int main (void)
{
  const int B = HM_DIMENSION * HM_MORTON_BITS;
  /* bits lo to hi - 1: whole keys, the cells of 2^10 and 2^4 a side,
     ranges inside (bits above hi set) */
  int ranges[][2] = {
    { 0, B }, { HM_DIMENSION * (HM_MORTON_BITS - 10), B },
    { HM_DIMENSION * (HM_MORTON_BITS - 4), B }, { 13, 37 }, { 5, 6 },
    { B - 1, B }, { 8, 24 }
  };
  uint32_t sizes[] = { 0, 1, 2, 3, 255, 256, 257, 511, 2048, 2049, 4097,
                       65537 };
  int nr = (int) (sizeof (ranges) / sizeof (ranges[0]));
  int ns = (int) (sizeof (sizes) / sizeof (sizes[0]));
  for (int r = 0; r < nr; r++)
    for (int s = 0; s < ns; s++)
      for (int kind = 0; kind < NKINDS; kind++) {
        uint32_t n = sizes[s];
        uint64_t * key = xmalloc (n * sizeof (uint64_t));
        uint32_t * id = xmalloc (n * sizeof (uint32_t));
        fill (kind, key, n, ranges[r][0], ranges[r][1]);
        for (uint32_t i = 0; i < n; i++)  /* ids not in order */
          id[i] = (uint32_t) (rnd64 () >> 32);
        check (kinds[kind], key, id, n, ranges[r][0], ranges[r][1]);
        free (key), free (id);
      }

  if (getenv ("HM_SORT_BENCH")) {
    printf ("the sort, ns a key, best of 5:\n");
    printf ("  %-9s %-14s %10s\n", "n", "bits", "ns a key");
    uint32_t bn[] = { 1u << 20, 1u << 22 };
    int br[][2] = { { HM_DIMENSION * (HM_MORTON_BITS - 10), B }, { 0, B } };
    for (int i = 0; i < 2; i++)
      for (int j = 0; j < 2; j++) {
        char bits[32];
        snprintf (bits, sizeof (bits), "%d to %d", br[j][0], br[j][1] - 1);
        printf ("  2^%-7d %-14s %10.2f\n", i ? 22 : 20, bits,
                bench (bn[i], br[j][0], br[j][1]));
      }
  }

  if (fails)
    printf ("%d of %d checks failed\n", fails, checks);
  else
    printf ("%d checks passed\n", checks);
  return fails != 0;
}
