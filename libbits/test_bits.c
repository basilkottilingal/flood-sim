/* bits.h against naive loops over the bits: u8 and u16 exhaustively, u32
 * and u64 at the edges (0, every single bit, every 2^k - 1, 2^k + 1, the
 * maximum) and at random. C89 style, so that it builds as gnu89 too. */
#include <stdio.h>
#include "bits.h"

static long checks, fails;

/* the naive answers, for x of n bits */
static int ref_popcount (uint64_t x, int n)
{
  int k, c = 0;
  for (k = 0; k < n; k++)
    c += (int) ((x >> k) & 1u);
  return c;
}

static int ref_ctz (uint64_t x, int n)
{
  int k;
  for (k = 0; k < n; k++)
    if ((x >> k) & 1u)
      return k;
  return n;
}

static int ref_clz (uint64_t x, int n)
{
  int k;
  for (k = n - 1; k >= 0; k--)
    if ((x >> k) & 1u)
      return n - 1 - k;
  return n;
}

static int ref_highest (uint64_t x, int n)
{
  return n - 1 - ref_clz (x, n);
}

static int ref_log2_ceil (uint64_t x, int n)
{
  int k;
  if (!x)
    return -1;
  for (k = 0; k < n; k++)
    if ((uint64_t) 1 << k >= x)
      return k;
  return n;                       /* x above 2^(n-1) */
}

static uint64_t ref_ceil_pow2 (uint64_t x, int n)
{
  int k = ref_log2_ceil (x, n);
  if (k <= 0)
    return 1;
  return k >= n ? 0 : (uint64_t) 1 << k;
}

#define CHECK(c, what, x) do { checks++; if (!(c)) { fails++;                  \
  if (fails < 20) printf ("FAIL %s at 0x%llx\n", what, (unsigned long long) (x)); } } while (0)

/* every function of width N at x, against the naive ones */
#define TEST(N, x0) do {                                                       \
  uint##N##_t x = (uint##N##_t) (x0), y = x;                                   \
  uint64_t w = (uint64_t) x;                                                   \
  int b;                                                                       \
  CHECK (bits_popcount_u##N (x) == ref_popcount (w, N), "popcount" #N, w);     \
  CHECK (bits_ctz_u##N (x) == ref_ctz (w, N), "ctz" #N, w);                    \
  CHECK (bits_clz_u##N (x) == ref_clz (w, N), "clz" #N, w);                    \
  CHECK (bits_lowest_u##N (x) == (w ? ref_ctz (w, N) : -1), "lowest" #N, w);   \
  CHECK (bits_highest_u##N (x) == ref_highest (w, N), "highest" #N, w);        \
  CHECK (bits_log2_floor_u##N (x) == ref_highest (w, N), "log2_floor" #N, w);  \
  CHECK (bits_log2_ceil_u##N (x) == ref_log2_ceil (w, N), "log2_ceil" #N, w);  \
  CHECK (bits_is_pow2_u##N (x) == (ref_popcount (w, N) == 1), "is_pow2" #N, w);\
  CHECK ((uint64_t) bits_ceil_pow2_u##N (x) == ref_ceil_pow2 (w, N),           \
         "ceil_pow2" #N, w);                                                   \
  CHECK (bits_parity_u##N (x) == (ref_popcount (w, N) & 1), "parity" #N, w);   \
  b = bits_pop_lowest_u##N (&y);                                               \
  CHECK (b == (w ? ref_ctz (w, N) : -1) &&                                     \
         (uint64_t) y == (w ? w & (w - 1) : 0), "pop_lowest" #N, w);           \
} while (0)

/* xorshift64: random values, spread over every magnitude by a shift */
static uint64_t state = UINT64_C (0x9E3779B97F4A7C15);
static uint64_t rnd (void)
{
  state ^= state << 13;
  state ^= state >> 7;
  state ^= state << 17;
  return state >> (state & 63);
}

int main (void)
{
  long i;
  int k;
  uint64_t r;

  /* 8 and 16 bits: every value */
  for (i = 0; i < 256; i++)
    TEST (8, i);
  for (i = 0; i < 65536; i++)
    TEST (16, i);

  /* 32 and 64 bits: the edges */
  TEST (32, 0);
  TEST (64, 0);
  TEST (32, UINT32_MAX);
  TEST (64, UINT64_MAX);
  for (k = 0; k < 64; k++) {
    r = (uint64_t) 1 << k;
    TEST (32, r);
    TEST (32, r - 1);
    TEST (32, r + 1);
    TEST (64, r);
    TEST (64, r - 1);
    TEST (64, r + 1);
    TEST (64, ~r);
  }

  /* and at random */
  for (i = 0; i < 2000000; i++) {
    r = rnd ();
    TEST (32, r);
    TEST (64, r);
    TEST (16, r);
    TEST (8, r);
  }

  /* a loop over the set bits */
  {
    uint64_t x = UINT64_C (0x8000000100000005);
    int seen [4], n = 0, b;
    while ((b = bits_pop_lowest_u64 (&x)) >= 0 && n < 4)
      seen [n++] = b;
    CHECK (n == 4 && seen [0] == 0 && seen [1] == 2 && seen [2] == 32
           && seen [3] == 63 && x == 0, "set bits loop", x);
  }

#if defined (__STDC_VERSION__) && __STDC_VERSION__ >= 201112L
  /* C11: the names without suffix, by type */
  {
    uint8_t a = 0x10;
    uint16_t c = 0x0100;
    uint32_t d = 0x00010000u;
    uint64_t e = UINT64_C (1) << 40;
    unsigned long long f = 1;
    CHECK (bits_clz (a) == 3 && bits_clz (c) == 7 && bits_clz (d) == 15
           && bits_clz (e) == 23, "generic clz", 0);
    CHECK (bits_highest (a) == 4 && bits_lowest (e) == 40
           && bits_popcount (d) == 1 && bits_log2_ceil (c) == 8, "generic", 0);
    CHECK (bits_clz (f) == 63 && bits_clz (1) == 63, "generic default", 0);
    CHECK (bits_ceil_pow2 ((uint8_t) 129) == 0 && bits_ceil_pow2 (c) == c,
           "generic ceil_pow2", 0);
    CHECK (bits_pop_lowest (&a) == 4 && a == 0 && bits_pop_lowest (&a) == -1,
           "generic pop_lowest", 0);
  }
#endif

  printf ("%s, %s: %ld checks, %ld failed\n",
          BITS_HAVE_BUILTINS ? "builtins" : "de Bruijn",
#if defined (__STDC_VERSION__)
          __STDC_VERSION__ >= 201112L ? "C11 or later" : "C99",
#else
          "C89",
#endif
          checks, fails);
  return fails != 0;
}
