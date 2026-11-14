/*
.. Morton keys and their sort: see hmesh-morton.h.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "hmesh-morton.h"

#if HM_DIMENSION == 2
/* the low 32 bits of x to the even bits */
static uint64_t spread (uint64_t x)
{
  x &= 0xffffffffu;
  x = (x | x << 16) & 0x0000ffff0000ffffull;
  x = (x | x << 8)  & 0x00ff00ff00ff00ffull;
  x = (x | x << 4)  & 0x0f0f0f0f0f0f0f0full;
  x = (x | x << 2)  & 0x3333333333333333ull;
  x = (x | x << 1)  & 0x5555555555555555ull;
  return x;
}
#elif HM_DIMENSION == 3
/* the low 21 bits of x to every third bit */
static uint64_t spread (uint64_t x)
{
  x &= 0x1fffff;
  x = (x | x << 32) & 0x001f00000000ffffull;
  x = (x | x << 16) & 0x001f0000ff0000ffull;
  x = (x | x << 8)  & 0x100f00f00f00f00full;
  x = (x | x << 4)  & 0x10c30c30c30c30c3ull;
  x = (x | x << 2)  & 0x1249249249249249ull;
  return x;
}
#endif

uint64_t hm_morton_encode (const uint64_t * q)
{
#if HM_DIMENSION == 1
  return q[0];
#elif HM_DIMENSION == 2
  return spread (q[0]) | spread (q[1]) << 1;
#else
  return spread (q[0]) | spread (q[1]) << 1 | spread (q[2]) << 2;
#endif
}

uint64_t hm_morton_key (const double * x, const double * origin, double size)
{
  const double s = (double) ((uint64_t) 1 << HM_MORTON_BITS);
  const uint64_t top = ((uint64_t) 1 << HM_MORTON_BITS) - 1;
  uint64_t q[HM_DIMENSION];
  for (int d = 0; d < HM_DIMENSION; d++) {
    double u = (x[d] - origin[d]) / size * s;
    q[d] = !(u > 0.) ? 0 : u >= s ? top : (uint64_t) u;
  }
  return hm_morton_encode (q);
}

uint64_t hm_morton_cell (const int * ijk, int level)
{
  uint64_t q[HM_DIMENSION];
  for (int d = 0; d < HM_DIMENSION; d++)
    q[d] = (uint64_t) ijk[d] << (HM_MORTON_BITS - level);
  return hm_morton_encode (q);
}

/* insertion by the bits from lo up, giving up (0) once it has moved budget
   entries */
static int insertion (uint64_t * key, uint32_t * id, uint32_t n, int lo,
                      uint64_t budget)
{
  uint64_t moved = 0;
  for (uint32_t i = 1; i < n; i++) {
    uint64_t k = key[i];
    uint32_t v = id[i], j = i;
    while (j > 0 && key[j - 1] >> lo > k >> lo) {
      key[j] = key[j - 1];
      id[j] = id[j - 1];
      j--;
      if (++moved > budget) {
        key[j] = k;            /* still a permutation: radix goes on */
        id[j] = v;
        return 0;
      }
    }
    key[j] = k;
    id[j] = v;
  }
  return 1;
}

/* 11 bits a digit: 6 passes over 64 bit keys, the counts in L1 */
#define RADIX_BITS 11
#define RADIX      (1 << RADIX_BITS)

static void radix (uint64_t * key, uint32_t * id, uint32_t n, int lo,
                   int hi)
{
  uint64_t * k2 = malloc ((size_t) n * sizeof (uint64_t));
  uint32_t * i2 = malloc ((size_t) n * sizeof (uint32_t));
  if (!k2 || !i2) {
    fprintf (stderr, "hmesh: out of memory\n");
    abort ();
  }
  uint64_t * k = key;
  uint32_t * v = id;
  for (int shift = lo; shift < hi; shift += RADIX_BITS) {
    /* the digit's bits below hi only */
    uint64_t mask = hi - shift >= RADIX_BITS ? RADIX - 1 :
      ((uint64_t) 1 << (hi - shift)) - 1;
    uint32_t c[RADIX] = { 0 };
    for (uint32_t i = 0; i < n; i++)
      c[(k[i] >> shift) & mask]++;
    if (c[(k[0] >> shift) & mask] == n)
      continue;
    uint32_t s = 0;
    for (int d = 0; d < RADIX; d++) {
      uint32_t t = c[d];
      c[d] = s;
      s += t;
    }
    for (uint32_t i = 0; i < n; i++) {
      uint32_t p = c[(k[i] >> shift) & mask]++;
      k2[p] = k[i];
      i2[p] = v[i];
    }
    uint64_t * tk = k;
    uint32_t * tv = v;
    k = k2, v = i2, k2 = tk, i2 = tv;
  }
  if (k != key) {
    memcpy (key, k, (size_t) n * sizeof (uint64_t));
    memcpy (id, v, (size_t) n * sizeof (uint32_t));
    free (k), free (v);
  }
  else
    free (k2), free (i2);
}

void hm_sort_keys (uint64_t * key, uint32_t * id, uint32_t n, int lo,
                   int hi)
{
  /* the descents: none, sorted; a few, insertion (which gives up if they
     go far); else radix */
  uint32_t descents = 0;
  for (uint32_t i = 1; i < n; i++)
    descents += key[i] >> lo < key[i - 1] >> lo;
  if (!descents)
    return;
  if (descents <= n / 256 && insertion (key, id, n, lo, 2 * (uint64_t) n))
    return;
  radix (key, id, n, lo, hi);
}
