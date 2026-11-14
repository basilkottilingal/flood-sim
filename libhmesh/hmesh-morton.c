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

/*
.. Notes on the radix sort
.. -----------------------
..
.. LSD (least significant digit first). Each pass is a counting sort on one
.. digit: count the keys of each digit value, turn the counts into starting
.. offsets (an exclusive prefix sum), then scatter every (key, id) pair to
.. its offset. A counting sort is stable, so after the pass on digit d the
.. keys are in order by digits 0 ... d, which is why the passes go from the
.. lowest digit up and why the whole sort is stable (Knuth, 5.2.5). Its
.. cost is the same whatever the input, (hi - lo) / RADIX_BITS passes of
.. two reads and a write per key, with no comparison and no branch on the
.. data: it does not get worse on adversarial input, and does not get
.. better on sorted input either, which is what the descents check in
.. hm_sort_keys () is for.
..
.. Why stable matters here, twice over:
..   - the passes compose (above);
..   - a sort on bits lo ... hi - 1 only (the tree's cells) keeps the
..     order the points had within a cell: the last cache's, so that the
..     next sort starts from an order nearly right.
..
.. The digit: 11 bits. The counts are 2^11 4 byte words, 8 KiB, which stay
.. in L1 while the keys stream past; 64 bits take 6 passes (8 with 8 bit
.. digits, 4 with 16 bit, whose 256 KiB of counts would not stay in L1).
.. This is the choice of Herf's "Radix Tricks" for 32 bit floats (three
.. passes of 11 bits). Measured here: about 4 to 5 ns a key a pass for 2^20
.. keys, memory bound.
..
.. A pass is skipped when every key has the same digit (the count of the
.. first key's digit is n): the high bits of keys in a domain mostly empty,
.. or a sort to coarse cells. Terdiman ("Radix Sort Revisited") also checks
.. for sorted input before sorting; here hm_sort_keys () does, and tries
.. insertion for a few descents.
..
.. Two buffers, swapped each pass ("ping-pong"); after an odd number of
.. passes the result is in the second, and is copied back. The histogram of
.. every digit could be made in one read of the keys before the first pass
.. (Terdiman, Herf), a further saving of up to a pass of reads; not done.
..
.. What else could be done:
..   - MSD radix, recursing in each bucket, stops early where buckets get
..     small and suits keys of many bits; it is not stable without care,
..     and recursion costs where buckets are large.
..   - On a GPU, or in parallel, the same LSD scheme with per-block
..     histograms and a scan: Merrill and Grimshaw (2011), CUB's
..     DeviceRadixSort.
..   - For the leaves alone, the points could be counted into the leaves
..     (one counting sort keyed by the leaf, nleaves buckets) instead of
..     sorted along the whole curve; the curve order is kept as it also
..     gives the order within a leaf and the duplicates.
..
.. References
..   D. E. Knuth, The Art of Computer Programming, Vol. 3, Sorting and
..     Searching, 2nd ed., Addison-Wesley, 1998, 5.2.5 (sorting by
..     distribution).
..   P. Terdiman, Radix Sort Revisited, 2000 (codercorner.com).
..   M. Herf, Radix Tricks, 2001 (stereopsis.com/radix.html).
..   D. Merrill, A. Grimshaw, High performance and scalable radix sorting,
..     Parallel Processing Letters 21 (2), 2011.
..   G. M. Morton, A computer oriented geodetic data base and a new
..     technique in file sequencing, IBM, 1966: the curve.
..   H. Sundar, R. S. Sampath, G. Biros, Bottom-up construction and 2:1
..     balance refinement of linear octrees in parallel, SIAM J. Sci.
..     Comput. 30 (5), 2008: an octree as its leaves' Morton keys, sorted.
..   T. Karras, Maximizing parallelism in the construction of BVHs, octrees,
..     and k-d trees, High Performance Graphics, 2012: the tree implicit in
..     sorted keys (each node a common prefix), the "one point a leaf" tree
..     of hmesh-tree.h.
..   S. E. Anderson, Bit Twiddling Hacks (graphics.stanford.edu/~seander/
..     bithacks.html): the interleaving by masks of spread () above.
*/

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
