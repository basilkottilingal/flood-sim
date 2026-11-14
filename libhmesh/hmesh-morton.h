#ifndef _HMESH_MORTON_H_
#define _HMESH_MORTON_H_

  /*
  .. Morton keys of points in a box, and their sort.
  ..
  .. A point x of the box [origin, origin + size]^D is quantised to
  .. HM_MORTON_BITS bits an axis, q = floor ((x - origin) / size 2^b),
  .. q = 2^b - 1 on the upper faces, and its key interleaves the bits of q,
  .. x in the lowest: the order of libb-amr's children (index_child (): bit
  .. 0 of c along x, 1 along y, 2 along z), so that sorted keys follow the
  .. tree's leaves (tree_caches ()->leaves). b is 52 in 1D (a double's
  .. mantissa), 32 in 2D and 21 in 3D (64 bits).
  ..
  .. The node of level l at (i, j, k) covers the keys whose top D l bits are
  .. those of hm_morton_cell (ijk, l), a range of 2^(D (b - l)) keys.
  ..
  .. hm_sort_keys (key, id, n, lo, hi): key[] ascending by its bits lo to
  .. hi - 1, id[] along, the order of keys equal there kept (stable): with
  .. lo = D (b - c), hi = D b, the order of the cells of 2^c a side, which
  .. is all the tree needs and costs (hi - lo) / 11 passes. The sort chooses
  .. by what it is given:
  ..
  ..   sorted already          one pass that counts the descents (a key
  ..                           below the one before), about 0.6 ns a key
  ..   a descent in 256 or     insertion, while it moves fewer than 2 n
  ..   fewer                   entries in all (a few points added at the
  ..                           end, a few moved)
  ..   anything else           LSD radix, 11 bits a pass from bit lo, a pass
  ..                           skipped where every key has the same digit:
  ..                           6 passes for 64 bits
  ..
  .. Points that all moved, even a little, are not nearly sorted at the
  .. keys' resolution (the curve's jumps), and take the radix sort.
  */

  #include <stdint.h>
  #include "hmesh.h"

  #if HM_DIMENSION == 1
    #define HM_MORTON_BITS 52
  #elif HM_DIMENSION == 2
    #define HM_MORTON_BITS 32
  #else
    #define HM_MORTON_BITS 21
  #endif

  uint64_t hm_morton_encode (const uint64_t * q);
  uint64_t hm_morton_key (const double * x, const double * origin,
                          double size);
  uint64_t hm_morton_cell (const int * ijk, int level);
  void     hm_sort_keys (uint64_t * key, uint32_t * id, uint32_t n, int lo,
                         int hi);

#endif
