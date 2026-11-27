/* tree_morton_key (): the interleaving by masks against the reference, one
 * bit at a time (the code it replaced). Exhaustive along each axis: every
 * value of one coordinate at the finest level, the others 0 (the key is the
 * OR of one part an axis, a copy of each bit); then every index of every
 * level whose whole grid has at most 2^18 cells; then random indices at
 * random levels, every axis at once. And the key of a child from its
 * parent's: key (p) | c << dimension * (BAMR_MAX_DEPTH - level - 1). */
#include <stdio.h>
#include <stdint.h>
#include "tree.h"

static int checks, fails;
#define CHECK(c) do { checks++; if (!(c)) { fails++; \
  printf ("FAIL line %d: %s\n", __LINE__, #c); } } while (0)

/* the reference: bit b of axis a to bit dimension * b + a */
static uint64_t key_ref (Index p)
{
  int s = BAMR_MAX_DEPTH - p.level;
  uint64_t key = 0, x [3] = { (uint64_t) p.i << s, 0, 0 };
  #if dimension > 1
  x [1] = (uint64_t) p.j << s;
  #endif
  #if dimension > 2
  x [2] = (uint64_t) p.k << s;
  #endif
  for (int b = 0; b < BAMR_MAX_DEPTH; b++)
    for (int a = 0; a < dimension; a++)
      key |= ((x [a] >> b) & 1) << (dimension * b + a);
  return key;
}

static Index make (int level, const int64_t x [3])
{
  Index p = { 0 };
  p.level = level;
  p.i = (int) x [0];
  #if dimension > 1
  p.j = (int) x [1];
  #endif
  #if dimension > 2
  p.k = (int) x [2];
  #endif
  return p;
}

static uint64_t rng = 88172645463325252ull;   /* xorshift64 */
static uint64_t next (void)
{
  rng ^= rng << 13, rng ^= rng >> 7, rng ^= rng << 17;
  return rng;
}

int main (void)
{
  /* 1. each axis alone, every value at the finest level */
  long wrong = 0;
  for (int a = 0; a < dimension; a++)
    for (int64_t v = 0; v < (int64_t) 1 << BAMR_MAX_DEPTH; v++)
    {
      int64_t x [3] = { 0, 0, 0 };
      x [a] = v;
      Index p = make (BAMR_MAX_DEPTH, x);
      wrong += tree_morton_key (p) != key_ref (p);
    }
  CHECK (wrong == 0);

  /* 2. every index of each level with at most 2^18 cells */
  wrong = 0;
  long cells = 0;
  for (int l = 0; l <= BAMR_MAX_DEPTH && dimension * l <= 18; l++)
  {
    int64_t n = (int64_t) 1 << l, x [3] = { 0, 0, 0 };
    for (x [0] = 0; x [0] < n; x [0]++)
    #if dimension > 1
    for (x [1] = 0; x [1] < n; x [1]++)
    #endif
    #if dimension > 2
    for (x [2] = 0; x [2] < n; x [2]++)
    #endif
    {
      Index p = make (l, x);
      wrong += tree_morton_key (p) != key_ref (p);
      cells++;
    }
  }
  CHECK (wrong == 0);

  /* 3. random indices, all axes at once, and the child's key */
  wrong = 0;
  long family = 0;
  for (int t = 0; t < 1000000; t++)
  {
    int l = (int) (next () % (BAMR_MAX_DEPTH + 1));
    int64_t x [3] = { 0, 0, 0 };
    for (int a = 0; a < dimension; a++)
      x [a] = (int64_t) (next () & (((uint64_t) 1 << l) - 1));
    Index p = make (l, x);
    uint64_t k = tree_morton_key (p);
    wrong += k != key_ref (p);
    if (l < BAMR_MAX_DEPTH)
    {
      int c = (int) (next () % TREE_NCHILDREN);
      int64_t y [3] = { 2 * x [0] + (c & 1), 0, 0 };
      #if dimension > 1
      y [1] = 2 * x [1] + ((c >> 1) & 1);
      #endif
      #if dimension > 2
      y [2] = 2 * x [2] + ((c >> 2) & 1);
      #endif
      family += tree_morton_key (make (l + 1, y))
        != (k | (uint64_t) c << (dimension * (BAMR_MAX_DEPTH - l - 1)));
    }
  }
  CHECK (wrong == 0);
  CHECK (family == 0);

  printf ("dimension %d: %ld cells of the small levels, 2^%d values an "
          "axis, 10^6 random; %d checks, %d failed\n", dimension, cells,
          BAMR_MAX_DEPTH, checks, fails);
  return fails != 0;
}
