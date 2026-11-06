/* Loops for the nesting tests: loops in loops' bodies, up to three levels */
#include <cwrap/cache.h>

@loop repeat (int n) : k
{
  for (int k = 0; k < n; k++)
    {...}
}

@loop range (int lo, int hi) : i
{
  int count = 0;
  for (int i = lo; i < hi; i++) {
    count++;
    {...}
  }
}

/* the same inner loop nested in itself: the inner i hides the outer, which
 * is saved first */
@loop grid (int n, int m) : r, c
{
  range (0, n) {
    int r = i;
    range (0, m) {
      int c = i;
      {...}
    }
  }
}

/* three levels: cube uses grid, which uses range. r and c come from grid,
 * s is cube's own, and range's i stays hidden */
@loop cube (int n) : r, c, s
{
  range (0, n) {
    int s = i;
    grid (n, n) {
      {...}
    }
  }
}

/* a name passed on from two levels down */
@loop cube_r (int n) : r
{
  cube (n) {
    (void) c;
    (void) s;
    {...}
  }
}

/* inner loops before {...}, around it and after it */
@loop pre_post (int n) : k
{
  int before = 0;
  range (0, n)
    before += i;
  repeat (before / (n * (n - 1) / 2))
    {...}
  range (0, before)
    (void) i;
}

/* a break in the inner of two nested inner loops */
@loop pairs (int n) : a, b
{
  range (0, n) {
    int a = i;
    range (0, n) {
      if (i > a)
        break;
      int b = i;
      {...}
    }
  }
}

typedef struct { int v [6]; int n; } Bag;

/* a parameter the block sees */
@loop each (const Bag * bag) : x, bag
{
  for (int idx = 0; idx < bag->n; idx++) {
    int x = bag->v [idx];
    {...}
  }
}

/* a static local: one per use */
@loop counted (int n) : k, calls
{
  static int calls;
  calls++;
  repeat (n)
    {...}
}

/* a cache and a foreach */
@loop squares (int n) : v
{
  cache (int) sq;
  range (0, n)
    cache_push (sq, i * i);
  foreach (v in sq)
    {...}
}

/* a scalar of the loop's own, given to the block */
@loop with_scalar (int n) : s, k
{
  scalar s;
  repeat (n)
    {...}
}
