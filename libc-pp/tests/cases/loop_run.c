#include <stdio.h>
#include "loop.h"

int other (void);

int main (void)
{
  /* names the loops also use: theirs are hidden */
  int count = 100, idx = -1, depth = 7, h = 9, n = 0, total = 0, sums = 0;
  int in_use = 0;
  repeat (3) {
    n += k + 1;
  }
  range (0, 10, 3) {
    total += i;
    count++;
  }
  Bag bag = { { 5, 6, 7 }, 3 };
  Bag * b = &bag;
  each (b) {
    total += x * 100;
    if (x == 6)
      break;
  }
  repeat (2) repeat (2) n += 10;
  walk ((P) { 3, 0 }) { total += level * 1000 + q.depth; }
  walk2 (4) { sums += sum; }
  squares (4) { total += v * 10000; }
  walk ((P) { 2, 0 }) {
    scalar s;
    in_use = scalar_count ();
    (void) s;
    (void) q;
  }
  printf ("%d %d %d %d %d %d %d %d %d %d\n", n, total, sums, count, idx,
          depth, h, in_use, scalar_count (), other ());
  return 0;
}
