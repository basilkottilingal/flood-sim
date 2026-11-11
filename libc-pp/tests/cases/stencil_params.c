#include <stdio.h>
#include "stencil.h"

/* the leaf i = l alone */
@loop leaf_at (int l) : i, j, level, sblocks
{
  for (int i = l; i <= l; i++) {
    int j = 0, level = 2;
    (void) j, (void) level;
    void *** sblocks = leaf_blocks (i);
    {...}
  }
}

/* a stencil loop with a parameter, given to its block loop */
@foreach foreach_in_leaf (int l) = leaf_at (l) foreach_cell (tree_block_shape (), 0);

scalar s;

int main (void)
{
  int calls = 0;
  foreach ()
    s[] = 0;
  foreach_in_leaf (1 + 0 * calls++)     /* the argument, evaluated once */
    s[] = x;
  double sum = 0;
  foreach ()
    sum += s[];
  printf ("%g %d\n", sum, calls);       /* x of the leaf i = 1 only */
  return 0;
}
