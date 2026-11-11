#include <stdio.h>
#include "stencil.h"

scalar s, t, r, q;

static double sum_of (scalar a)
{
  double sum = 0;
  foreach ()
    sum += a[];
  return sum;
}

int main (void)
{
  /* t: its storage index, ghosts included */
  foreach_leaf () {
    Real * a = block_data (sblocks, t);
    for (int k = 0; k < NX * NX; k++)
      a [k] = k + 100 * i;
  }

  /* the 3 x 3 box around each cell, against its offsets written out */
  foreach () {
    double sum = 0;
    foreach_neighbor (1)
      sum += t[];
    r[] = sum - (t[-1,-1] + t[0,-1] + t[1,-1] + t[-1,0] + t[] + t[1,0]
                 + t[-1,1] + t[0,1] + t[1,1]);
  }

  /* an offset in it adds to the neighbour's; x and Delta are the
   * neighbour's: their sum over the box is 9 times the cell's */
  foreach () {
    double shifted = 0, xs = 0, ds = 0;
    foreach_neighbor (1) {
      shifted += t[1,0];
      xs += x;
      ds += Delta;
    }
    s[] = shifted - (t[0,-1] + t[1,-1] + t[2,-1] + t[0,0] + t[1,0] + t[2,0]
                     + t[0,1] + t[1,1] + t[2,1]) + xs - 9 * x + ds - 9 * Delta;
  }

  /* q written, then read by foreach_neighbor (): its ghosts updated first */
  foreach ()
    q[] = 1;
  foreach () {
    double n = 0;
    foreach_neighbor (1)
      n += q[];
    r[] += n;
  }
  printf ("%g %g | %s\n", sum_of (r), sum_of (s), ghosts_log);
  return 0;
}
