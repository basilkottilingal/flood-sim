#include <stdio.h>
#include "stencil.h"

scalar s, t;
vector u;

static double sum_of (scalar a)
{
  double sum = 0;
  foreach ()
    sum += a[];
  return sum;
}

int main (void)
{
  /* t: its storage index, and 100 per leaf, ghosts included */
  foreach_leaf () {
    Real * a = block_data (sblocks, t);
    for (int k = 0; k < NX * NX; k++)
      a [k] = k + 100 * level;
  }
  int n = 0;
  foreach () {
    s[] = t[1] - t[-1] + t[0,1] - t[0,-1];  /* 2 + 12 */
    n++;
  }
  double c = 0.5;
  foreach ()
    u.x[] = c + t[] - t[];                  /* the body's c, not the loop's */
  printf ("%d %g %g\n", n, sum_of (s), sum_of (u.x));
  return 0;
}
