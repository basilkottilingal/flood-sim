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
      a [k] = k + 100 * i;
  }
  int n = 0;
  foreach () {
    s[] = t[1] - t[-1] + t[0,1] - t[0,-1];  /* 2 + 12 */
    n++;
  }
  double c = 0.5;
  foreach ()
    u.x[] = c + t[] - t[];                  /* the body's c, not the loop's */
  printf ("%d %g %g", n, sum_of (s), sum_of (u.x));

  /* x, y and Delta are the cell's, an outer x is hidden, the body's own is
   * its own */
  double x = 99;
  foreach () {
    s[] = x;
    u.x[] = y + Delta;
  }
  foreach () {
    double x = 7;
    u.y[] = x;
  }
  printf (" %g %g %g %g\n", sum_of (s), sum_of (u.x), sum_of (u.y), x);
  return 0;
}
