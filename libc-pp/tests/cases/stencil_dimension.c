#include <stdio.h>
#include <math.h>
#include "stencil.h"

scalar p, q;
vector u;

int main (void)
{
  foreach ()
    p[] = (Real) (x + 10 * y);
  /* the gradient, written for x: u.x from p[1] - p[-1], u.y from p[0,1] -
   * p[0,-1]; read at cells away from the block's edges (no ghosts here) */
  foreach ()
    foreach_dimension ()
      u.x[] = (p[1] - p[-1]) / (2 * Delta);
  int good = 0, bad = 0;
  foreach () {
    int ix = (int) floor (x / Delta), iy = (int) floor (y / Delta);
    if (ix % 4 == 0 || ix % 4 == 3 || iy % 4 == 0 || iy % 4 == 3)
      continue;
    if (fabs (u.x[] - 1) < 1e-12 && fabs (u.y[] - 10) < 1e-12)
      good++;
    else
      bad++;
  }
  /* names ending in _x, and x, y themselves */
  double g_x = 1, g_y = 2;
  foreach_dimension ()
    g_x *= 3;
  double worst = 0;
  foreach () {
    double s = 0;
    foreach_dimension ()
      s += x;
    q[] = (Real) s;
    worst = fmax (worst, fabs (q[] - (x + y)));
  }
  printf ("gradient %d good %d bad; g %g %g; coordinates %g\n", good, bad,
          g_x, g_y, worst);
  return 0;
}
