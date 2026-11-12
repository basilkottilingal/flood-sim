#include <stdio.h>
#include "stencil.h"

scalar p, q;
enum { ONE = 1 };

int main (void)
{
  foreach ()
    p[] = 0;
  /* the red cells, from their black neighbours: no warning */
  foreach_mock_colour (0)
    p[] = p[1] + p[-1] + p[0,1] + p[0,-ONE] + 1;
  /* an odd sum of offsets reaches the other colour too */
  foreach_mock_colour (1)
    p[] = p[2,1] + p[ONE,0] + 2;
  int red = 0, black = 0, other = 0;
  foreach () {
    if (p[] == 1)
      red++;
    else if (p[] >= 2)
      black++;
    else
      other++;
  }
  printf ("red %d black %d other %d\n", red, black, other);
  return 0;
}
