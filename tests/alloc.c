#include "grid.h"
#include <stdio.h>

int main ()
{
  int ** g = grid (int, 100, 105);
  GridData gd = grid_data (g);
  printf ("Grid {%d %d %d, %d}", gd.width, gd.height, gd.nghosts, (int) gd.size);
  grid_free (g);
}
