#include "stencil.h"
scalar p;
void f (void)
{
  foreach_mock_colour (0) {
    double s = 0;
    foreach_neighbor (1)
      s += p[];
    p[] = s;
  }
}
