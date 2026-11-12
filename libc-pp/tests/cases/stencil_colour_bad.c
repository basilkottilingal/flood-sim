#include "stencil.h"
scalar p, q;
enum { TWO = 2 };
void f (void)
{
  foreach_mock_colour (0)
    p[] = p[1,1];
  foreach_mock_colour (1)
    p[] = p[TWO];
  foreach_mock_colour (0)
    p[] = q[1,1] + p[1];
}
