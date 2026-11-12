#include "stencil.h"
scalar p, q;
void f (void)
{
  foreach_mock_level (2)
    p[1] = 0;
  foreach_mock_level (2) {
    p[] = 1;
    q[] = p[1];
  }
}
