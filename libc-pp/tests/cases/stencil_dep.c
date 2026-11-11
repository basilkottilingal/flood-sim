#include "stencil.h"
scalar p, q;
void f (void)
{
  foreach () {
    p[] = x;
    q[] = p[1,0] - p[-1,0];
  }
}
