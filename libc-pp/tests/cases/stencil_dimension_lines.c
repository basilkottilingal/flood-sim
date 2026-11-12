#include "stencil.h"
scalar p;
vector u;
int f (void)
{
  foreach ()
    foreach_dimension ()
      u.x[] = p[1] -
              p[-1];
  int unused;
  return 0;
}
