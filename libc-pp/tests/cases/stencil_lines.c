#include "stencil.h"
scalar s;
void f (void)
{
  foreach () {
    s[] = 1;
    s[] = nope1;
  }
  int x = nope2;
  (void) x;
}
