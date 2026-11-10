#include "stencil.h"
scalar s;
int a [4];
void f (void)
{
  s[] = 1;
  foreach () a[] = 1;
  foreach () s[1,2,3] = 0;
  foreach () foreach () s[] = 0;
  foreach () a[1] = s[];
}
scalar Delta;
void g (void) { foreach () Delta[] = 0; }
