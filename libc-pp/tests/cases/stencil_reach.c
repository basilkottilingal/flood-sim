#include "stencil.h"
scalar s, t;
void f (int n)
{
  foreach () s[] = t[3];
  foreach () s[] = t[n];
  foreach () { double a = 0; foreach_neighbor (3) a += t[]; s[] = a; }
  foreach () { double a = 0; foreach_neighbor (2) a += t[0,1]; s[] = a; }
  foreach () { double a = 0; foreach_neighbor (1) a += t[1,-1]; s[] = a + t[-2,2]; }
  foreach_up () s[] = coarse (t, 2) + coarse (t, 0, -1);
}
