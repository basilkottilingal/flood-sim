#include "stencil.h"
scalar s, t;
void f (void)
{
  foreach_up () s[] = t[1];
  foreach_up () coarse (s) = 1;
  foreach_up () foreach_child () s[] += t[];
  foreach_up () { double a = 0; foreach_neighbor (1) a += t[]; s[] = a; }
  foreach () s[] = coarse (t);
  foreach_down () foreach_child () t[] = 0;
}
