#include "stencil.h"
scalar s, t;
void f (void)
{
  foreach_neighbor (1) s[] = 0;
  foreach () foreach_neighbor (1) foreach_neighbor (1) s[] += t[];
  foreach () foreach_neighbor (1) s[] = t[];
}
