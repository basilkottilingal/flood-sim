#include <stdio.h>
#include "stencil.h"
scalar s, t;
int main (void)
{
  foreach () {
    double sum = 0;
    foreach_neighbor (2)              /* the mock has one ghost layer */
      sum += t[];
    s[] = sum;
  }
  printf ("not reached\n");
  return 0;
}
