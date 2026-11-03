#include <stdio.h>
#include "scalar.h"
scalar m1, m2;
extern scalar o1, o2, o3;
extern int n_global_scalars;
scalar hidden_value (void);
int main (void)
{
  scalar all[] = { m1, m2, o1, o2, o3, hidden_value () };
  unsigned seen = 0;

  for (int i = 0; i < 6; i++)
    seen |= 1u << all[i];
  printf ("%d %x\n", n_global_scalars, seen);
  return 0;
}
