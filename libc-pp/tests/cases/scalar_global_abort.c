#include <stdio.h>
scalar zb;
int main (void)
{
  scalar mine = scalar_new (), * p = &zb;
  scalar_free (mine);               /* a scalar of its own: freed */
  scalar_free (*p);                 /* a global, through a pointer: stops */
  printf ("not reached\n");
  return 0;
}
