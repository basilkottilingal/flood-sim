#include <stdio.h>

#include "../libb-amr/morton-curve.h"

int main ()
{
  tree (3, 6);
  morton ((Index) {0, 0, 0})
  do {
    printf("\n[%d %d %d]", i, j, level);
    fflush (stdout);
  } while (0);

  scalar s;
  vector v;
  printf ("\nscalar s %u", s);
  printf ("\nvector v %u %u", v.x, v.y);
  vector u;
  printf ("\nvector u %u %u", u.x, u.y);

  tree_destroy ();  
}
