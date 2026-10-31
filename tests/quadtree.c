#include <stdio.h>

#include "libb-amr/tree.h"
#include "libb-amr/morton-curve.h"
#include "libb-amr/scalar.h"

int main ()
{
  tree (3, 6);
  morton_start (0, 0, 0)
  {
    printf("\n[%d %d %d]", i, j, level);
    fflush (stdout);
  } morton_end ();

  scalar s = scalar_new ();
  vector v = vector_new ();
  printf ("\nscalar s %u. vector v %u %u", s, v.x, v.y);
  scalar_free (s);
  vector u = vector_new ();
  printf ("\nvector u %u %u", u.x, u.y);

  tree_destroy ();  
}
