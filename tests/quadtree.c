#include <stdio.h>

#include "libb-amr/tree.h"
#include "libb-amr/morton-curve.h"

int main ()
{
  tree (3, 1<<6);
  morton_start (0, 0, 0)
  {
    printf("\n[%d %d %d]", i, j, level);
    fflush (stdout);
  } morton_end ();
  tree_destroy ();  
}
