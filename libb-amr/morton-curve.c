#include "morton-curve.h"
#include "tree.h"

Array * tree_block_indices ()
{
  Array * a = array_new ();
  int i = 0, j = 0, level = 0;
  morton_start (i,j,level)
    array_append (a, & ((Index) {i,j,level}), sizeof (Index));
  morton_end ();
  return a;
}
