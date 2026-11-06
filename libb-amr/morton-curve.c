#include "morton-curve.h"
#include "tree.h"

Array * tree_block_indices (void)
{
  Array * a = array_new ();
  morton_start (TREE_ROOT)
    array_append (a, & ((Index) {IJKL}), sizeof (Index));
  morton_end ();
  return a;
}
