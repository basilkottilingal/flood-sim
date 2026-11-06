#include "morton-curve.h"
#include "tree.h"

Array * tree_block_indices (void)
{
  Array * a = array_new ();
  morton (TREE_ROOT_INDEX)
    array_append (a, & ((Index) {IJKL}), sizeof (Index));
  return a;
}
