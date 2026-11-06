#include "morton-curve.h"
#include "tree.h"

void tree_block_indices (cache (Index) * c)
{
  cache_clear (*c);
  cache_reserve (*c, tree_get_nblocks ());
  morton (TREE_ROOT_INDEX)
    cache_push (*c, (Index) {IJKL});
}
