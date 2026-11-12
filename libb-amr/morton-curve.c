#include "morton-curve.h"
#include "tree.h"

void tree_block_indices (cache (Index) * c)
{
  cache_clear (*c);
  cache_reserve (*c, tree_get_nblocks ());
  foreach_node ()
    cache_push (*c, (Index) {IJKL});
}
