#ifndef _QUADTREE_MORTON_CURVE_H_
#define _QUADTREE_MORTON_CURVE_H_

  #include <cwrap/cache.h>
  #include "tree.h"

  /*
  .. morton (p) { ... }
  .. runs the block once for every node of the subtree at the Index p
  .. (TREE_ROOT_INDEX for the whole tree), in any dimension. A cwrap @loop: the
  .. block sees IJKL, node and sblocks of the node visited, and nothing else of
  .. the traversal.
  ..
  .. Children come in Morton (Z) order, index_child (p, c) for c = 0, 1, ...,
  .. i varying fastest. A leaf is visited when it is reached, a parent after
  .. all of its children. 'break' and 'continue' go to the next node.
  ..
  .. Each step down pushes the parent and all but the first child, so the
  .. stack never holds more than TREE_NCHILDREN per level.
  */
  #define TREE_ROOT_INDEX  ((Index) {TREE_ROOT})

  @loop morton (Index root) : IJKL, node, sblocks
  {
    struct mstack { Index p; int visited; }
      stack [TREE_NCHILDREN * BAMR_MAX_DEPTH + 1];
    int depth = 0, visited = 0;
    Index cur = root;
    if (index_is_allocated (cur))
      for (;;)
      {
        while (!visited)
        {
          BAMRNode * n = index_node (cur);
          if (!tree_has_children (n))
            break;
          stack [depth++] = (struct mstack) {cur, 1};
          for (int c = TREE_NCHILDREN - 1; c > 0; c--)
            stack [depth++] = (struct mstack) {index_child (cur, c), 0};
          cur = index_child (cur, 0);
        }
        do {
          INDEX_UNPACK (cur);
          BAMR_NODE_ENV ();
          {...}
        } while (0);
        if (!depth)
          break;
        depth--;
        cur     = stack [depth].p;
        visited = stack [depth].visited;
      }
  }

  /*
  .. API
  .. tree_block_indices (&c) puts the Index of every node in c, in the order
  .. of morton (), after emptying it. c is the caller's: a cache cannot be
  .. returned, as the copy would share its array. A local cache is freed as
  .. its scope ends:
  ..
  ..   cache (Index) blocks;
  ..   tree_block_indices (&blocks);
  ..   foreach (b in blocks)
  ..     ... b.level ...
  */
  void tree_block_indices (cache (Index) * c);

#endif
