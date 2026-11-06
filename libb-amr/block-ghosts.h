#ifndef _BAMR_BLOCK_GHOSTS_H_
#define _BAMR_BLOCK_GHOSTS_H_

  /*
  .. Ghost cells, for blocks stored without them.
  ..
  .. A block stores m^dimension cells (tree_block_m ()). A stencil that reads
  .. neighbours runs on a scratch block with g ghost layers, of
  .. block_shape (m, g) or any layout block-loops.h knows:
  ..
  ..   block_gather (w, shape, p, s)   fills w: the interior from the block of
  ..                                   scalar s at p, the ghosts from its
  ..                                   neighbours
  ..   ... the stencil, reading w, writing a stored block, by block_interior ()
  ..   block_scatter (p, s, w, shape)  if the result is in w: writes w's
  ..                                   interior back to the block
  ..
  .. A ghost cell's value:
  ..
  ..   - outside the domain, the cell reflected back inside: a zero gradient
  ..     boundary;
  ..   - else the cell of the node of the same level that holds it, leaf or
  ..     not: a parent's block holds the average of its children, once
  ..     tree_restrict () has run;
  ..   - else, where the neighbour is coarser, the value of the coarser cell
  ..     that covers it: piecewise constant. The root covers the domain, so
  ..     there always is one.
  ..
  .. g must be at most m. Not thread safe while blocks change; many threads
  .. may gather at once, each in a scratch block of its own.
  ..
  .. Why no ghosts in storage: in 3D with m = 16 and g = 2, a block with
  .. ghosts is 20^3 cells, 1.95 times the 16^3 it computes. Stored without
  .. them, memory halves, and the scratch block, one per thread, stays in
  .. cache. The gather does the work a ghost update would do anyway.
  */

  #include "block-loops.h"

  /* a stored block: no ghosts */
  static inline BlockShape tree_block_shape (void)
  {
    return block_shape (tree_block_m (), 0);
  }

  /* a scratch block with g ghost layers, compact */
  static inline BlockShape tree_scratch_shape (int g)
  {
    return block_shape (tree_block_m (), g);
  }

  /* api */
  void block_gather  (Real * w, BlockShape shape, Index p, int s);
  void block_scatter (Index p, int s, const Real * w, BlockShape shape);
  void tree_restrict (int s);   /* each parent: the mean of its children */

  /*
  .. foreach_leaf_gathered (s, g) { ... }: every leaf, with w, a scratch
  .. block of scalar s with g ghost layers already gathered, and its shape,
  .. besides IJKL, node and sblocks. w is the same memory for every leaf,
  .. freed at the end however the loop ends.
  ..
  ..   foreach_leaf_gathered (u, 1) {
  ..     Real * out = block_data (sblocks, v);
  ..     block_interior (shape)
  ..       out [o] = cell_at (w, 1, 0) + cell_at (w, -1, 0) - 2 * cell_at (w, 0, 0);
  ..   }
  */
  static inline void block_ghosts_free (Real ** w)
  {
    free (*w);
  }

  @loop foreach_leaf_gathered (int scalar_id, int ghosts) : IJKL, node, sblocks, w, shape
  {
    BlockShape shape = tree_scratch_shape (ghosts);
    Real * w __attribute__ ((cleanup (block_ghosts_free))) =
      malloc ((size_t) shape.cells * sizeof (Real));
    foreach_leaf () {
      block_gather (w, shape, (Index) {IJKL}, scalar_id);
      {...}
    }
  }

#endif
