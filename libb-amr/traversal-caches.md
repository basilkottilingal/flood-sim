# Traversals from caches, not from morton (): a plan

*Deferred (2026-11-11). Nothing here is implemented yet.*

The aim: loops over the tree read a list of nodes kept in a cwrap `cache`,
made again only when the tree changes. `morton ()` walks the tree only to
fill those caches.

## Now

`morton (TREE_ROOT_INDEX)` walks the whole tree on every use (with a stack
of `TREE_NCHILDREN * BAMR_MAX_DEPTH` entries):

| Where | What, at each call |
|---|---|
| `foreach_leaf ()` (`block-loops.h`) | a walk, to visit the leaves |
| `foreach_leaf_parallel ()` | a walk, gathering the leaves into a new `cache (Index)`, then the loop |
| `levels_new ()` (`block-ghosts-impl.h`) | two walks, counting then filling the nodes by level, at each `ghosts_update ()` |
| `tree_restrict ()` | a walk |
| `adapt.c`, `grid-scalar.c` | walks, to flag, to allocate and free blocks |
| `tree_block_indices ()` (`morton-curve.c`) | a walk, into a `cache (Index)`: the one that already caches |

So every `foreach ()` (built on `foreach_leaf`) and every ghost update pays
one or more walks, whose cost grows with the tree, while the tree changes
only at a refinement or a coarsening.

## Proposed

- **Caches of indices, kept with the tree:**
  - the leaves, in Morton order: `cache (Index) tree_leaves`;
  - the nodes by level, leaves only and parents too (what `Levels` holds
    now): one cache of all, with the first of each level;
  - later, for the plans of `multigrid.md`: the level lists that
    `foreach_level (l)`, `foreach_level_or_leaf (l)` and
    `foreach_coarse_level (l)` read.
- **Made by `morton ()`, through the cache's own update:** cwrap's caches
  have `update` and `is_dirty` (`include/cwrap/cache.h`). The update is the
  walk (as `tree_block_indices ()` is), and a cache is made again only when
  read while dirty.
- **Marked dirty when the tree changes:** `tree_refine ()`, `tree_coarsen ()`,
  `tree ()` and `tree_destroy ()` mark every tree cache dirty. Nothing else
  changes the tree, so nothing else needs to.
- **The loops read them:**
  - `foreach_leaf ()` becomes `foreach (p in tree_leaves)`, with `IJKL`,
    `node` and `sblocks` from `p`, as now;
  - `foreach_leaf_parallel ()` shares out the same cache among threads, with
    no gathering of its own;
  - `ghosts_update_list ()` reads the levels cache instead of
    `levels_new ()`.
- **`morton ()` stays,** as the one walk, for the caches' updates, and for
  code that changes the tree as it walks (adapt.c's flagging, `grid-scalar.c`).

## What to keep in mind

- **The order is the same:** a cache filled by `morton ()` has the nodes in
  Morton order, so every loop visits them as now, and results stay bit for
  bit the same: the check to make.
- **A loop that changes the tree** (refines inside `foreach_leaf ()`) would
  read a cache made dirty under it. Either it is an error (checked: the
  cache dirty during its own loop), or such code keeps `morton ()`.
- **`node` from an index** is a lookup in the level index (`level-index.h`)
  per node, where the walk had the pointer at hand. Keeping the node
  pointers in the cache with the indices (a `cache` of `{ Index, BAMRNode * }`)
  avoids it; a coarsening frees nodes, so these too are made again on a
  change, which the dirty flag already covers.
- **Memory:** one `Index` (and pointer) per node, a few tens of bytes per
  block of thousands of cells.
- **Threads:** the update must not run in two threads at once: the caches
  are made on the first read after a change, outside parallel regions, or
  under a lock.

## Order

1. The leaves cache, and `foreach_leaf ()` and `foreach_leaf_parallel ()`
   over it; marking dirty in `tree_refine ()`, `tree_coarsen ()`. Checked bit
   for bit (the checksum program of the operator work), and timed: a
   `foreach ()` with no walk.
2. The levels cache, for `ghosts_update_list ()` and `tree_restrict ()`.
3. The level lists of `multigrid.md`, when that comes.
