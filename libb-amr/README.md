# block based amr

## Ghost cells

Every block stores its m^d interior cells and `GHOSTS` (default 2) ghost
layers around them, as PARAMESH/FLASH, Athena++ and Parthenon do. A stencil
reads its neighbours in place:

```c
ghosts_update (u);                     /* block-ghosts.h */
BlockShape shape = tree_block_shape ();
foreach_leaf_parallel () {
  const Real * a = block_data (sblocks, u);
  Real * b = block_data (sblocks, v);
  block_cells_simd (shape, 0)
    b [c] = cell_at (a, 1, 0) + cell_at (a, -1, 0) + cell_at (a, 0, 1)
          + cell_at (a, 0, -1) - 4 * cell_at (a, 0, 0);
}
```

`ghosts_update (s, ...)`:
1. **restricts**, finest level first, every parent from its children;
2. **fills the ghosts**, coarsest level first, each block's from
   - the node of its own level beside it (a leaf, or a restricted parent):
     rows copied;
   - otherwise the leaf one level coarser (2:1 balance), by the scalar's
     prolongation. That reads the coarse leaf's interior and first ghost
     layer, which are already valid because that level was done first;
   - outside the domain, the scalar's boundary condition, one axis after
     the other, so that corners come out right.

Within a level the blocks are independent and run on all threads.
`ghosts_update_levels ()` fills parents too, for coarse-level stencils.

**Per-scalar operations** (`ScalarOps`, `scalar_set_ops ()`): a
prolongation, a restriction and a boundary condition. `NULL` takes the
default: `prolong_linear`, `restriction_mean`, `bc_neumann`. Also given:
`prolong_minmod` (no new extrema, for a water depth), `prolong_injection`
and `bc_odd`. One prolongation fills both the children of a refined block
and the ghosts at a level jump, so the two agree. Callbacks are written with
`prolong_cells ()` and `boundary_cells ()`, and may read other scalars of
the same nodes, for coupled operators such as `h` with `zb`.

**Derived scalars** are not prolonged. Fill the ghosts of the scalars they
come from, then compute them on the halo, `block_cells (shape, GHOSTS)`, as
Parthenon's `FillDerived` does.

**Refining and coarsening keep ghosts valid.**
- `tree_refine ()` prolongs the children from the parent (interior and
  ghosts), then fills the children's ghosts.
- `tree_coarsen ()` fills the ghosts of the new leaf, and of the finer
  leaves that touch it, which read its children before.

Every other leaf reads the same values as before. So a balance cascade
never reads a stale halo, as long as the ghosts were valid before. Call
`ghosts_update ()` after changing values; `adapt_wavelet ()` does it for
every scalar first.

**Layout:** compact by default; `-DBAMR_ALIGN=V` starts every interior row
on a multiple of V Reals, for aligned vector loads. Code that needs the
interior coordinates of a cell takes `ci - shape.ox`, `cj - shape.g`, and
`block_at ()` goes the other way.

**Cost** (3D, m = 16, g = 2, 568 leaves with one level jump; serial, -O2,
on a loaded machine, so rough):

| | stored halos | gather (before) |
|---|---|---|
| fill, restriction included | 23–31 µs/leaf | 14–23 µs/leaf (restriction only) |
| 7-point stencil sweep | 12–16 µs/leaf | 16–25 µs/leaf, gather included |
| one fill + one sweep | 35–46 µs/leaf | 30–47 µs/leaf |

The fill costs about what the gather did, and buys linear (not piecewise
constant) ghosts at level jumps. Every further sweep that reads the same
ghosts then runs at stencil speed. The price is memory: 1.95× the interior
in 3D for m = 16 (1.42× for m = 32).

Not done yet: periodic boundaries, and flux correction (refluxing) at
coarse-fine faces, which conservation needs whichever way ghosts are
stored.

## Tests

    sh tests/run.sh

It builds each test in `tests/` with cwrap (`CWRAP`, by default
`../libc-pp/build/cwrap`), in 1, 2 and 3 dimensions, in double and float
(`-DBAMR_SINGLE`), under AddressSanitizer and UBSan, and fails on any warning:
- `test_pool.c`: the hot and cold free lists, with page faults counted;
- `test_tree.c`: counts, the scalar hooks, Morton order, refinement, and the
  loops of `tree.h`;
- `test_blocks.c`: every loop of `block-loops.h` against `block_cells`, in
  both layouts and every halo, and `block_at ()`;
- `test_ghosts.c`:
  - a linear field, extrapolated at the walls: every ghost is exact, at
    level jumps, walls and corners, and on parents with
    `ghosts_update_levels ()`;
  - `bc_neumann` and `bc_odd` mirror over every wall;
  - injection gives the covering coarse cell;
  - minmod is exact for a linear field and stays positive at a kink;
  - a custom restriction runs once per parent;
  - a derived scalar computed on the halo matches filling its ghosts;
  - refining and coarsening keep every leaf's ghosts equal to a fresh
    `ghosts_update ()`;
  - a freed id gets the default operations back.
- `test_adapt.c`: families allocated at once, prolongation (exact for a
  linear field, and keeping the mean), coarsening, 2:1 balance through
  random refining and coarsening, and `adapt_wavelet ()` on a front: it
  settles, keeps the integral, and coarsens a constant field to minlevel.

It also runs `blocks`, `ghosts` and `adapt` with OpenMP (2D, 3D) and with
`-DBAMR_ALIGN=4` (1D to 3D).



# Future plans

## Block allocation

Ideas to make allocating blocks, and sweeping over them, faster. They are
ranked by the gain expected for libb-amr as it is now. The gains are
estimates, not measurements.

How allocation works now:
- Every block comes from a pool (`pool.c`) of 1 MiB `mmap` chunks.
- A large block that is freed gives its pages back to the system at once
  (`madvise`).
- Each node holds a two-level table of its blocks, `s[id/8][id%8]`.

### 1. Stop handing pages back on every free

**Done.** `pool.c` keeps freed objects on a hot list, and
`pool_alloc_uninit ()` skips the zeroing; `pool_trim ()` (or
`grid_scalar_trim ()`) gives pages back after a coarsening. Measured by
`tests/test_pool.c`, reusing 62 KiB blocks: about 26 µs a block before (15
page faults each), 2.4 µs zeroed and 1.3 µs uninitialised now.

Refinement fills children in full, so their blocks use
`pool_alloc_uninit ()`.

**Lazy release, to be decided** (`_BAMR_MADV_DONTNEED_RECOMMENDED_` in
`pool.h`):
- **1, the default:** objects over the hot limit go cold with
  `MADV_DONTNEED`.
- **0:** they go to a "muzzy" list with `MADV_FREE`, all pages but the first
  (which holds the list's link). The system takes those pages only under
  memory pressure, so reuse costs no fault while memory is plentiful.
  `pool_trim ()` still gives everything back for real. The drawback: muzzy
  pages count in the RSS until the system takes them.

Measured by `tests/test_pool.c` (62 KiB blocks, over the hot limit, then
reused uninitialised):

| | faults on reuse | free | alloc |
|---|---|---|---|
| `MADV_DONTNEED` | 15 per block | 8.6 µs | 27.5 µs |
| `MADV_FREE` | 0 | 6.3 µs | 3.9 µs |

That is without memory pressure. Under pressure (a cgroup with a small
`memory.max`), `MADV_FREE` should fault like `MADV_DONTNEED`; this has not
been measured yet.

The reasoning, as planned:

`pool_free` gives a large block's pages back to the system at once. The next
`pool_alloc` of that block then takes a page fault for every 4 KiB page, and
the kernel zeroes each one. A 3D block with m = 16 is 20³ doubles, 62 KiB:
that is 16 faults, very roughly 10–20 µs. A stencil over that block's 4096
cells takes about 4 µs. A grid that refines and coarsens every few steps could
spend more time in page faults than computing.

The fix:
- **A hot free list.** Freed blocks are kept as they are, up to a limit, such
  as the most blocks freed in one regrid. Only blocks beyond the limit go back
  to the system.
- **`pool_alloc_uninit ()`** for blocks that are overwritten straight away, as
  when refinement fills children from their parent. These skip the zeroing
  too.

### 2. Allocate the children of a node as one group

**Done.** `tree_refine ()` allocates a family at once: one object of
TREE_NCHILDREN nodes, and for each scalar one object of TREE_NCHILDREN
blocks, which `tree_coarsen ()` frees together. Children are filled by linear
prolongation, which writes every cell, so their blocks are not zeroed
(`pool_alloc_uninit ()`).

The reasoning, as planned:

`add_children` always creates all 2^d children at once, and in a 2:1-balanced
tree a family is also coarsened at once. So allocate one object of 2^d blocks
per scalar, instead of 2^d separate ones:
- one pool operation instead of 8 in 3D;
- siblings are adjacent in memory, and they are also adjacent in Morton
  order, so a traversal walks memory nearly sequentially;
- stencils that cross between siblings stay within the same pages.

This is the "quadrant family" of p4est and similar codes. `grid-scalar` would
allocate per family, and each child's pointer would be an offset into the
family's object.

### 3. Ghost cells are a large share of memory in 3D

**Tried, then reverted: blocks are stored with their ghosts again** (see
[Ghost cells](#ghost-cells)).

For a while, blocks stored m^d cells, and `block_gather ()` filled a
scratch block with ghosts before every stencil. It halved memory in 3D. It
was reverted because:
- each sweep repeated the gather, and a coarse-neighbour ghost had to be
  rebuilt every time it was read;
- interpolating at level jumps, which second-order stencils need, made that
  repeated work larger;
- custom prolongations and derived scalars fit a single fill per step much
  better.

Storage per block, relative to the interior alone:

| | m = 16, g = 2 | m = 32, g = 2 |
|---|---|---|
| 1D | 1.25× | 1.13× |
| 2D | 1.56× | 1.27× |
| 3D | **1.95×** | 1.42× |

Ways to reduce it, if memory becomes the limit:
- **Larger m**: cheaper ghosts, but coarser refinement.
- **Fewer ghost layers** for scalars whose stencils need only one: a
  per-scalar g would need a block pool per size.
- **Face ghosts only**: a 5- or 7-point stencil never reads edge or corner
  ghosts.

### 4. Huge pages

**Gain: moderate. Effort: small.**

A sweep over many blocks of several scalars touches far more 4 KiB pages than
the TLB can hold: 300 MB of blocks is about 75,000 pages. With 2 MiB pages it
is about 150.
- Make pool chunks 2 MiB-aligned multiples of 2 MiB.
- Call `madvise (MADV_HUGEPAGE)` on them.

Bandwidth-bound codes often gain 5–20%.

This conflicts with 1: giving back part of a huge page splits it. So for
pools of huge pages, rely on the hot free list, and give memory back only in
an explicit `pool_trim ()`.

### 5. Contiguous storage in Morton order

**Gain: large for a GPU and for prefetching. Effort: medium.**

After each regrid:
- number the leaves in Morton order;
- keep each scalar's blocks in one slab, `u[rank * cells + c]`.

A traversal then reads memory sequentially, and the hardware prefetcher
follows blocks across their boundaries. The slab is also the contiguous device
array that `gpu_cells` (`block-loops.h`) wants. It needs no gathering through
`block_table`, and no table of pointers to map.

The cost is copying all blocks at each regrid: one `memcpy` per block, cheap
next to many computation steps. Combined with 2, it becomes "families in
Morton order".

### 6. A flat scalar table per node

**Gain: small, and simpler code. Effort: small.**

`node->s[id/8][id%8]` follows two pointers. It also needs a group pool, and
the group bookkeeping in `grid_scalar_check`.

A flat table of `NSCALAR_MAX` pointers is 512 bytes per node, under 1% of a
single 3D block. It removes one pointer per access and the group pool, and
makes `grid_scalar_node_free` simpler. The groups only pay off when blocks are
tiny, and `tree ()` already requires m ≥ 16.

### 7. NUMA and threads

**Matters on machines with more than one socket.**

Fresh `mmap` pages are placed in the memory of the first thread that writes
them. If blocks are first written in the same `foreach_leaf_parallel` pattern
later used to compute on them, each thread's blocks end up in its own socket's
memory. The uninitialised allocation of 1 helps, because the allocator itself
then never writes the memory.

`pool.c` is not thread-safe. If refinement ever runs in parallel, it needs
per-thread free lists, or a lock around `pool_alloc` and `pool_free`.

### 8. Further options

- **Mixed precision.** Done in part: `real.h` gives `Real`, a double, or a
  float with `-DBAMR_SINGLE`, and blocks store `Real`. Geometry stays double.
  What remains:
  - **Store in float, compute in double.** Blocks hold float, and kernels
    convert each value to double, compute, and convert back. Stencils are
    limited by bandwidth, so this is nearly as fast as float, with double
    arithmetic. It needs a computation type beside `Real`.
  - **Precision per field.** Pressure in double and a tracer in float, in one
    run. It needs a size per scalar, a block pool per size, and a typed
    accessor.
- **Reserving ahead.** Before a regrid the tree knows how many blocks it will
  create, so `pool_reserve (n)` could map them in one call.

### Order of work

| Order | Change | Effort | Expected gain |
|---|---|---|---|
| 1 | hot free list + uninitialised allocation | small | large, if the grid adapts often |
| 2 | flat scalar table | small | small, simpler code |
| 3 | huge pages, aligned chunks | small | moderate |
| 4 | children allocated as a group | medium | moderate (locality, fewer allocations) |
| 5 | Morton-ordered slabs, renumbered at each regrid | medium | large for a GPU and for prefetching |
| 6 | ~~blocks without stored ghosts~~ | large | tried and reverted: see 3 |

**Measure before changing anything:** time a stencil sweep, and a refine and
coarsen cycle. Use `grid_scalar_check ()` and `pool_used ()` to make sure no
change leaks blocks.

# Comparison of AMR implementations

How to store the tree itself: its nodes, and the way to find a node or its
neighbours. This is about the tree's structure, not the blocks of field data.

## The starting idea: dense tiles of N levels

One idea is to store the tree densely, N levels at a time:
- a complete tree of depth N is laid out contiguously, as a dense array;
- the same is done again for each subtree whose root is at a depth that is a
  multiple of N.

OpenVDB, NanoVDB and other sparse-voxel trees are built this way.

**The risk is sparseness.** In 3D, a tile of N = 4 levels has
1 + 8 + 64 + 512 + 4096 = 4681 slots. A tile holding just one refined corner
path uses about 33 of them, under 1%.

## Keep the tree's structure separate from its data

A node's structure is a few bytes: flags, an index, a pointer to its scalar
table. A block of field data is about 62 KiB per scalar in 3D with m = 16.
So:
- the structure is roughly 0.01% of memory;
- walking it costs a few nanoseconds per block, against about 4 µs of stencil
  work on the block's 4096 cells.

What follows from that:
- **Sparse structure is cheap.** 37 KiB of mostly empty slots in a tile costs
  about what one block costs.
- **Sparse data storage would be ruinous.** Laying out blocks densely in
  tiles would reserve 62 KiB per empty slot, per scalar.
- **Contiguity matters most for the blocks.** Items 2 and 5 of
  [Block allocation](#block-allocation), sibling groups and Morton-ordered
  slabs, address that. The structure can be organised for fast lookups
  instead.

## Alternatives for the tree's structure

### A. Per-level sparse arrays (what libb-amr has now)

`memindex` gives each level a virtual `[i][j][k]` array, and the operating
system only commits the pages that are actually written.

- **Good:** finding a node or a neighbour is pure index arithmetic, O(1),
  with no pointer chasing. Basilisk uses this scheme.
- **Cost:** sparseness becomes virtual address space, plus one 4 KiB page per
  isolated row of nodes, as a page holds 512 pointers. Refined regions are
  usually clustered, so pages fill well in practice.
- **Verdict:** this is already a version of the dense-tile idea in which the
  operating system handles the sparseness. Measure before replacing it.

### B. Dense tiles, with occupancy bitmasks and rank indexing (NanoVDB-style)

Each tile stores:
- a bitmask of its occupied slots: 8 bits per family in 3D, or 4096 bits
  (512 bytes) for a whole tile's bottom level;
- an array of only the occupied nodes, in Morton order.

To find slot `s`, count the set bits of the mask below `s` (`popcount`). With
a small table of running counts, one per 64-bit word, that is O(1). This
removes the sparseness completely, while keeping the tiles contiguous.

libb-amr allows a simplification. In a 2:1 tree whose families are always
created whole, a family's mask is all ones or all zeros. The per-family mask
is then a single bit, and a tile is just a record of which families exist.

- **Good:** contiguous and compact, with O(1) lookup. A tile copies to a GPU
  in one piece.
- **Cost:** refining inserts into a tile's compact array, so tiles should be
  rebuilt in a batch at each regrid, rather than edited one node at a time.

### C. Linear octree: a sorted array of Morton keys (p4est, Dendro)

Store only the leaves, as keys made of the level and the interleaved
coordinates, in sorted order.

- **Good:** fully contiguous, and never sparse. A traversal is a scan.
  Splitting the work across MPI ranks means cutting the array into pieces. It
  goes to a GPU as it is.
- **Cost:** finding a neighbour needs a binary search, O(log n), or a hash on
  the side. Each regrid rebuilds the array, an O(n) merge, which is fine when
  regrids are batched.
- **Verdict:** the standard choice once MPI is needed.

### D. Hashed octree (Warren–Salmon): Morton key → node, in a hash table

- **Good:** O(1) lookup at any level. Memory is proportional to the number of
  nodes, and never sparse.
- **Cost:** a traversal in hash order jumps around memory. It needs a sorted
  list of the leaves beside it for traversals.

### E. Child families plus compaction, keeping a pointer tree

Each node points to one group of 2^d children, allocated together. Whole
families are always created, so nothing is wasted. After each regrid, the
families are copied into Morton order.

- **Good:** the smallest change from the code as it is.
- **Cost:** finding a neighbour goes back to walking up and down the tree, so
  the arithmetic lookups of A are lost.

## Comparison

| | Lookup / neighbour | Traversal locality | Wasted space | Refining | GPU / MPI |
|---|---|---|---|---|---|
| A. per-level arrays (now) | O(1) arithmetic | good, in clustered regions | virtual, per page | O(1) per node | awkward |
| dense tiles | O(1) | excellent | high when refinement is sparse | O(1) | good |
| B. tiles + bitmask/rank | O(1) with `popcount` | excellent | none | batched rebuild | very good |
| C. linear octree | O(log n), or a hash | excellent | none | batched O(n) | best |
| D. hashed | O(1) expected | poor alone | none | O(1) | fair |
| E. families + compaction | walk the tree | good after compaction | none | O(1) | fair |

## Recommendation

1. **For now, keep A for the structure.** Spend the effort on the data
   instead: sibling-group allocation and Morton-ordered slabs for blocks.
   That is where the memory traffic is.
2. **When a GPU or MPI is needed, add C beside A.** C is a sorted array of
   leaf keys, rebuilt at each regrid. A stays the fast way to find
   neighbours, and C becomes the contiguous list for GPU launches and for
   splitting the work across ranks. Each regrid already traverses the whole
   tree, so keeping the two in step is cheap.
3. **For tiles, build them as B,** with bitmasks and rank indexing. Pick N by
   dimension: N = 2 in 3D (a tile of 73 slots), and N = 3 or 4 in 2D. Tiles
   suit regions that are refined fairly evenly, and the bitmask removes the
   cost of those that are not.

**To decide, measure:** prototype C beside the current tree, and compare the
two on a benchmark that refines and traverses.
