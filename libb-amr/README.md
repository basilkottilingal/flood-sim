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

**The same stencil with cwrap's `foreach ()`** (see `../libc-pp/README.md`,
"Stencils"):

```c
ghosts_update (u);
foreach ()
  v[] = u[1,0] + u[-1,0] + u[0,1] + u[0,-1] - 4 * u[];
```

cwrap expands it into `foreach_leaf ()` and `foreach_cell (shape, 0)`,
with one `block_data ()` per leaf for each scalar used (`const` for `u`, only
read) and `u[i,j]` as `cell_at`. It is serial over the leaves. In the body,
`x`, `y`, `z` and `Delta` are the cell's centre and size, from
`tree_get_origin (axis)` and `tree_get_size ()`, the domain that
`tree_shift_resize ()` sets.

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

- `test_stencil.c`: cwrap's `foreach ()` and `s[i,j,k]`: every interior cell
  once and no ghost written, stencils in every dimension equal to `cell_at`
  bit for bit, offsets that are expressions, `+=` and `++` on a cell, a
  vector's components, an array of scalars and a pointer to one, a scalar
  parameter, and names in the body (`c`, `level`, `i`, `node`, `sblocks`)
  that the loops' own names do not hide; `x`, `y`, `z` and `Delta` on a
  moved, resized domain with two levels of leaves, against the centre of
  each cell computed apart. Also with OpenMP and
  `-DBAMR_ALIGN=4`.

- `test_index.c`: the three level indices side by side (memindex,
  voxel-memindex, nano-vdb, dense and compact) against a reference hash
  map: random assigns, overwrites and frees, clustered and scattered; every
  entry and the indices around it, by get and by cursor, out of range too;
  emptying (nano-vdb gives back every node), destroying and reusing; levels
  on both sides of nano-vdb's hashed root (above 4096).

It also runs `blocks`, `ghosts` and `adapt` with OpenMP (2D, 3D) and with
`-DBAMR_ALIGN=4` (1D to 3D), and `tree`, `blocks`, `ghosts` and `adapt` on
the other level indices: `-DBAMR_INDEX=2` (voxel-memindex), `3` (nano-vdb)
and `3 -DNVDB_COMPACT` (2D, 3D).

### Benchmark: the level indices

    sh tests/bench_memindex.sh [RAW]

It builds `tests/bench_memindex.c` with gcc, in 2D and 3D, and runs each case
(dense, shell, river, random) once per implementation (memindex,
voxel-memindex, nano-vdb, nano-vdb compact), each in its own process. It
prints a table, and checks that every implementation returns the same values
as memindex (checksum). nano-vdb's dense nodes skip the scattered 3D cases:
about 300 KB per isolated point. Results and their reading are in
`tests/bench_memindex.md`, summed up below in
[Measured: memindex against voxel-memindex](#measured-memindex-against-voxel-memindex).
The sparse cases need about 1.5 GB for memindex. Cache and dTLB misses are
printed only if `kernel.perf_event_paranoid` is 2 or less.

    sh tests/bench_balanced.sh [RAW]

It does the same on 2:1 balanced trees, built as `tree.c` builds them. It
times a depth-first Morton traversal, a pass over a level list, and the
neighbour box of every leaf. Results are in `tests/bench_balanced.md`.


## Level index: memindex, voxel-memindex, nano-vdb

Each level of the tree has an index: a sparse (2^level)^d array of node
pointers. `level-index.h` picks the implementation when compiling, and
`tree.c` and `tree.h` use only its interface:

| `-DBAMR_INDEX=` | files | structure |
|---|---|---|
| `1` (default) | `memindex.[ch]` | per-level virtual arrays, as Basilisk: `d` loads, the system commits pages as they are written |
| `2` | `voxel-memindex.[ch]` | a tree of 64-way nodes (4^3, 8^2, 2^6), 512 bytes each, from a pool |
| `3` | `nano-vdb.[ch]` | a VDB tree: root hash, upper 32^d, lower 16^d, leaf 8^d nodes, each with a mask of its occupied slots; `-DNVDB_COMPACT` stores only the occupied children, found by rank (`popcount`) |

The interface, for a `LevelIndex * m` (IJK: `i`, `j`, `k` as the dimension
has them; `len`: `1 << level`):

| | |
|---|---|
| `lidx_get (m, IJK)` | the pointer, or `NULL`, out of range too: one walk |
| `lidx_allocated (m, IJK)` | `lidx_get () != NULL` |
| `lidx_assign (m, IJK, len, p)` | `p`, not `NULL`, at IJK |
| `lidx_free (m, IJK, len, p)` | IJK empty again; empty nodes go back at once |
| `lidx_destroy (m)` | everything freed; `m` empty and reusable |
| `lidx_get_cursor (m, &c, IJK)` | `lidx_get ()` through a `LevelCursor c` |

A cursor makes nearby lookups cheaper, and is reset (`{0}`) after any
`lidx_free ()` or `lidx_destroy ()` of its index:
- memindex: nothing to keep, `mem_get ()`;
- voxel-memindex: the last bottom node (an index in the same 64 is one
  load);
- nano-vdb: the last leaf, lower and upper node, as VDB's `ValueAccessor`
  (the same leaf is one load, the same lower node two).

Each implementation also keeps its own names (`mem_`, `vmem_`, `nvdb_`),
so that the three can be used side by side, as `tests/test_index.c` and the
benchmarks do. `BAMR_INDEX_NAME` names the one chosen.

**nano-vdb, in short.** An index's bits are cut as VDB cuts a coordinate:
bits 2..0 pick a leaf slot, 6..3 a lower-node slot, 11..7 an upper-node slot,
the first index most significant. While `len <= 4096` there is a single
upper node; above, the root is a hash table of upper nodes keyed by the
indices `>> 12` (open addressing, deletion by backward shift). Dense nodes
are tables of every slot, as OpenVDB's internal nodes: one load per level,
but 32^3 slots (256 KiB) for a 3D upper node. Compact nodes hold a mask, a
running count of set bits per 64-bit word, and only the occupied children;
slot s is at `prefix [s/64] + popcount (mask [s/64] & ((1 << s%64) - 1))`.
Build compact nodes with `-mpopcnt` or `-march=native`: otherwise gcc calls
a library function for each popcount.

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

The reverse holds for `memindex`: huge pages only hurt its sparse arrays. With
transparent huge pages set to `[always]`, one pointer written can commit 2 MiB
instead of 4 KiB, up to 512 times the memory. `mymap ()` therefore calls
`madvise (MADV_NOHUGEPAGE)` on every array it maps. Under `[always]`,
khugepaged can also collapse a pool chunk whose cold objects were given back.
If that shows up, give pool chunks `MADV_NOHUGEPAGE` too, or adopt the
huge-page pools described above.

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

## Tree structure: voxel-memindex

`voxel-memindex.[ch]` (F, compared under
[Measured: memindex against voxel-memindex](#measured-memindex-against-voxel-memindex))
beat memindex on memory and build time for every level that is not full. Next
steps, in order:

1. **Wire it into `tree.c` behind a switch.** Done, with nano-vdb too:
   `-DBAMR_INDEX=2` or `3` (see
   [Level index](#level-index-memindex-voxel-memindex-nano-vdb)).
   `BAMR_NODE_ENV ()` makes a single lookup. The tests pass on every index.
   What remains: time `adapt` and `ghosts_update ()` on a sparse 3D case
   with each index. That is the measure that decides, rather than isolated
   lookups.
2. **Use cursors in the hot loops.** `block-ghosts` gathers neighbours and
   `adapt` checks 2:1 balance: both visit nearby indices, so one cursor per
   level makes most lookups a single load. Cursors must be reset after a
   `vmem_free ()`, so reset them at the end of each regrid pass. Each OpenMP
   thread needs its own cursors.
   - For the neighbour box, add a box cursor that holds the up to 2^d bottom
     nodes a leaf's box can touch. On balanced trees the box is the one place
     where memindex is still faster, by up to 36% (see On 2:1 balanced
     trees).
3. **Compact nodes for sparse levels.** A node holds a 64-bit mask and only
   its occupied children, found by rank (`popcount` from libutils/bits.h).
   - It would cut the 2–4 kB per isolated point by up to 8 times.
   - Children sit in arrays of varying sizes, from a few size-classed pools.
   - Inserting moves the children of one node, so batch it per regrid.
   - Possibly mixed: full nodes near the bottom, compact ones above.
4. **Memory after freeing.** Freed nodes stay hot in the pool. Call
   `pool_trim ()` after a coarsening pass, or give voxel its own hot limit
   (`pool_set_hot_max ()`), so that memory returns after a large coarsening.
5. **Measure the misses.** Rerun the benchmark with
   `kernel.perf_event_paranoid=1`, to get LLC and dTLB misses per lookup.
   Rerun it on a cluster node with transparent huge pages `[always]`, and with
   memindex without `MADV_NOHUGEPAGE`, to quantify the blow-up seen elsewhere.
6. **Later, for MPI or GPU:** the nodes already give Morton order, so a
   linear octree (C) can be built in one traversal of the leaves.

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
  usually clustered, so pages fill well in practice. Each touched page also
  needs its share of a page-table page, which stays after the page is given
  back. It depends on the system: page size, overcommit, and transparent huge
  pages (opted out with `MADV_NOHUGEPAGE`, see Huge pages).
- **Verdict:** this is already a version of the dense-tile idea in which the
  operating system handles the sparseness. Measured against a tree of 64-way
  nodes below: it is cheap only on levels refined fairly evenly.

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
| F. 64-way node tree (`voxel-memindex`) | ceil(L/B) loads; one with a cursor | good | per node of 64 | O(depth) | fair |
| G. VDB tree (`nano-vdb`) | 4 loads (root hash, upper, lower, leaf); one with the accessor | good | dense: 32^d-slot tables; compact: masks only | O(1) per level; compact: a shift in a node | good (fixed depth) |

## Measured: memindex against voxel-memindex

`voxel-memindex.[ch]` is a prototype of the packed-node idea (F above). It is
not wired into `tree.c` or the `Makefile`.
- Each node is 64 pointers, 512 bytes, from a pool. It takes B bits of every
  index: B = 2 in 3D (two tree levels), 3 in 2D, 6 in 1D. A level-L array is
  ceil (L / B) nodes deep.
- Its interface is the memindex one, prefixed `vmem_`, plus:
  - `vmem_get ()`, which finds and reads an entry in one walk, where
    memindex needs `mem_allocated ()` and then `mem_data ()`;
  - `vmem_get_cursor ()`, which remembers the last bottom node, so that an
    index in the same group of 64 costs a single load.
- Empty nodes go back to the pool as soon as their last entry is freed.

`tests/bench_memindex.sh` compares them (see Tests). The numbers below come
from an i7-8650U laptop with about 2.3 GB free, Linux 6.8, transparent huge
pages set to `[madvise]`, and memindex with `MADV_NOHUGEPAGE`, on a quiet
machine (2026-11-08). Each case ran once. `tests/bench_memindex.md` now holds
the later run with nano-vdb too (see
[With nano-vdb](#with-nano-vdb)).

| case (points) | B/point mem / voxel | build ns/point mem / voxel | shuffled lookup ns mem / voxel | neighbours ns mem / voxel+cursor | absent entry ns mem / voxel |
|---|--:|--:|--:|--:|--:|
| 2D dense L11 (4.2M) | 8 / 8 | 9 / 13 | 16 / 56 | 4.2 / 4.8 | — |
| 2D river L16 (2.1M) | 168 / 19 | 184 / 23 | 62 / 105 | 6.8 / 4.8 | 934 / 8 |
| 2D random L20 (50k) | 8370 / 1907 | 6398 / 1351 | 118 / 268 | 44 / 67 | 103 / 49 |
| 3D dense L7 (2.1M) | 32 / 8 | 48 / 15 | 32 / 69 | 6.1 / 6.3 | — |
| 3D shell L9 (0.5M) | 1027 / 47 | 1915 / 33 | 54 / 72 | 13.6 / 9.8 | 42 / 19 |
| 3D shell L10 (2.1M) | out of memory / 47 | — / 33 | — / 114 | — / 10.2 | — / 17 |
| 3D river L12 (0.85M) | 344 / 26 | 594 / 34 | 86 / 137 | 8.7 / 12.2 | 26 / 15 |
| 3D random L20 (30k) | 16715 / 3636 | 12325 / 2974 | 291 / 629 | 76 / 172 | 72 / 37 |

What it shows:
- **Memory.** voxel takes 3–60 times less wherever a level is not full, and
  about the same, or less in 3D, on a full level. memindex pays a page, and
  part of a page-table page, per row touched. On the 3D shell at level 10 that
  came to more than the 2.3 GB free, where voxel needed 99 MB.
- **Building.** voxel is 5–60 times faster on sparse levels, because it takes
  far fewer page faults. Both are about the same on full levels.
- **Finding an entry that exists, in random order.** memindex is 1.3–3 times
  faster, more so on deeper levels. voxel makes ceil (L / B) dependent loads,
  where memindex makes one per dimension. On a full level in Morton order,
  memindex is about 3 times faster (2.5 against 7.4 ns in 2D).
- **Neighbours.** With the cursor, voxel is about as fast as memindex on full
  and clustered levels, and faster on the 3D shell and the 2D river. It is
  1.5–2.3 times slower on isolated points (random, deep levels).
- **Absent entries** (indices drawn uniformly over a sparse level). voxel is
  2–120 times faster. memindex reads untouched pages, and each read is a page
  fault on the zero page.
- **After freeing everything.** memindex keeps its page tables and its upper
  rows (about 570–590 MB on the random level-20 cases). voxel keeps its freed
  nodes, hot, in the pool, until `pool_trim ()` or `vmem_destroy ()`.
- **voxel's weak case:** isolated points on deep levels, at 2–4 kB per point.
  Each point owns most of the nodes on its path.

In libb-amr the tree indexes blocks, not cells. Stencils run inside padded
blocks, and `block-ghosts` gathers each block's neighbourhood once. So the
slower lookups of voxel are mostly off the hot path, while its memory and
build costs are on it, at every regrid. Confirming that needs voxel wired into
`tree.c` (see Future plans).

### On 2:1 balanced trees

The random cases above have isolated points, which a 2:1 balanced tree never
has. In such a tree, families are whole and lie inside one bottom node, and
every refined cell has neighbours around it. `tests/bench_balanced.sh` builds
balanced trees as `tree.c` does, each level in its own array, and times what
libb-amr actually does with them. The table below is from 2026-11-09, on a
quiet machine; `tests/bench_balanced.md` now holds the later run with
nano-vdb too.

| case (nodes) | B/node mem / voxel | Morton dfs ns/node: mem / voxel / voxel+cursor | neighbour box ns/leaf: mem / voxel / voxel+cursor |
|---|--:|--:|--:|
| 2D uniform, depth 10 (1.4M) | 9 / 8 | 5.5 / 8.1 / 4.7 | 38 / 61 / 52 |
| 2D circle, depth 15 (0.64M) | 615 / 25 | 20.4 / 13.8 / 7.1 | 106 / 101 / 77 |
| 2D river, depth 15 (1.3M) | 234 / 23 | 11.7 / 12.5 / 6.6 | 71 / 92 / 73 |
| 3D uniform, depth 6 (0.3M) | 77 / 8 | 32.9 / 25.2 / 7.3 | 272 / 310 / 274 |
| 3D sphere, depth 9 (2.3M) | 252 / 14 | 44.5 / 28.3 / 8.5 | 310 / 462 / 392 |
| 3D river, depth 10 (0.37M) | 446 / 18 | 38.9 / 28.7 / 8.6 | 321 / 472 / 393 |

- **Morton traversal.** With a cursor per level, voxel is faster in every
  case: 1.2x on the uniform 2D tree, 1.8-2.9x on the 2D fronts, and 4.5-5.3x
  in 3D. A family's children are in one voxel node, while memindex reads them
  from 2 (2D) or 4 (3D) rows. Over a level list, as `ghosts_update ()` walks,
  voxel with a cursor is 1.2-3x faster.
- **Neighbour box** (3^d - 1 neighbours per leaf, about a fifth of them
  falling back to a coarser level). memindex keeps an edge here. With
  cursors, voxel ranges from 28% faster (2D circle) to 36% slower (2D
  uniform), and is up to 26% slower in 3D. The worst gap is about 80 ns per
  leaf. That is an estimated percent or less of the ghost copies that
  follow, for blocks of 16^3.
- **Memory and build** are as above: 10-25x less memory, and 12-25x faster
  builds, on every non-uniform tree.

So in a Morton-ordered pass over a balanced tree, voxel's slower lookups do
not show; they show only in the neighbour box, and there by less than the
work done per block. A box cursor would narrow that gap: it would hold the
up to 2^d bottom nodes that a leaf's box can touch, fetched once per leaf.

### With nano-vdb

Both benchmarks were run again (2026-11-10) with nano-vdb, dense and compact
(`-DNVDB_COMPACT`), all built with `-mpopcnt`. A browser was busy on the
machine, so the times are slower than above; compare implementations within a
case. Full tables and readings: `tests/bench_memindex.md`,
`tests/bench_balanced.md`.

| balanced tree (nodes) | B/node mem / voxel / nvdb / nvdb-compact | dfs+cursor ns/node mem / voxel / nvdb / compact | neighbour box+cursor ns/leaf mem / voxel / nvdb / compact |
|---|--:|--:|--:|
| 2D circle, depth 15 (0.64M) | 615 / 25 / 30 / 13 | 54 / 11.6 / 9.8 / 11.4 | 228 / 135 / 103 / 107 |
| 2D river, depth 15 (1.3M) | 234 / 23 / 28 / 13 | 46 / 10.3 / 12.5 / 10.5 | 207 / 95 / 100 / 98 |
| 3D sphere, depth 9 (2.3M) | 252 / 14 / 26 / 12 | 53 / 20.7 / 27.2 / 25.0 | 629 / 614 / 561 / 604 |
| 3D river, depth 10 (0.37M) | 446 / 18 / 44 / 14 | 45 / 21.0 / 26.7 / 24.5 | 657 / 624 / 587 / 612 |

(memindex has no cursor: its dfs and box columns are plain lookups.)

- **Compact nano-vdb takes the least memory** on every non-uniform tree,
  12-14 B per node, and on sparse single levels (12-22 B per point on
  shells and rivers). Its one weak case is scattered 3D points on deep
  levels (7.4 KB a point at level 20: the 4 KiB of mask bits of a 3D upper
  node), which a balanced tree does not have.
- **Dense nano-vdb** costs 26-44 B per node (its 32^3-slot upper tables),
  frees fastest, and is hopeless for scattered 3D points (about 300 KB
  each).
- **The whole-path accessor makes nano-vdb the fastest on the 3D neighbour
  box**: 11-19% faster than memindex, where voxel's bottom-only cursor is
  2-13% faster.
- **Traversal** with cursors: the three trees are close (voxel usually
  first), all 1.7-5.5 times faster than memindex.

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

**Measured since:** F (`voxel-memindex`) against A, in isolation (see
[Measured: memindex against voxel-memindex](#measured-memindex-against-voxel-memindex)).
F wins on memory and build time on every level that is not full. A wins on
random lookups, by 1.3–3 times. This revises point 1: replace A by F, once a
run of `adapt` and the ghosts on a sparse case confirms that the lookups stay
off the hot path (Future plans, Tree structure).
