# block based amr



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

**Gain: large, if the grid adapts often. Effort: small.**

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

**Gain: moderate. Effort: medium.**

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

**Gain: about 2× memory and ghost bandwidth in 3D. Effort: large.**

Storage per block, relative to the interior alone:

| | m = 16, g = 2 | m = 32, g = 2 |
|---|---|---|
| 1D | 1.25× | 1.13× |
| 2D | 1.56× | 1.27× |
| 3D | **1.95×** | 1.42× |

In 3D with m = 16, nearly half of every block is ghost cells, and every sweep
carries them through the cache. The options:
- **No ghosts in storage.** Blocks store m^d cells. Before computing on a
  block, a thread copies the block and its neighbours' faces into a scratch
  block with ghosts, which stays in cache. It computes there, then writes the
  interior back. That copy is the ghost filling that is needed anyway. This
  halves memory in 3D, and the bandwidth spent on ghosts, but touches every
  kernel.
- **Smaller ghost regions.** A 5- or 7-point stencil never reads edge or
  corner ghosts. Storing only face ghosts saves less, but changes less.
- **Larger m.** Cheaper ghosts, but coarser refinement.

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

- **Mixed precision.** Storing some fields as `float` halves their memory and
  bandwidth. It is a numerical decision, not only an allocation one.
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
| 6 | blocks without stored ghosts | large | ~2× memory and ghost bandwidth in 3D |

**Measure before changing anything:** time a stencil sweep, and a refine and
coarsen cycle. Use `grid_scalar_check ()` and `pool_used ()` to make sure no
change leaks blocks.
