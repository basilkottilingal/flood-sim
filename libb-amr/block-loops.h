#ifndef _BAMR_BLOCK_LOOPS_H_
#define _BAMR_BLOCK_LOOPS_H_

  /*
  ..  Traversals of the cells of a block, as cwrap @loops: a reference.
  ..
  ..  A block has m cells per side, m a power of 2 (2^blockdepth). libb-amr
  ..  stores m^dimension of them, without ghosts. A stencil runs on a scratch
  ..  block of (m + 2 g)^dimension cells, g ghost layers (GHOSTS) on every
  ..  side, which block_gather () fills: see block-ghosts.h. The loops here
  ..  work on either, as a BlockShape says which.
  ..
  ..
  ..  1. Layout
  ..  ---------
  ..
  ..  A block is a plain array of Reals (real.h: double, or float with
  ..  -DBAMR_SINGLE), i fastest:
  ..
  ..      cell (ci, cj, ck)  =  block [ci + cj * sy + ck * sz]
  ..
  ..  so that a row along i is contiguous: what SIMD units, hardware
  ..  prefetchers and GPU memory coalescing all want. ci, cj, ck are storage
  ..  indices, ghosts included. Along j and k the interior is [g, g + m). Along
  ..  i it is [ox, ox + m), where ox depends on the layout:
  ..
  ..    compact   ox = g, rows of m + 2 g. Rows start anywhere, so vector
  ..              loads are unaligned.
  ..
  ..    aligned   ox and the row length are rounded up to the vector width V,
  ..              in Reals: 4 doubles or 8 floats fill a 32-byte AVX vector.
  ..              Every interior row starts on a V boundary, and as m is a
  ..              power of 2, m / V vectors cover it with no remainder. It
  ..              costs some padding, (ox - g) cells before each row and the
  ..              rest after, and the block must then be allocated with
  ..              shape.cells Reals. The pool aligns blocks of 64 bytes or more
  ..              to 64 bytes (a page, from 16 KiB), and block_scratch_new ()
  ..              aligns scratch blocks the same, so vectors up to 64 bytes
  ..              stay aligned.
  ..
  ..  A stored block has g = 0, so ox = 0 and rows of m: as m is a power of 2,
  ..  a multiple of V, both layouts are the same there, aligned, with no
  ..  padding. The choice only matters for scratch blocks, whose ghosts move
  ..  the interior off the vector boundary in the compact layout.
  ..
  ..  libb-amr stores blocks without ghosts: block_shape (m, 0), m^dimension
  ..  cells. A stencil that reads neighbours runs on a scratch block with
  ..  ghosts, block_shape (m, g), that block_gather () (block-ghosts.h) fills
  ..  from the stored block and its neighbours; block_interior () then gives
  ..  each interior cell its index in both.
  ..
  ..  Every loop here reads the layout from a BlockShape, so code written with
  ..  them works with either, and with any other that keeps i fastest.
  ..
  ..
  ..  2. The loops
  ..  ------------
  ..
  ..  Each gives its block the storage indices CIJK (ci, cj, ck as the
  ..  dimension has them), the linear index c of the cell, and the shape, so
  ..  that cell_at (a, di, dj, dk) reads a neighbour. halo is how many ghost
  ..  layers to visit too: 0 for the interior, up to g for the whole block.
  ..
  ..    loop                   unit of work   break      the block must ...
  ..    ---------------------  -------------  ---------  ---------------------
  ..    block_cells            a cell         leaves     -
  ..    block_cells_simd       a cell         no         have no dependence
  ..                                                     along i
  ..    block_rows             a row          leaves     loop over the row
  ..    block_cells_tiled      a cell (3D     no         -
  ..                           tiles)
  ..    block_cells_threads    a cell         no         be thread safe, write
  ..                                                     only cell c
  ..    gpu_cells              a cell of      no         run on the device: no
  ..                           every block               host pointers, no I/O
  ..    foreach_cell           whichever BLOCK_BACKEND chooses
  ..
  ..  'continue' goes to the next cell (or row) in all of them.
  ..
  ..  Over blocks: foreach_leaf () and foreach_leaf_parallel () visit every
  ..  leaf block; block_table () gathers a scalar's blocks for gpu_cells ().
  ..
  ..
  ..  3. Which one
  ..  ------------
  ..
  ..  CPU, one core. block_cells_simd. A cell update that reads neighbours and
  ..  writes cell c of another array vectorises along i. A Gauss-Seidel sweep,
  ..  which reads cells it has just written along i, does not: use
  ..  block_cells, or red-black ordering with block_cells_simd.
  ..
  ..  CPU, explicit SIMD. block_rows gives a row at a time, to be handled with
  ..  GCC vector extensions or intrinsics:
  ..
  ..      typedef Real v4 __attribute__ ((vector_size (32)));  // 4 doubles
  ..                                                           // or 8 floats
  ..      block_rows (shape, 0)
  ..        for (int x = 0; x < n; x += 32 / sizeof (Real)) {  // n = m: no
  ..                                                           // remainder
  ..          v4 a; memcpy (&a, &u [row + x], sizeof a);
  ..          ...
  ..        }
  ..
  ..  With the aligned layout and halo 0, every row starts aligned: the
  ..  memcpy becomes an aligned load.
  ..
  ..  CPU, many cores. Parallelise over blocks, not inside them: there are
  ..  many blocks, and each is a natural unit of cache. A stencil reads
  ..  neighbours, which a stored block does not hold, so each thread gathers
  ..  a scratch block of its own: foreach_leaf_gathered_parallel
  ..  (block-ghosts.h) runs blocks on all threads, with w the scratch block
  ..  and shape its shape, and block_interior inside each:
  ..
  ..      foreach_leaf_gathered_parallel (s, 1) {
  ..        Real * v = block_data (sblocks, t);
  ..        block_interior (shape)
  ..          v [o] = cell_at (w, 1, 0) + cell_at (w, -1, 0) + ... ;
  ..      }
  ..
  ..  A pointwise update, which reads no neighbour, needs no gather:
  ..  foreach_leaf_parallel with block_cells_simd (tree_block_shape (), 0) on
  ..  the stored block.
  ..
  ..  block_cells_threads, threads inside one block, only pays when there are
  ..  fewer blocks than cores and the blocks are large (3D, m >= 32).
  ..
  ..  Cache. With m = 16 in 3D, a stored block is 16^3 doubles, 32 KiB, and a
  ..  scratch block with 2 ghost layers 20^3, 62 KiB (half each in float):
  ..  the scratch block, one per thread, stays in L2 while the stencil runs on
  ..  it. From m = 64 (68^3 doubles with ghosts, 2.4 MiB) planes of k stop
  ..  fitting, and block_cells_tiled, which walks (j, k) in tiles, keeps the
  ..  k-1, k, k+1 planes of a tile in cache.
  ..
  ..  GPU. Launch once for all blocks, not once per block: a launch costs
  ..  microseconds, a block of 16^3 cells is a few microseconds of work.
  ..  gpu_cells offloads blocks x cells as one OpenMP target region, one
  ..  thread per cell; with i fastest, neighbouring threads read neighbouring
  ..  values, and loads coalesce. The block's data must be on the device:
  ..
  ..    - gather the blocks of each scalar into a table, block_table (), and
  ..      index them as u [blk][c]; with unified shared memory
  ..      (#pragma omp requires unified_shared_memory, on recent NVIDIA and
  ..      AMD GPUs) host pool memory is usable as it is;
  ..    - otherwise copy the blocks of a scalar into one device array of
  ..      nblocks * shape.cells Reals, index it u [blk * shape.cells + c],
  ..      and map it with #pragma omp target data map (...) around the loops.
  ..
  ..  Stored blocks have no ghosts, so a stencil on the GPU needs them
  ..  gathered first: block_gather () each block into its slot of a device
  ..  array of nblocks scratch blocks (on the host, then copied; or on the
  ..  device, with the neighbours of each block found beforehand), then
  ..  gpu_cells over that array with the scratch shape, halo 0. A pointwise
  ..  update runs on the stored blocks directly.
  ..
  ..  CUDA or HIP, which cwrap cannot compile (its output goes to gcc), map
  ..  the same way: one CUDA block per AMR block, a thread per cell of a
  ..  (m, m/T) slab, the slab and its halo staged in shared memory and the
  ..  block stepping along k. The aligned layout makes each row a whole
  ..  number of 128-byte transactions when m >= 16.
  ..
  ..  Without OpenMP every pragma here expands to nothing (or GCC's ivdep),
  ..  and every loop runs serially on the host, with the same results: code
  ..  written for gpu_cells or block_cells_threads still compiles and runs.
  ..
  ..
  ..  4. Building
  ..  -----------
  ..
  ..    -O2 or -O3                 so that the simd loops vectorise
  ..    -fopenmp                   threads, simd and offload pragmas
  ..    -fopenmp-simd and          simd pragmas only, no threads
  ..      -DBLOCK_OPENMP_SIMD
  ..    -foffload=nvptx-none       a GPU with gcc (or clang -fopenmp-targets=)
  ..    -march=native              the widest vectors of this machine
  ..    -DBLOCK_BACKEND=...        BLOCK_BACKEND_SERIAL, _SIMD (default) or
  ..                               _THREADS, for foreach_cell
  */

  #include <assert.h>
  #include "tree.h"
  #include "grid-scalar.h"
  #include "morton-curve.h"

  /* the shape of a block: see 1. Layout */
  typedef struct
  {
    int m;        /* interior cells per side, a power of 2             */
    int g;        /* ghost layers per side                             */
    int ox;       /* storage ci of the first interior cell: g, or more */
    int nx;       /* storage cells along i, padding included           */
    int sy, sz;   /* strides of cj and ck: 0 where the dimension has none */
    long cells;   /* Reals in a block, to allocate                     */
  } BlockShape;

  static inline int block_round_up (int n, int a)
  {
    return (n + a - 1) / a * a;
  }

  /* vectors of v Reals: 1 gives the compact layout */
  static inline BlockShape block_shape_aligned (int m, int g, int v)
  {
    assert (m > 0 && (m & (m - 1)) == 0 && g >= 0 && v > 0);
    BlockShape b;
    b.m = m;
    b.g = g;
    b.ox = block_round_up (g, v);
    b.nx = block_round_up (b.ox + m + g, v);
    int ny = m + 2 * g;
    b.sy = dimension > 1 ? b.nx : 0;
    b.sz = dimension > 2 ? b.nx * ny : 0;
    b.cells = (long) b.nx * (dimension > 1 ? ny : 1) * (dimension > 2 ? ny : 1);
    return b;
  }

  /* compact: (m + 2 g)^dimension, no padding. With g = 0, the layout of a
  .. stored block: tree_block_shape () in block-ghosts.h */
  static inline BlockShape block_shape (int m, int g)
  {
    return block_shape_aligned (m, g, 1);
  }

  /* a scalar's block in a node, as BAMR_NODE_ENV () gives sblocks */
  static inline Real * block_data (void *** sblocks, int s)
  {
    return (Real *) sblocks [s / SCALAR_ARRAY_WIDTH][s % SCALAR_ARRAY_WIDTH];
  }

  /* the storage indices a cell loop gives its block, and a row loop */
  #if dimension == 1
    #define CIJK               ci
    #define ROW_NAMES          row, n
    #define cell_at(a, di)     ((a) [c + (di)])
  #elif dimension == 2
    #define CIJK               ci, cj
    #define ROW_NAMES          cj, row, n
    #define cell_at(a, di, dj) ((a) [c + (di) + (dj) * shape.sy])
  #else
    #define CIJK               ci, cj, ck
    #define ROW_NAMES          cj, ck, row, n
    #define cell_at(a, di, dj, dk)                                      \
      ((a) [c + (di) + (dj) * shape.sy + (dk) * shape.sz])
  #endif

  /*
  .. Pragmas: a pragma when the compiler takes it, nothing otherwise, so
  .. that every loop also builds without OpenMP.
  */
  #ifdef _OPENMP
    #define BLOCK_PRAGMA(...)  _Pragma (#__VA_ARGS__)
  #else
    #define BLOCK_PRAGMA(...)
  #endif

  #if defined (_OPENMP) || defined (BLOCK_OPENMP_SIMD)
    #define BLOCK_SIMD         _Pragma ("omp simd")
  #elif defined (__clang__)
    #define BLOCK_SIMD         _Pragma ("clang loop vectorize (enable)")
  #elif defined (__GNUC__)
    #define BLOCK_SIMD         _Pragma ("GCC ivdep")
  #else
    #define BLOCK_SIMD
  #endif

  /* target offload from OpenMP 4.5 */
  #if defined (_OPENMP) && _OPENMP >= 201511 && !defined (BLOCK_NO_OFFLOAD)
    #define BLOCK_OFFLOAD 1
  #else
    #define BLOCK_OFFLOAD 0
  #endif

  /* how many ghost layers to visit, checked */
  static inline int block_halo (BlockShape shape, int halo)
  {
    assert (halo >= 0 && halo <= shape.g);
    return halo;
  }

  /*
  .. block_cells (shape, halo): one cell at a time, in storage order.
  .. break leaves the block: each row sets 'stop' while a cell runs and
  .. clears it in its step, which break skips.
  */
  @loop block_cells (BlockShape shape, int halo) : CIJK, c, shape
  {
    int h = block_halo (shape, halo);
    int ilo = shape.ox - h, ihi = shape.ox + shape.m + h;
    int jlo = shape.g - h, jhi = shape.g + shape.m + h;
    int stop = 0;
    (void) jlo, (void) jhi, (void) stop;
  #if dimension > 2
    for (int ck = jlo; ck < jhi && !stop; ck++)
  #endif
  #if dimension > 1
    for (int cj = jlo; cj < jhi && !stop; cj++)
  #endif
    {
      int base = 0
  #if dimension > 1
        + cj * shape.sy
  #endif
  #if dimension > 2
        + ck * shape.sz
  #endif
        ;
      for (int ci = ilo; ci < ihi; stop = 0, ci++)
      {
        int c = base + ci;
        stop = 1;
        (void) c;
        {...}
      }
    }
  }

  /*
  .. block_cells_simd (shape, halo): rows vectorised along i. The cells of a
  .. row may run in any order, or at once: no break, and no reading along i
  .. what the row writes.
  */
  @loop block_cells_simd (BlockShape shape, int halo) : CIJK, c, shape
  {
    int h = block_halo (shape, halo);
    int ilo = shape.ox - h, ihi = shape.ox + shape.m + h;
    int jlo = shape.g - h, jhi = shape.g + shape.m + h;
    (void) jlo, (void) jhi;
  #if dimension > 2
    for (int ck = jlo; ck < jhi; ck++)
  #endif
  #if dimension > 1
    for (int cj = jlo; cj < jhi; cj++)
  #endif
    {
      int base = 0
  #if dimension > 1
        + cj * shape.sy
  #endif
  #if dimension > 2
        + ck * shape.sz
  #endif
        ;
      BLOCK_SIMD
      for (int ci = ilo; ci < ihi; ci++)
      {
        int c = base + ci;
        (void) c;
        {...}
      }
    }
  }

  /*
  .. block_rows (shape, halo): one row along i at a time, for explicit
  .. vector code: the block gets the linear index 'row' of the first cell
  .. and the length 'n', m + 2 halo, and cj, ck of the row. With halo 0 n is
  .. m, a power of 2. break leaves the block.
  */
  @loop block_rows (BlockShape shape, int halo) : ROW_NAMES, shape
  {
    int h = block_halo (shape, halo);
    int n = shape.m + 2 * h;
    int jlo = shape.g - h, jhi = shape.g + shape.m + h;
    int stop = 0;
    (void) jlo, (void) jhi, (void) stop;
  #if dimension > 2
    for (int ck = jlo; ck < jhi && !stop; ck++)
  #endif
  #if dimension > 1
    for (int cj = jlo; cj < jhi && !stop; stop = 0, cj++)
  #else
    for (int once = 1; once; stop = 0, once = 0)
  #endif
    {
      int row = shape.ox - h
  #if dimension > 1
        + cj * shape.sy
  #endif
  #if dimension > 2
        + ck * shape.sz
  #endif
        ;
      stop = 1;
      (void) row, (void) n;
      {...}
    }
  }

  /*
  .. block_interior (shape): the interior cells of a block of that shape,
  .. rows vectorised, with c their index in it and o their index in a stored
  .. block of m^dimension cells, without ghosts: to read neighbours in a
  .. scratch block and write the result to the stored one,
  ..
  ..   block_interior (shape)
  ..     out [o] = cell_at (w, 1, 0) + ... ;
  ..
  .. Same rules as block_cells_simd.
  */
  @loop block_interior (BlockShape shape) : CIJK, c, o, shape
  {
    int m = shape.m;
    (void) m;
  #if dimension > 2
    for (int ck = shape.g; ck < shape.g + m; ck++)
  #endif
  #if dimension > 1
    for (int cj = shape.g; cj < shape.g + m; cj++)
  #endif
    {
      int base = 0, obase = - shape.ox;   /* c and o of the row at ci = 0 */
  #if dimension > 1
      base += cj * shape.sy;
      obase += (cj - shape.g) * m;
  #endif
  #if dimension > 2
      base += ck * shape.sz;
      obase += (ck - shape.g) * m * m;
  #endif
      BLOCK_SIMD
      for (int ci = shape.ox; ci < shape.ox + m; ci++)
      {
        int c = base + ci, o = obase + ci;
        (void) c, (void) o;
        {...}
      }
    }
  }

  /*
  .. block_cells_tiled (shape, halo, tile): in 3D, (j, k) in tiles of
  .. tile x tile rows, rows vectorised; in 1D and 2D a block is small enough
  .. that it is block_cells_simd. Same rules as block_cells_simd.
  */
  @loop block_cells_tiled (BlockShape shape, int halo, int tile) : CIJK, c, shape
  {
    int h = block_halo (shape, halo);
    int ilo = shape.ox - h, ihi = shape.ox + shape.m + h;
    int jlo = shape.g - h, jhi = shape.g + shape.m + h;
    int t = tile > 0 ? tile : 1;
    (void) jlo, (void) jhi, (void) t;
  #if dimension > 2
    for (int tk = jlo; tk < jhi; tk += t)
    for (int tj = jlo; tj < jhi; tj += t)
    for (int ck = tk; ck < tk + t && ck < jhi; ck++)
    for (int cj = tj; cj < tj + t && cj < jhi; cj++)
  #elif dimension > 1
    for (int cj = jlo; cj < jhi; cj++)
  #endif
    {
      int base = 0
  #if dimension > 1
        + cj * shape.sy
  #endif
  #if dimension > 2
        + ck * shape.sz
  #endif
        ;
      BLOCK_SIMD
      for (int ci = ilo; ci < ihi; ci++)
      {
        int c = base + ci;
        (void) c;
        {...}
      }
    }
  }

  /*
  .. block_cells_threads (shape, halo): the rows of one block shared out
  .. among threads, each row vectorised. The cells run at the same time: no
  .. break, no shared variable written, and each cell writes only cell c
  .. (of any array). A sum needs a reduction of its own.
  */
  @loop block_cells_threads (BlockShape shape, int halo) : CIJK, c, shape
  {
    int h = block_halo (shape, halo);
    int ilo = shape.ox - h, ihi = shape.ox + shape.m + h;
    int jlo = shape.g - h, jhi = shape.g + shape.m + h;
    (void) jlo, (void) jhi;
  #if dimension == 1
    BLOCK_PRAGMA (omp parallel for simd)
    for (int ci = ilo; ci < ihi; ci++)
    {
      int c = ci;
      (void) c;
      {...}
    }
  #else
    #if dimension == 2
    BLOCK_PRAGMA (omp parallel for)
    #else
    BLOCK_PRAGMA (omp parallel for collapse (2))
    for (int ck = jlo; ck < jhi; ck++)
    #endif
    for (int cj = jlo; cj < jhi; cj++)
    {
      int base = cj * shape.sy
    #if dimension > 2
        + ck * shape.sz
    #endif
        ;
      BLOCK_SIMD
      for (int ci = ilo; ci < ihi; ci++)
      {
        int c = base + ci;
        (void) c;
        {...}
      }
    }
  #endif
  }

  /*
  .. gpu_cells (shape, nblocks, halo): every cell of blocks 0 to nblocks - 1
  .. as one offloaded region, one thread per cell; the block gets the block
  .. number blk, to index its data: u [blk][c] from block_table (), or
  .. u [blk * shape.cells + c] from one device array. See 3. Which one.
  .. Without offload, it runs on the host, in the same order as block_cells.
  */
  @loop gpu_cells (BlockShape shape, int nblocks, int halo) : blk, CIJK, c, shape
  {
    int h = block_halo (shape, halo);
    int ilo = shape.ox - h, ihi = shape.ox + shape.m + h;
    int jlo = shape.g - h, jhi = shape.g + shape.m + h;
    (void) jlo, (void) jhi;
  #if BLOCK_OFFLOAD
    #if dimension == 1
    BLOCK_PRAGMA (omp target teams distribute parallel for collapse (2))
    #elif dimension == 2
    BLOCK_PRAGMA (omp target teams distribute parallel for collapse (3))
    #else
    BLOCK_PRAGMA (omp target teams distribute parallel for collapse (4))
    #endif
  #endif
    for (int blk = 0; blk < nblocks; blk++)
  #if dimension > 2
    for (int ck = jlo; ck < jhi; ck++)
  #endif
  #if dimension > 1
    for (int cj = jlo; cj < jhi; cj++)
  #endif
    for (int ci = ilo; ci < ihi; ci++)
    {
      int c = ci
  #if dimension > 1
        + cj * shape.sy
  #endif
  #if dimension > 2
        + ck * shape.sz
  #endif
        ;
      (void) c;
      {...}
    }
  }

  /* foreach_cell (shape, halo): the back end -DBLOCK_BACKEND chooses */
  #define BLOCK_BACKEND_SERIAL   0
  #define BLOCK_BACKEND_SIMD     1
  #define BLOCK_BACKEND_THREADS  2
  #ifndef BLOCK_BACKEND
  #define BLOCK_BACKEND  BLOCK_BACKEND_SIMD
  #endif

  @loop foreach_cell (BlockShape shape, int halo) : CIJK, c, shape
  {
  #if BLOCK_BACKEND == BLOCK_BACKEND_SERIAL
    block_cells (shape, halo)
  #elif BLOCK_BACKEND == BLOCK_BACKEND_THREADS
    block_cells_threads (shape, halo)
  #else
    block_cells_simd (shape, halo)
  #endif
      {...}
  }

  /*
  .. foreach_leaf () { ... }: every leaf block, in Morton order, with IJKL,
  .. node and sblocks. 'continue' goes to the next leaf.
  */
  @loop foreach_leaf () : IJKL, node, sblocks
  {
    morton (TREE_ROOT_INDEX)
      if (node->flags & BAMR_NODE_LEAF)
        {...}
  }

  /*
  .. foreach_leaf_parallel () { ... }: the leaf blocks shared out among
  .. threads, dynamically, as blocks differ in cost. The leaves are gathered
  .. first, in Morton order, so that neighbouring blocks tend to go to the
  .. same thread. Each block must be thread safe: no break, and nothing
  .. shared written. Without OpenMP, foreach_leaf () in Morton order.
  ..
  .. The gathering is a traversal of the tree for each use. Code that runs
  .. many of these between two refinements can keep the list in a cache
  .. whose update gathers it and mark it dirty when the tree changes.
  */
  @loop foreach_leaf_parallel () : IJKL, node, sblocks
  {
    cache (Index) leaves;
    cache_reserve (leaves, tree_nleaves ());
    morton (TREE_ROOT_INDEX)
      if (node->flags & BAMR_NODE_LEAF)
        cache_push (leaves, (Index) {IJKL});
    Index * list = leaves.a;
    int n = leaves.n;
    BLOCK_PRAGMA (omp parallel for schedule (dynamic))
    for (int b = 0; b < n; b++)
    {
      INDEX_UNPACK (list [b]);
      BAMR_NODE_ENV ();
      {...}
    }
  }

  /*
  .. block_table (s, table, where): the blocks of scalar s in the leaves, in
  .. Morton order, into table [0..], and their indices into where [0..] if
  .. not NULL. Returns how many: tree_nleaves (). For gpu_cells (): with
  .. table on the device, u [blk][c] is cell c of leaf blk.
  */
  static inline int block_table (int s, Real ** table, Index * where)
  {
    int n = 0;
    foreach_leaf () {
      if (table)
        table [n] = block_data (sblocks, s);
      if (where)
        where [n] = (Index) {IJKL};
      n++;
    }
    return n;
  }

#endif
