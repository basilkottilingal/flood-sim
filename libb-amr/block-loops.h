#ifndef _BAMR_BLOCK_LOOPS_H_
#define _BAMR_BLOCK_LOOPS_H_

  /*
  ..  Traversals of the cells of a block, as cwrap @loops: a reference.
  ..
  ..  A block has m interior cells per side, m a power of 2 (2^blockdepth),
  ..  and g = GHOSTS ghost layers on every side, stored with it: libb-amr
  ..  stores (m + 2 g)^dimension cells a block (more with padding), and
  ..  ghosts_update () (block-ghosts.h) fills the ghosts from the neighbours.
  ..  A stencil reads its neighbours in place. The loops here work on any
  ..  block, as a BlockShape says what it holds.
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
  ..              shape.cells Reals. The pool aligns a family's object of
  ..              blocks to 64 bytes (a page, from 16 KiB), and shape.cells
  ..              is then a multiple of V, so every block of the family
  ..              starts on a V boundary too.
  ..
  ..  Stored blocks have the shape tree_block_shape (): compact by default,
  ..  aligned to BAMR_ALIGN Reals with -DBAMR_ALIGN=V (tree.h). In 3D with
  ..  m = 16, g = 2, compact is 20^3 cells; aligned to 4 doubles, rows of 24
  ..  and 24 x 20 x 20 cells, 20% more, for aligned vector loads.
  ..
  ..  Every loop here reads the layout from a BlockShape, so code written with
  ..  them works with either, and with any other that keeps i fastest. Code
  ..  that needs the interior coordinates of a cell (0 to m - 1, ghosts
  ..  outside) takes ci - shape.ox, cj - shape.g, ck - shape.g; block_at ()
  ..  goes the other way.
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
  ..  A stencil reads the ghosts in place: cell_at (a, -1, 0) of a cell on
  ..  the edge of the interior is a ghost, valid after ghosts_update ().
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
  ..  many blocks, and each is a natural unit of cache. Once ghosts_update ()
  ..  has filled the ghosts, a stencil reads them in place, and writes
  ..  another scalar (a block written while its neighbours read it would
  ..  race):
  ..
  ..      ghosts_update (u);
  ..      BlockShape shape = tree_block_shape ();
  ..      foreach_leaf_parallel () {
  ..        const Real * a = block_data (sblocks, u);
  ..        Real * v = block_data (sblocks, t);
  ..        block_cells_simd (shape, 0)
  ..          v [c] = cell_at (a, 1, 0) + cell_at (a, -1, 0) + ... ;
  ..      }
  ..
  ..  block_cells_threads, threads inside one block, only pays when there are
  ..  fewer blocks than cores and the blocks are large (3D, m >= 32).
  ..
  ..  Cache. With m = 16 in 3D, a stored block with 2 ghost layers is 20^3
  ..  doubles, 62 KiB (half in float): a block, or two of them for a stencil
  ..  from one scalar into another, stays in L2 while the stencil runs on it.
  ..  From m = 64 (68^3 doubles, 2.4 MiB) planes of k stop fitting, and
  ..  block_cells_tiled, which walks (j, k) in tiles, keeps the k-1, k, k+1
  ..  planes of a tile in cache.
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
  ..  Stored blocks hold their ghosts, so after ghosts_update () on the host
  ..  a stencil runs on the stored blocks as they are: gpu_cells over the
  ..  table, with tree_block_shape (), halo 0.
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
  #include "libutils/bits.h"

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
    assert (m > 0 && bits_is_pow2_u32 ((uint32_t) m) && g >= 0 && v > 0);
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

  /* compact: (m + 2 g)^dimension, no padding */
  static inline BlockShape block_shape (int m, int g)
  {
    return block_shape_aligned (m, g, 1);
  }

  /* the shape of every stored block: m, GHOSTS layers, BAMR_ALIGN */
  static inline BlockShape tree_block_shape (void)
  {
    return block_shape_aligned (tree_get_block_m (), GHOSTS, BAMR_ALIGN);
  }

  /* the stride of axis a (0, 1, 2) in a block: 0 beyond the dimension */
  static inline long block_stride (BlockShape b, int a)
  {
    return a == 0 ? 1 : a == 1 ? b.sy : b.sz;
  }

  /*
  .. block_at (shape, x): the index of the cell at interior coordinates
  .. x [0..dimension) (0 to m - 1 inside, -g to -1 and m to m + g - 1 the
  .. ghosts); x [a] for a >= dimension is not read.
  */
  static inline long block_at (BlockShape b, const int x [3])
  {
    long c = x [0] + b.ox;
  #if dimension > 1
    c += (long) (x [1] + b.g) * b.sy;
  #endif
  #if dimension > 2
    c += (long) (x [2] + b.g) * b.sz;
  #endif
    return c;
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
  .. node and sblocks. 'continue' goes to the next leaf. It reads the tree's
  .. list of leaves (tree_caches (), tree.h), made again only after the tree
  .. changes: no walk. The tree must not change in it.
  */
  @loop foreach_leaf () : IJKL, node, sblocks
  {
    const TreeCaches * tc = tree_caches ();
    const TreeItem * items = tc->leaves;
    unsigned generation = tree_generation;
    for (int q = 0, n = tc->nleaves; q < n; q++)
    {
      if (tree_generation != generation)
        tree_changed_under ("foreach_leaf");
      INDEX_UNPACK (items [q].p);
      BAMRNode * node = items [q].node;
      void *** sblocks = node->s;
      (void) node, (void) sblocks;
      {...}
    }
  }

  /*
  .. foreach_leaf_parallel () { ... }: the leaf blocks shared out among
  .. threads, dynamically, as blocks differ in cost, in Morton order, so
  .. that neighbouring blocks tend to go to the same thread. Each block must
  .. be thread safe: no break, and nothing shared written. Without OpenMP,
  .. foreach_leaf () in Morton order.
  ..
  .. The list is the tree's, kept between changes (tree_caches ()): no walk,
  .. no gathering of its own.
  */
  @loop foreach_leaf_parallel () : IJKL, node, sblocks
  {
    const TreeCaches * tc = tree_caches ();
    const TreeItem * items = tc->leaves;
    int n = tc->nleaves;
    BLOCK_PRAGMA (omp parallel for schedule (dynamic))
    for (int b = 0; b < n; b++)
    {
      INDEX_UNPACK (items [b].p);
      BAMRNode * node = items [b].node;
      void *** sblocks = node->s;
      (void) node, (void) sblocks;
      {...}
    }
  }

  /*
  .. The blocks of one level l of the tree, for multigrid, from the tree's
  .. lists (tree_caches ()), with IJKL, node and sblocks:
  ..
  ..   foreach_level_block (l)          every node of level l, leaves and
  ..                                    parents
  ..   foreach_level_or_leaf_block (l)  the grid at depth l: the leaves
  ..                                    coarser than l, coarsest first, then
  ..                                    every node of level l
  ..
  .. Below tree_get_depth_min () the two are the same, every level being
  .. full. They give ghosts_by_caller: in the stencil loops over them
  .. (foreach_level (l), foreach_level_or_leaf (l), below) cwrap writes no
  .. ghost update or mark; the caller updates the ghosts of the level
  .. (ghosts_update_level (), block-ghosts.h) and marks what it wrote. The
  .. tree must not change in them.
  */
  static inline int block_level_clamp (int l)
  {
    return l < 0 ? -1 : l > BAMR_MAX_DEPTH ? BAMR_MAX_DEPTH : l;
  }

  @loop foreach_level_block (int l) : IJKL, node, sblocks, ghosts_by_caller
  {
    const TreeCaches * tc = tree_caches ();
    unsigned generation = tree_generation;
    int ghosts_by_caller = 1, L = block_level_clamp (l);
    (void) ghosts_by_caller;
    int q0 = L < 0 ? 0 : tc->level_first [L], q1 = L < 0 ? 0 : tc->level_first [L + 1];
    for (int q = q0; q < q1; q++)
    {
      if (tree_generation != generation)
        tree_changed_under ("foreach_level");
      INDEX_UNPACK (tc->by_level [q].p);
      BAMRNode * node = tc->by_level [q].node;
      void *** sblocks = node->s;
      (void) node, (void) sblocks;
      {...}
    }
  }

  @loop foreach_level_or_leaf_block (int l) : IJKL, node, sblocks, ghosts_by_caller
  {
    const TreeCaches * tc = tree_caches ();
    unsigned generation = tree_generation;
    int ghosts_by_caller = 1, L = block_level_clamp (l);
    (void) ghosts_by_caller;
    /* the leaves of levels 0 .. L - 1, then the nodes of level L */
    int nleaf = L < 0 ? 0 : tc->leaf_first [L] - tc->leaf_first [0];
    int nlevel = L < 0 ? 0 : tc->level_first [L + 1] - tc->level_first [L];
    for (int q = 0; q < nleaf + nlevel; q++)
    {
      if (tree_generation != generation)
        tree_changed_under ("foreach_level_or_leaf");
      const TreeItem * it = q < nleaf ? &tc->leaves_by_level [q]
                                      : &tc->by_level [tc->level_first [L] + q - nleaf];
      INDEX_UNPACK (it->p);
      BAMRNode * node = it->node;
      void *** sblocks = node->s;
      (void) node, (void) sblocks;
      {...}
    }
  }

  /*
  .. Stencil loops, for cwrap's s[] (see libc-pp/README.md, "Stencils"):
  ..
  ..   foreach () { ... }    every interior cell of every leaf, the ghosts
  ..                         kept up to date as cells are written
  ..
  .. A stencil loop does not nest, in itself or another, nor through a call.
  .. Loops that fill ghosts depend on what is beyond the block (a neighbour,
  .. a coarser leaf, the domain's boundary, another rank): block-ghosts.h
  .. does that, not a stencil loop.
  */
  @foreach foreach () = foreach_leaf () foreach_cell (tree_block_shape (), 0);

  /*
  .. On one level, for multigrid (see foreach_level_block above): every
  .. interior cell of every node of level l, or of the grid at depth l. The
  .. rules of foreach () hold, offsets and their reach checked, but the
  .. ghosts are the caller's: ghosts_update_level (s, l) before reading a
  .. neighbour, and ghosts_mark_dirty (s) after writing s, if a foreach ()
  .. reads its neighbours later.
  */
  /*
  .. The colours of a checkerboard, for red-black Gauss-Seidel: a cell is
  .. red if the sum of its global indices (at its level) is even, black if
  .. odd; BAMR_BOTH is every cell, in a fixed order.
  ..
  .. foreach_cell_colour (shape, halo, colour, parity): the cells of one
  .. colour, every other cell of a row, or every cell, i fastest. parity is
  .. that of the block's first interior cell, m (i + j + k) mod 2
  .. (block_parity ()): the colour is the global one whatever m. It gives
  .. checkerboard: in the stencil loops over it, cwrap lets a scalar the loop
  .. writes be read at an odd sum of offsets only, the other colour, which
  .. the loop does not change (a static assertion, as the reach). With
  .. BAMR_BOTH the same reads are in place, Gauss-Seidel in the order of the
  .. cells. No vector pragma: the order within a block is kept, so the
  .. result is the same for any threads.
  */
  enum { BAMR_RED = 0, BAMR_BLACK = 1, BAMR_BOTH = 2 };

  #if dimension == 1
  static inline int block_parity (int m, int i) { return (m * i) & 1; }
  #elif dimension == 2
  static inline int block_parity (int m, int i, int j) { return (m * (i + j)) & 1; }
  #else
  static inline int block_parity (int m, int i, int j, int k)
  {
    return (m * (i + j + k)) & 1;
  }
  #endif

  @loop foreach_cell_colour (BlockShape shape, int halo, int colour, int parity)
    : CIJK, c, shape, checkerboard
  {
    int checkerboard = 1;
    (void) checkerboard;
    assert (colour >= BAMR_RED && colour <= BAMR_BOTH);
    int h = block_halo (shape, halo), step = colour == BAMR_BOTH ? 1 : 2;
  #if dimension > 2
    for (int ck = shape.g - h; ck < shape.g + shape.m + h; ck++)
  #endif
  #if dimension > 1
    for (int cj = shape.g - h; cj < shape.g + shape.m + h; cj++)
  #endif
    {
      int row = parity, base = 0;
  #if dimension > 1
      row += cj - shape.g, base += cj * shape.sy;
  #endif
  #if dimension > 2
      row += ck - shape.g, base += ck * shape.sz;
  #endif
      /* the first x of the colour: x + row has its parity */
      int first = step == 1 ? -h : -h + (((-h + row - colour) % 2) + 2) % 2;
      for (int x = first; x < shape.m + h; x += step)
      {
        int ci = shape.ox + x;
        int c = base + ci;
        (void) c;
        {...}
      }
    }
  }

  /*
  .. Faces (libc-pp's face_vector: F.x at a cell is its lower face along x).
  .. A block's faces along axis a are at a = 0 .. m, face m in the ghost
  .. layer. Each face is made once, by its owner: a block makes its faces 0
  .. to m - 1, and face m too where no node of its level is above it along a
  .. (a wall, or a coarser leaf). So at a level jump the fine side makes the
  .. faces; the coarse side's face there is made too, and replaced by their
  .. mean when the ghosts of the face vector are made (block_ghosts ()):
  .. the fluxes agree, and what leaves one cell enters the next.
  ..
  ..   foreach_leaf_faces (axis)   the leaves, with face_last: 1 if the
  ..                               block makes its face m along axis
  ..   foreach_face_cells (shape, axis, last)  the faces of a block along
  ..                               axis, ci, cj, ck and c those of the cell
  ..                               whose lower face it is; gives face_axis
  ..
  .. foreach_face_x () (_y, _z), and foreach_face (): every face of the
  .. leaves, by axis, the body written for x and rotated (cwrap's
  .. stencil_by_axis ()); x, y, z are the face's centre.
  */
  static inline int block_neighbour_at_level (Index p, int axis, int dir)
  {
    int n = 1 << p.level, in = 1;
    if (axis == 0)
      p.i += dir, in = p.i >= 0 && p.i < n;
  #if dimension > 1
    else if (axis == 1)
      p.j += dir, in = p.j >= 0 && p.j < n;
  #endif
  #if dimension > 2
    else
      p.k += dir, in = p.k >= 0 && p.k < n;
  #endif
    return in && index_is_allocated (p);
  }

  @loop foreach_leaf_faces (int axis) : IJKL, node, sblocks, face_last
  {
    const TreeCaches * tc = tree_caches ();
    const TreeItem * items = tc->leaves;
    unsigned generation = tree_generation;
    for (int q = 0, n = tc->nleaves; q < n; q++)
    {
      if (tree_generation != generation)
        tree_changed_under ("foreach_face");
      INDEX_UNPACK (items [q].p);
      BAMRNode * node = items [q].node;
      void *** sblocks = node->s;
      int face_last = !block_neighbour_at_level (items [q].p, axis, 1);
      (void) node, (void) sblocks, (void) face_last;
      {...}
    }
  }

  @loop foreach_face_cells (BlockShape shape, int axis, int last)
    : CIJK, c, shape, face_axis
  {
    int face_axis = axis;
    (void) face_axis;
    int hi [3] = { 1, 1, 1 };
    for (int a = 0; a < dimension; a++)
      hi [a] = shape.m + (a == axis ? last : 0);
    int x [3] = { 0, 0, 0 };
  #if dimension > 2
    for (x [2] = 0; x [2] < hi [2]; x [2]++)
  #endif
  #if dimension > 1
    for (x [1] = 0; x [1] < hi [1]; x [1]++)
  #endif
    for (x [0] = 0; x [0] < hi [0]; x [0]++)
    {
      int ci = shape.ox + x [0];
  #if dimension > 1
      int cj = shape.g + x [1];
      (void) cj;
  #endif
  #if dimension > 2
      int ck = shape.g + x [2];
      (void) ck;
  #endif
      int c = (int) block_at (shape, x);
      (void) ci, (void) c;
      {...}
    }
  }

  @foreach foreach_face_x ()
    = foreach_leaf_faces (0) foreach_face_cells (tree_block_shape (), 0, face_last);
  #if dimension > 1
  @foreach foreach_face_y ()
    = foreach_leaf_faces (1) foreach_face_cells (tree_block_shape (), 1, face_last);
  #endif
  #if dimension > 2
  @foreach foreach_face_z ()
    = foreach_leaf_faces (2) foreach_face_cells (tree_block_shape (), 2, face_last);
  #endif

  /*
  .. Vertices (libc-pp's vertex_scalar: v at a cell is its lower corner).
  .. A block's vertices are at 0 .. m along each axis, m in the ghost layer.
  .. Each has one owner: the block where its indices are all below m; and at
  .. the domain's upper walls, the block there for its vertices at m along
  .. them. foreach_vertex () visits those. The ghosts of a vertex scalar
  .. then give every other copy the owner's value. Within a level that is a
  .. copy. Where a coarse leaf is above fine blocks, their vertices on it are
  .. the coarse ones: copied where they coincide, interpolated between
  .. (prolong_vertex). Where fine blocks are above, they are the parent's,
  .. the fine ones (restriction_vertex). So a vertex field is continuous,
  .. level jumps included. x, y, z are the vertex's position.
  */
  static inline int block_upper_walls (Index p)
  {
    int n = 1 << p.level, w = (p.i == n - 1);
  #if dimension > 1
    w |= (p.j == n - 1) << 1;
  #endif
  #if dimension > 2
    w |= (p.k == n - 1) << 2;
  #endif
    return w;
  }

  @loop foreach_leaf_vertices () : IJKL, node, sblocks, vertex_walls
  {
    const TreeCaches * tc = tree_caches ();
    const TreeItem * items = tc->leaves;
    unsigned generation = tree_generation;
    for (int q = 0, n = tc->nleaves; q < n; q++)
    {
      if (tree_generation != generation)
        tree_changed_under ("foreach_vertex");
      INDEX_UNPACK (items [q].p);
      BAMRNode * node = items [q].node;
      void *** sblocks = node->s;
      int vertex_walls = block_upper_walls (items [q].p);
      (void) node, (void) sblocks, (void) vertex_walls;
      {...}
    }
  }

  @loop foreach_vertex_cells (BlockShape shape, int walls)
    : CIJK, c, shape, at_vertex
  {
    int at_vertex = 1;
    (void) at_vertex;
    int hi [3] = { 1, 1, 1 };
    for (int a = 0; a < dimension; a++)
      hi [a] = shape.m + ((walls >> a) & 1);
    int x [3] = { 0, 0, 0 };
  #if dimension > 2
    for (x [2] = 0; x [2] < hi [2]; x [2]++)
  #endif
  #if dimension > 1
    for (x [1] = 0; x [1] < hi [1]; x [1]++)
  #endif
    for (x [0] = 0; x [0] < hi [0]; x [0]++)
    {
      int ci = shape.ox + x [0];
  #if dimension > 1
      int cj = shape.g + x [1];
      (void) cj;
  #endif
  #if dimension > 2
      int ck = shape.g + x [2];
      (void) ck;
  #endif
      int c = (int) block_at (shape, x);
      (void) ci, (void) c;
      {...}
    }
  }

  @foreach foreach_vertex ()
    = foreach_leaf_vertices () foreach_vertex_cells (tree_block_shape (), vertex_walls);

  @foreach foreach_level (int l)
    = foreach_level_block (l) foreach_cell (tree_block_shape (), 0);
  @foreach foreach_level_or_leaf (int l)
    = foreach_level_or_leaf_block (l) foreach_cell (tree_block_shape (), 0);

  /*
  .. One colour (BAMR_RED, BAMR_BLACK) of the grid at depth l, or every cell
  .. of it in order (BAMR_BOTH): see foreach_cell_colour above. The leaves
  .. are the grid at depth tree_get_depth (). Ghosts are the caller's, as in
  .. foreach_level (l): ghosts_update_level (s, l) between the colours.
  */
  @foreach foreach_level_or_leaf_colour (int l, int colour)
    = foreach_level_or_leaf_block (l)
      foreach_cell_colour (tree_block_shape (), 0, colour,
                           block_parity (tree_block_shape ().m, IJK));

  /*
  .. The ghost layers, as a constant that outlives the preprocessor: cwrap
  .. checks against it, at compile time, that every offset of s[...] and the
  .. radius of foreach_neighbor () is an integer constant that stays within
  .. them, so that no stencil access leaves its block.
  */
  enum { BAMR_GHOSTS = GHOSTS };

  /*
  .. block_table (s, table, where): the blocks of scalar s in the leaves, in
  .. Morton order, into table [0..], and their indices into where [0..] if
  .. not NULL. Returns how many: tree_get_nleaves (). For gpu_cells (): with
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
