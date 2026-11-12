/* A grid for the foreach () tests: what libb-amr's <block-loops.h> gives
 * cwrap, on NLEAF fixed blocks of 4 x 4 cells with two ghost layers, stored
 * 8 x 8, i fastest: the leaves i = 0, 1, 2, j = 0 of level 2 in the unit
 * square. Build with -dimension=2. */
typedef double Real;
typedef struct { int m, g, ox, nx, sy, sz; long cells; } BlockShape;

#define NLEAF 3
#define NS 8
#define NX 8
enum { BAMR_GHOSTS = 2 };

static Real store [NLEAF][NS][NX * NX];
static void * slots [NLEAF][NS];
static void ** groups [NLEAF][1];

static inline double tree_get_size (void) { return 1; }
static inline double tree_get_origin (int axis) { (void) axis; return 0; }
static inline int tree_get_block_m (void) { return 4; }

static inline BlockShape tree_block_shape (void)
{
  return (BlockShape) { 4, 2, 2, NX, NX, 0, NX * NX };
}

/* sblocks [s / NS][s % NS] is the block of scalar s, as in libb-amr */
static inline Real * block_data (void *** sblocks, int s)
{
  return (Real *) sblocks [s / NS][s % NS];
}

static inline void *** leaf_blocks (int l)
{
  for (int s = 0; s < NS; s++)
    slots [l][s] = store [l][s];
  groups [l][0] = slots [l];
  return groups [l];
}

/* Stale ghosts, as libb-amr's <block-ghosts.h> keeps them: here an update
 * only records which scalars it updated, in ghosts_log */
static unsigned long long ghosts_dirty;
static char ghosts_log [256];

static inline void ghosts_mark_dirty (int s) { ghosts_dirty |= 1ull << s; }
static inline void ghosts_mark_all_dirty (void) { ghosts_dirty = ~0ull; }

static inline void ghosts_update_dirty (int n, const int * s)
{
  for (int t = 0; t < 64; t++) {
    int listed = n < 0;
    for (int q = 0; q < n; q++)
      listed |= s [q] == t;
    if (listed && (ghosts_dirty >> t & 1)) {
      char * end = ghosts_log;
      while (*end)
        end++;
      end [0] = (char) ('0' + t), end [1] = ' ', end [2] = '\0';
      ghosts_dirty &= ~(1ull << t);
    }
  }
}

@loop foreach_leaf () : i, j, level, sblocks
{
  for (int i = 0; i < NLEAF; i++) {
    int j = 0, level = 2;
    (void) j, (void) level;
    void *** sblocks = leaf_blocks (i);
    {...}
  }
}

@loop foreach_cell (BlockShape shape, int halo) : ci, cj, c, shape
{
  for (int cj = shape.g - halo; cj < shape.g + shape.m + halo; cj++)
    for (int ci = shape.ox - halo; ci < shape.ox + shape.m + halo; ci++) {
      int c = ci + cj * shape.sy;
      {...}
    }
}

/* the stencil loop, as <block-loops.h> declares it */
@foreach foreach () = foreach_leaf () foreach_cell (tree_block_shape (), 0);

/* Two levels, as libb-amr's foreach_prolongation and foreach_restriction:
 * foreach_up () fills leaf 0's interior from leaf 1, as if it were the
 * coarser; foreach_down () fills leaf 2's from its "children", leaves 0, 1,
 * 2, 0, a quadrant each. */
enum { BAMR_COARSE_REACH = 1 };

@loop mock_up () : sblocks, coarse_sblocks
{
  void *** sblocks = leaf_blocks (0);
  void *** coarse_sblocks = leaf_blocks (1);
  {...}
}

@loop mock_up_cells (BlockShape shape) : c, shape, cc, sgn
{
  for (int cj = shape.g; cj < shape.g + shape.m; cj++)
    for (int ci = shape.ox; ci < shape.ox + shape.m; ci++) {
      int c = ci + cj * shape.sy;
      int cc = shape.ox + (ci - shape.ox) / 2 + (shape.g + (cj - shape.g) / 2) * shape.sy;
      int sgn [3] = { (ci - shape.ox) % 2 ? 1 : -1, (cj - shape.g) % 2 ? 1 : -1, 0 };
      (void) cc, (void) sgn;
      {...}
    }
}

@loop mock_down () : sblocks, child_sblocks, child
{
  void *** sblocks = leaf_blocks (2);
  for (int child = 0; child < 4; child++) {
    void *** child_sblocks = leaf_blocks (child % 3);
    {...}
  }
}

@loop mock_down_cells (BlockShape shape, int child) : c, shape, fc
{
  int h = shape.m / 2;
  for (int Y = 0; Y < h; Y++)
    for (int X = 0; X < h; X++) {
      int c = shape.ox + (child & 1) * h + X + (shape.g + (child >> 1) * h + Y) * shape.sy;
      int fc = shape.ox + 2 * X + (shape.g + 2 * Y) * shape.sy;
      (void) fc;
      {...}
    }
}

@foreach foreach_up () = mock_up () mock_up_cells (tree_block_shape ());
@foreach foreach_down () = mock_down () mock_down_cells (tree_block_shape (), child);

/* A loop on a level whose ghosts its caller keeps, as libb-amr's
 * foreach_level (l): the leaves here, at level l; its block loop gives
 * ghosts_by_caller, so no update or mark is written around it */
@loop mock_level (int l) : i, j, level, sblocks, ghosts_by_caller
{
  int ghosts_by_caller = 1;
  (void) ghosts_by_caller;
  for (int i = 0; i < NLEAF; i++) {
    int j = 0, level = l;
    (void) j, (void) level;
    void *** sblocks = leaf_blocks (i);
    {...}
  }
}

@foreach foreach_mock_level (int l) = mock_level (l) foreach_cell (tree_block_shape (), 0);

/* One colour of a checkerboard, as libb-amr's foreach_level_or_leaf_colour
 * (l, colour): the cells whose (i + j) has the parity colour, every other
 * cell of a row; the cell loop gives checkerboard */
@loop mock_colour_cells (BlockShape shape, int colour) : ci, cj, c, shape, checkerboard
{
  int checkerboard = 1;
  (void) checkerboard;
  for (int cj = shape.g; cj < shape.g + shape.m; cj++)
    for (int ci = shape.ox + ((cj - shape.g + colour) & 1); ci < shape.ox + shape.m; ci += 2) {
      int c = ci + cj * shape.sy;
      {...}
    }
}

@foreach foreach_mock_colour (int colour) = foreach_leaf () mock_colour_cells (tree_block_shape (), colour);
