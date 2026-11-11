/* A grid for the foreach () tests: what libb-amr's <block-loops.h> gives
 * cwrap, on NLEAF fixed blocks of 4 x 4 cells with one ghost layer, stored
 * 6 x 6, i fastest: the leaves i = 0, 1, 2, j = 0 of level 2 in the unit
 * square. Build with -dimension=2. */
typedef double Real;
typedef struct { int m, g, ox, nx, sy, sz; long cells; } BlockShape;

#define NLEAF 3
#define NS 8
#define NX 6

static Real store [NLEAF][NS][NX * NX];
static void * slots [NLEAF][NS];
static void ** groups [NLEAF][1];

static inline double tree_get_size (void) { return 1; }
static inline double tree_get_origin (int axis) { (void) axis; return 0; }
static inline int tree_get_block_m (void) { return 4; }

static inline BlockShape tree_block_shape (void)
{
  return (BlockShape) { 4, 1, 1, NX, NX, 0, NX * NX };
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
