/* A grid for the foreach () tests: what libb-amr's <block-loops.h> gives
 * cwrap, on NLEAF fixed blocks of 4 x 4 cells with one ghost layer, stored
 * 6 x 6, i fastest. Build with -dimension=2. */
typedef double Real;
typedef struct { int m, g, ox, nx, sy, sz; long cells; } BlockShape;

#define NLEAF 3
#define NS 8
#define NX 6

static Real store [NLEAF][NS][NX * NX];
static void * slots [NLEAF][NS];
static void ** groups [NLEAF][1];

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

@loop foreach_leaf () : level, sblocks
{
  for (int level = 0; level < NLEAF; level++) {
    void *** sblocks = leaf_blocks (level);
    {...}
  }
}

@loop foreach_cell (BlockShape shape, int halo) : c, shape
{
  for (int cj = shape.g - halo; cj < shape.g + shape.m + halo; cj++)
    for (int ci = shape.ox - halo; ci < shape.ox + shape.m + halo; ci++) {
      int c = ci + cj * shape.sy;
      {...}
    }
}
