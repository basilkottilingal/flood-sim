/* The cost of a ghost update, piece by piece, on the CPU and on the GPU,
 * each as fast as it goes: tests/run2.sh builds and runs it.
 *
 *   bench_gpu cpu|dev [base depth blockdepth]
 *
 * The tree: uniform at level base, refined to depth on a shell (a circle
 * or sphere of radius 0.3 about the centre), so it has level jumps and
 * walls. Every leaf holds u, smooth.
 *
 * cpu: the library's own functions, on every thread (block_restrict (),
 * block_ghosts (), block_prolong (), level by level, as ghosts_update ()
 * runs them). The ghost fill is timed three ways: with the walls off (a
 * boundary condition that does nothing) and the level jumps off
 * (prolong_none), the copies alone; with the jumps on; with everything:
 * the jumps and the walls are the differences.
 *
 * dev: every node's block in one device array, U [node * cells + c], put
 * there once; kernels of this file over it, one launch per level and kind:
 *   - restriction: a parent cell the mean of its children's;
 *   - the ghosts: precomputed tables, one entry per ghost cell of every
 *     leaf (built on the host from the tree, as block_ghosts () would fill
 *     it): a copy from a neighbour of the same level, a linear
 *     prolongation from the coarser leaf at a level jump, or a wall
 *     (zero gradient), followed back to the copy or prolongation it ends
 *     in, so that the walls need no pass of their own;
 *   - prolongation: every parent's children's interiors from it, as a
 *     refinement fills them;
 *   - the stencil: v = u + k lap (u) on every leaf.
 * The upload and download of every block are timed too: what a ghost
 * update on the host would cost the GPU each step. Before any timing, the
 * results are checked against the library's on the host, to a tolerance
 * (not bit for bit: the device accumulates in Real, the library in double,
 * and fma may fuse).
 *
 * Output: lines "R dim prec mode phase ms cells", ms the median per call,
 * cells the cells written per call; and "T ..." lines about the tree. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include <time.h>
#include "block-ghosts.h"
#if dimension < 2
# error "bench_gpu: 2 or 3 dimensions"
#endif
/* an operator of its own: this file dispatches, in place of block-ghosts.c */
#include "block-ops.h"

/* a wall condition that does nothing: the ghost fill without its walls */
@boundary bc_skip (const Boundary * bc)
{
  (void) bc;
}

#include "block-ghosts-impl.h"

scalar u, v, w;

static const char * prec = sizeof (Real) == 8 ? "double" : "float";
static const char * mode;
static BlockShape store;          /* tree_block_shape (): each function
                                       takes it as a local, shape, as the
                                       device reads no file-scope variable */
static long cells;                  /* Reals a stored block */
static int depth, dmin;

#if dimension == 2
# define LAP(a) (cell_at (a, 1, 0) + cell_at (a, -1, 0) + cell_at (a, 0, 1) \
                 + cell_at (a, 0, -1) - 4 * cell_at (a, 0, 0))
#else
# define LAP(a) (cell_at (a, 1, 0, 0) + cell_at (a, -1, 0, 0) + cell_at (a, 0, 1, 0) \
                 + cell_at (a, 0, -1, 0) + cell_at (a, 0, 0, 1) + cell_at (a, 0, 0, -1) \
                 - 6 * cell_at (a, 0, 0, 0))
#endif
static const Real kdiff = (Real) 0.1;

static double now (void)
{
  struct timespec t;
  clock_gettime (CLOCK_MONOTONIC, &t);
  return t.tv_sec + 1e-9 * t.tv_nsec;
}

static int cmp_double (const void * a, const void * b)
{
  double x = *(const double *) a, y = *(const double *) b;
  return x < y ? -1 : x > y;
}

/* f once to warm up (the GPU compiles its kernels on the first launch),
.. then at least 7 times and 0.5 s, at most 200 times: the median, in ms */
static double timeit (void (* f) (void))
{
  enum { NMAX = 200 };
  double t [NMAX];
  int n = 0;
  f ();
  double t0 = now ();
  while (n < NMAX && (n < 7 || now () - t0 < 0.5)) {
    double a = now ();
    f ();
    t [n++] = now () - a;
  }
  qsort (t, n, sizeof (double), cmp_double);
  return 1e3 * t [n / 2];
}

static void report (const char * phase, double ms, long ncells)
{
  printf ("R %d %s %s %s %.4f %ld\n", dimension, prec, mode, phase, ms, ncells);
  fflush (stdout);
}

static void coords (Index p, int P [3])
{
  P [0] = p.i, P [1] = p.j, P [2] = 0;
#if dimension > 2
  P [2] = p.k;
#endif
}

/* ---------------------------------------------------------------------- */
/* the tree                                                               */
/* ---------------------------------------------------------------------- */

static const double radius = 0.3;

static int touches_shell (Index p)
{
  int P [3];
  coords (p, P);
  double h = 1. / (1 << p.level), near = 0, far = 0;
  for (int a = 0; a < dimension; a++) {
    double lo = P [a] * h - 0.5, hi = lo + h;
    double n = lo > 0 ? lo : hi < 0 ? -hi : 0, f = fmax (fabs (lo), fabs (hi));
    near += n * n, far += f * f;
  }
  return near <= radius * radius && far >= radius * radius;
}

static void build (int base, int dmax, int blockdepth)
{
  tree (dmax, blockdepth);
  for (int l = 0; l < dmax; l++) {
    const TreeCaches * tc = tree_caches ();
    int n = 0;
    Index * list = malloc ((size_t) tc->nleaves * sizeof (Index));
    for (int q = 0; q < tc->nleaves; q++) {
      Index p = tc->leaves [q].p;
      if (p.level == l && (l < base || touches_shell (p)))
        list [n++] = p;
    }
    for (int q = 0; q < n; q++)
      tree_refine (list [q]);     /* 0 if balance refined it already */
    free (list);
  }
}

/* u smooth on every leaf; every node and ghost of u, v, w defined */
static void fill (void)
{
  BlockShape shape = store;
  int m = shape.m;
  foreach_leaf () {
    Real * a = block_data (sblocks, u);
    int P [3];
    coords ((Index) {IJKL}, P);
    double n = (double) (m << level);
    block_cells (shape, 0) {
      double x = (P [0] * m + ci - shape.ox + 0.5) / n;
      double y = (P [1] * m + cj - shape.g + 0.5) / n;
      double r = sin (6.3 * x) * cos (4.1 * y);
#if dimension > 2
      double z = (P [2] * m + ck - shape.g + 0.5) / n;
      r *= cos (2.7 * z);
#endif
      a [c] = (Real) r;
    }
  }
  tree_restrict (u);
  ghosts_update_levels (u);
  const TreeCaches * tc = tree_caches ();
  for (int q = 0; q < tc->nall; q++) {
    void *** s = tc->all [q].node->s;
    memcpy (block_data (s, v), block_data (s, u), cells * sizeof (Real));
    memcpy (block_data (s, w), block_data (s, u), cells * sizeof (Real));
  }
}

/* ---------------------------------------------------------------------- */
/* cpu: the library's functions, on every thread                          */
/* ---------------------------------------------------------------------- */

static void cpu_restrict (void)
{
  const TreeCaches * tc = tree_caches ();
  for (int l = depth - 1; l >= dmin; l--) {
    BLOCK_PRAGMA (omp parallel for schedule (dynamic))
    for (int q = tc->parent_first [l]; q < tc->parent_first [l + 1]; q++)
      block_restrict (tc->parents_by_level [q].p, u);
  }
}

static void cpu_fill (void)
{
  const TreeCaches * tc = tree_caches ();
  for (int l = 0; l <= depth; l++) {
    BLOCK_PRAGMA (omp parallel for schedule (dynamic))
    for (int q = tc->leaf_first [l]; q < tc->leaf_first [l + 1]; q++)
      block_ghosts (tc->leaves_by_level [q].p, u);
  }
}

static void cpu_update (void)
{
  ghosts_update (u);
}

static void cpu_prolong (void)
{
  const TreeCaches * tc = tree_caches ();
  for (int l = dmin; l < depth; l++) {
    BLOCK_PRAGMA (omp parallel for schedule (dynamic))
    for (int q = tc->parent_first [l]; q < tc->parent_first [l + 1]; q++)
      block_prolong (tc->parents_by_level [q].p, w);
  }
}

static void cpu_stencil (void)
{
  BlockShape shape = store;
  foreach_leaf_parallel () {
    const Real * a = block_data (sblocks, u);
    Real * b = block_data (sblocks, v);
    block_cells_simd (shape, 0)
      b [c] = a [c] + kdiff * LAP (a);
  }
}

static void cpu_step (void)
{
  ghosts_update (u);
  cpu_stencil ();
}

/* ---------------------------------------------------------------------- */
/* dev: the blocks on the device                                          */
/* ---------------------------------------------------------------------- */

enum { NC = TREE_NCHILDREN };
enum { K_COPY, K_WALLC, K_JUMP, K_WALLP, NK };   /* copy kinds, then
                                                    prolongation kinds */

static int nall, nleaves;
static const TreeItem * nodes;      /* tc->by_level: node n is nodes [n] */
static Real * U, * V, * W;          /* nall * cells each                 */
static long NU;
static int * LF;                    /* the leaves' node numbers          */
static int * PL, * CH;              /* the parents', by level; their
                                       children's, NC each               */
static int * ED, * ES, * EG;        /* ghost entries: destination, source
                                       (global cell indices), signs      */
static int eoff [BAMR_MAX_DEPTH + 2][NK + 1];   /* entries of level l and
                                                   kind k: eoff [l][k] ..
                                                   eoff [l][k + 1]       */
static long nentries [NK];

/* node numbers from pointers: a sorted table */
typedef struct { const BAMRNode * node; int n; } NodeId;
static NodeId * ids;

static int cmp_id (const void * a, const void * b)
{
  const BAMRNode * x = ((const NodeId *) a)->node, * y = ((const NodeId *) b)->node;
  return x < y ? -1 : x > y;
}

static int node_number (Index p)
{
  NodeId key = { index_node (p), 0 };
  NodeId * f = bsearch (&key, ids, nall, sizeof (NodeId), cmp_id);
  assert (f);
  return f->n;
}

/* where a ghost cell's value comes from: a cell (global index), and the
.. signs of a linear prolongation from it (0: a copy) */
typedef struct { int src, sg, kind, wall; } Rec;

static int pack_signs (const int sgn [3])
{
  int k = 0;
  for (int a = 0; a < dimension; a++)
    k |= (sgn [a] + 1) << (2 * a);
  return k;
}

static int is_interior (const int x [3])
{
  BlockShape shape = store;
  for (int a = 0; a < dimension; a++)
    if (x [a] < 0 || x [a] >= shape.m)
      return 0;
  return 1;
}

/* every ghost cell of leaf p, as block_ghosts_cells () and block_walls ()
.. would fill it with the default operators */
static void leaf_recipes (Index p, int n, Rec * rec)
{
  BlockShape shape = store;
  int m = shape.m, g = shape.g, l = p.level, nb = 1 << l, P [3];
  coords (p, P);
  for (long c = 0; c < cells; c++)
    rec [c].kind = -1, rec [c].wall = 0;
  enum { ndir = dimension == 2 ? 9 : 27 };
  for (int dir = 0; dir < ndir; dir++) {
    int d [3] = { 0, 0, 0 }, q [3] = { 0, 0, 0 };
    int lo [3] = { 0, 0, 0 }, hi [3] = { 1, 1, 1 }, t = dir, centre = 1, outside = 0;
    for (int a = 0; a < dimension; a++) {
      d [a] = t % 3 - 1, t /= 3;
      centre &= d [a] == 0;
      q [a] = P [a] + d [a];
      outside |= q [a] < 0 || q [a] >= nb;
      lo [a] = d [a] < 0 ? -g : d [a] == 0 ? 0 : m;
      hi [a] = d [a] < 0 ? 0 : d [a] == 0 ? m : m + g;
    }
    if (centre || outside)
      continue;
    Index b = { q [0], q [1],
#if dimension > 2
                q [2],
#endif
                l };
    int x [3] = { 0, 0, 0 };
    if (index_is_allocated (b)) {
      long base = (long) node_number (b) * cells;
      for (x [2] = lo [2]; x [2] < hi [2]; x [2]++)
      for (x [1] = lo [1]; x [1] < hi [1]; x [1]++)
      for (x [0] = lo [0]; x [0] < hi [0]; x [0]++) {
        int y [3];
        for (int a = 0; a < 3; a++)
          y [a] = a < dimension ? x [a] - d [a] * m : 0;
        Rec * r = &rec [block_at (shape, x)];
        r->src = (int) (base + block_at (shape, y)), r->sg = 0, r->kind = K_COPY;
      }
      continue;
    }
    Index C = b;
    C.level = l - 1, C.i >>= 1, C.j >>= 1;
#if dimension > 2
    C.k >>= 1;
#endif
    assert (l > 0 && index_is_allocated (C));
    int Cq [3];
    coords (C, Cq);
    long base = (long) node_number (C) * cells;
    for (x [2] = lo [2]; x [2] < hi [2]; x [2]++)
    for (x [1] = lo [1]; x [1] < hi [1]; x [1]++)
    for (x [0] = lo [0]; x [0] < hi [0]; x [0]++) {
      int X [3] = { 0, 0, 0 }, sgn [3] = { 0, 0, 0 };
      for (int a = 0; a < dimension; a++) {
        int F = (P [a] - 2 * Cq [a]) * m + x [a];
        assert (F >= 0 && F < 2 * m);
        X [a] = F >> 1, sgn [a] = F & 1 ? 1 : -1;
      }
      Rec * r = &rec [block_at (shape, x)];
      r->src = (int) (base + block_at (shape, X)), r->sg = pack_signs (sgn);
      r->kind = K_JUMP;
    }
  }
  /* the walls, one axis after the other: a ghost takes what its mirror
   * inside has; one whose mirror is not made yet gets it from a later axis */
  long own = (long) n * cells;
  for (int a = 0; a < dimension; a++)
    for (int side = 0; side < 2; side++) {
      if (side ? P [a] != nb - 1 : P [a] != 0)
        continue;
      int lo [3] = { 0, 0, 0 }, hi [3] = { 1, 1, 1 };
      for (int e = 0; e < dimension; e++)
        lo [e] = -g, hi [e] = m + g;
      for (int layer = 0; layer < g; layer++) {
        lo [a] = side ? m + layer : -1 - layer, hi [a] = lo [a] + 1;
        int x [3] = { 0, 0, 0 };
        for (x [2] = lo [2]; x [2] < hi [2]; x [2]++)
        for (x [1] = lo [1]; x [1] < hi [1]; x [1]++)
        for (x [0] = lo [0]; x [0] < hi [0]; x [0]++) {
          int y [3] = { x [0], x [1], x [2] };
          y [a] = side ? 2 * m - 1 - x [a] : -1 - x [a];
          Rec * r = &rec [block_at (shape, x)];
          if (is_interior (y))
            r->src = (int) (own + block_at (shape, y)), r->sg = 0, r->kind = K_COPY;
          else
            *r = rec [block_at (shape, y)];
          r->wall = 1;
        }
      }
    }
}

/* growable lists of entries, one per level and kind */
typedef struct { int * d, * s, * g; long n, cap; } List;

static void list_add (List * L, int d, int s, int g)
{
  if (L->n == L->cap) {
    L->cap = L->cap ? 2 * L->cap : 1024;
    L->d = realloc (L->d, L->cap * sizeof (int));
    L->s = realloc (L->s, L->cap * sizeof (int));
    L->g = realloc (L->g, L->cap * sizeof (int));
  }
  L->d [L->n] = d, L->s [L->n] = s, L->g [L->n] = g, L->n++;
}

static void dev_setup (void)
{
  BlockShape shape = store;
  const TreeCaches * tc = tree_caches ();
  nall = tc->nall, nleaves = tc->nleaves, nodes = tc->by_level;
  NU = (long) nall * cells;
  assert (NU < INT_MAX);
  ids = malloc (nall * sizeof (NodeId));
  for (int n = 0; n < nall; n++)
    ids [n] = (NodeId) { nodes [n].node, n };
  qsort (ids, nall, sizeof (NodeId), cmp_id);

  U = malloc (NU * sizeof (Real)), V = malloc (NU * sizeof (Real));
  W = malloc (NU * sizeof (Real));
  for (int n = 0; n < nall; n++) {
    void *** s = nodes [n].node->s;
    memcpy (U + n * cells, block_data (s, u), cells * sizeof (Real));
    memcpy (V + n * cells, block_data (s, v), cells * sizeof (Real));
    memcpy (W + n * cells, block_data (s, w), cells * sizeof (Real));
  }
  LF = malloc (nleaves * sizeof (int));
  int nl = 0, np = tc->parent_first [depth + 1] - tc->parent_first [0];
  for (int n = 0; n < nall; n++)
    if (nodes [n].node->flags & BAMR_NODE_LEAF)
      LF [nl++] = n;
  assert (nl == nleaves);
  PL = malloc ((np + 1) * sizeof (int)), CH = malloc ((np + 1) * NC * sizeof (int));
  for (int q = 0; q < np; q++) {
    Index p = tc->parents_by_level [q].p;
    PL [q] = node_number (p);
    for (int k = 0; k < NC; k++)
      CH [q * NC + k] = node_number (index_child (p, k));
  }

  /* the ghost entries, level by level */
  List lists [BAMR_MAX_DEPTH + 1][NK];
  memset (lists, 0, sizeof lists);
  Rec * rec = malloc (cells * sizeof (Rec));
  for (int q = 0; q < nleaves; q++) {
    int n = LF [q], l = nodes [n].p.level;
    leaf_recipes (nodes [n].p, n, rec);
    int x [3] = { 0, 0, 0 }, G = shape.g, m = shape.m;
    for (x [2] = dimension > 2 ? -G : 0; x [2] < (dimension > 2 ? m + G : 1); x [2]++)
    for (x [1] = -G; x [1] < m + G; x [1]++)
    for (x [0] = -G; x [0] < m + G; x [0]++) {
      if (is_interior (x))
        continue;
      long c = block_at (shape, x);
      Rec r = rec [c];
      assert (r.kind >= 0);
      int k = r.wall ? (r.kind == K_JUMP ? K_WALLP : K_WALLC) : r.kind;
      list_add (&lists [l][k], (int) (n * cells + c), r.src, r.sg);
    }
  }
  free (rec);
  long total = 0;
  for (int l = 0; l <= depth; l++)
    for (int k = 0; k < NK; k++)
      total += lists [l][k].n, nentries [k] += lists [l][k].n;
  ED = malloc ((total + 1) * sizeof (int)), ES = malloc ((total + 1) * sizeof (int));
  EG = malloc ((total + 1) * sizeof (int));
  long e = 0;
  for (int l = 0; l <= depth + 1; l++)
    for (int k = 0; k < NK; k++) {
      eoff [l][k] = (int) e;
      if (l > depth)
        continue;
      List * L = &lists [l][k];
      memcpy (ED + e, L->d, L->n * sizeof (int));
      memcpy (ES + e, L->s, L->n * sizeof (int));
      memcpy (EG + e, L->g, L->n * sizeof (int));
      e += L->n;
      free (L->d), free (L->s), free (L->g);
    }
  for (int l = 0; l <= depth; l++)
    eoff [l][NK] = eoff [l + 1][0];
  long ne = total, npc = (long) np * NC;
  (void) ne, (void) npc;
#pragma omp target enter data map (to: U [0:NU], V [0:NU], W [0:NU], LF [0:nleaves], \
                                   PL [0:np + 1], CH [0:npc + NC], ED [0:ne + 1], \
                                   ES [0:ne + 1], EG [0:ne + 1])
}

static void dev_upload (void)
{
#pragma omp target update to (U [0:NU])
}

static void dev_download (void)
{
#pragma omp target update from (U [0:NU])
}

static void dev_restrict (void)
{
  BlockShape shape = store;
  const TreeCaches * tc = tree_caches ();
  Real * u_ = U;
  const int * pl = PL, * ch = CH;
  long cl = cells;
  int h = shape.m / 2, m = shape.m, sy = shape.sy, ox = shape.ox, g = shape.g;
  (void) g;
#if dimension > 2
  int sz = shape.sz;
#endif
  for (int l = depth - 1; l >= dmin; l--) {
    int q0 = tc->parent_first [l], np = tc->parent_first [l + 1] - q0;
    if (!np)
      continue;
    gpu_cells (shape, np, 0) {
      int x0 = ci - ox, x1 = cj - g, k = (x0 >= h) | (x1 >= h) << 1;
      int f0 = 2 * x0 - (x0 >= h) * m, f1 = 2 * x1 - (x1 >= h) * m;
#if dimension > 2
      int x2 = ck - g;
      k |= (x2 >= h) << 2;
      int f2 = 2 * x2 - (x2 >= h) * m;
#endif
      const Real * f = u_ + ch [(q0 + blk) * NC + k] * cl + f0 + ox + (f1 + g) * sy;
#if dimension == 2
      Real s = (f [0] + f [1]) + (f [sy] + f [sy + 1]);
      u_ [pl [q0 + blk] * cl + c] = s * (Real) 0.25;
#else
      f += (f2 + g) * sz;
      Real s = ((f [0] + f [1]) + (f [sy] + f [sy + 1]))
               + ((f [sz] + f [sz + 1]) + (f [sz + sy] + f [sz + sy + 1]));
      u_ [pl [q0 + blk] * cl + c] = s * (Real) 0.125;
#endif
    }
  }
}

/* entries e0 .. e1 - 1: copies */
static void dev_copies (int e0, int e1)
{
  Real * u_ = U;
  const int * ed = ED, * es = ES;
  if (e1 <= e0)
    return;
#pragma omp target teams distribute parallel for simd
  for (int e = e0; e < e1; e++)
    u_ [ed [e]] = u_ [es [e]];
}

/* entries e0 .. e1 - 1: linear prolongations */
static void dev_jumps (int e0, int e1)
{
  BlockShape shape = store;
  Real * u_ = U;
  const int * ed = ED, * es = ES, * eg = EG;
  int sy = shape.sy;
#if dimension > 2
  int sz = shape.sz;
#endif
  if (e1 <= e0)
    return;
#pragma omp target teams distribute parallel for simd
  for (int e = e0; e < e1; e++) {
    int s = es [e], k = eg [e];
    const Real * a = u_ + s;
    Real r = (Real) ((k & 3) - 1) * (a [1] - a [-1])
             + (Real) (((k >> 2) & 3) - 1) * (a [sy] - a [-sy]);
#if dimension > 2
    r += (Real) (((k >> 4) & 3) - 1) * (a [sz] - a [-sz]);
#endif
    u_ [ed [e]] = a [0] + r * (Real) 0.125;
  }
}

/* the whole fill: two launches a level, the copies (walls ending in a
.. copy too), then the prolongations (walls ending in one too) */
static void dev_fill (void)
{
  for (int l = 0; l <= depth; l++) {
    dev_copies (eoff [l][K_COPY], eoff [l][K_JUMP]);
    dev_jumps (eoff [l][K_JUMP], eoff [l][NK]);
  }
}

/* each kind alone, for the split */
static void dev_fill_copy (void)
{
  for (int l = 0; l <= depth; l++)
    dev_copies (eoff [l][K_COPY], eoff [l][K_WALLC]);
}

static void dev_fill_jump (void)
{
  for (int l = 0; l <= depth; l++)
    dev_jumps (eoff [l][K_JUMP], eoff [l][K_WALLP]);
}

static void dev_fill_wall (void)
{
  for (int l = 0; l <= depth; l++) {
    dev_copies (eoff [l][K_WALLC], eoff [l][K_JUMP]);
    dev_jumps (eoff [l][K_WALLP], eoff [l][NK]);
  }
}

static void dev_update (void)
{
  dev_restrict ();
  dev_fill ();
}

static void dev_prolong (void)
{
  BlockShape shape = store;
  const TreeCaches * tc = tree_caches ();
  Real * w_ = W;
  const int * pl = PL, * ch = CH;
  long cl = cells;
  int m = shape.m, sy = shape.sy, ox = shape.ox, g = shape.g;
#if dimension > 2
  int sz = shape.sz;
#endif
  for (int l = dmin; l < depth; l++) {
    int q0 = tc->parent_first [l], np = tc->parent_first [l + 1] - q0;
    if (!np)
      continue;
    /* block blk: child blk % NC of parent blk / NC */
    gpu_cells (shape, np * NC, 0) {
      int q = q0 + blk / NC, k = blk % NC;
      int X0 = (k & 1) * m + ci - ox, X1 = ((k >> 1) & 1) * m + cj - g;
      int C = (X0 >> 1) + ox + ((X1 >> 1) + g) * sy;
#if dimension > 2
      int X2 = ((k >> 2) & 1) * m + ck - g;
      C += ((X2 >> 1) + g) * sz;
#endif
      const Real * a = w_ + pl [q] * cl + C;
      Real r = (Real) (X0 & 1 ? 1 : -1) * (a [1] - a [-1])
               + (Real) (X1 & 1 ? 1 : -1) * (a [sy] - a [-sy]);
#if dimension > 2
      r += (Real) (X2 & 1 ? 1 : -1) * (a [sz] - a [-sz]);
#endif
      w_ [ch [q * NC + k] * cl + c] = a [0] + r * (Real) 0.125;
    }
  }
}

static void dev_stencil (void)
{
  BlockShape shape = store;
  const Real * u_ = U;
  Real * v_ = V;
  const int * lf = LF;
  long cl = cells;
  Real kd = kdiff;
  gpu_cells (shape, nleaves, 0) {
    const Real * a = u_ + lf [blk] * cl;
    v_ [lf [blk] * cl + c] = a [c] + kd * LAP (a);
  }
}

static void dev_step (void)
{
  dev_update ();
  dev_stencil ();
}

/* ---------------------------------------------------------------------- */
/* the device against the library, before timing                          */
/* ---------------------------------------------------------------------- */

static int fails;

/* the largest |device - host| of scalar s's blocks in D: over the leaves
.. (leaves 1), the parents (0), or every node finer than dmin (-1); the
.. whole block or the interior */
static double diff (const Real * D, scalar s, int leaves, int interior)
{
  BlockShape shape = store;
  double e = 0;
  for (int n = 0; n < nall; n++) {
    int leaf = nodes [n].node->flags & BAMR_NODE_LEAF;
    if (leaves >= 0 && leaf != (leaves ? BAMR_NODE_LEAF : 0))
      continue;
    if (leaves < 0 && nodes [n].p.level <= dmin)
      continue;
    const Real * h = block_data (nodes [n].node->s, s), * d = D + n * cells;
    block_cells (shape, interior ? 0 : shape.g)
      e = fmax (e, fabs ((double) d [c] - h [c]));
  }
  return e;
}

static void check (const char * what, double e)
{
  double tol = 64 * REAL_EPSILON;   /* |u| <= 1 */
  printf ("C %d %s %s: max |device - host| %.3g%s\n", dimension, prec, what, e,
          e <= tol ? "" : " FAIL");
  fails += e > tol;
}

static void dev_check (void)
{
  BlockShape shape = store;
  /* the parents and the leaves' ghosts: u spoilt first on the device,
   * every ghost and parent, so that nothing passes by being left as is */
  for (int n = 0; n < nall; n++) {
    int leaf = nodes [n].node->flags & BAMR_NODE_LEAF;
    Real * d = U + n * cells;
    block_cells (shape, shape.g) {
      int x [3] = { ci - shape.ox, cj - shape.g, 0 };
#if dimension > 2
      x [2] = ck - shape.g;
#endif
      if (!leaf || !is_interior (x))
        d [c] = (Real) 7;
    }
  }
  dev_upload ();
  dev_update ();
  dev_download ();
  ghosts_update (u);
  check ("ghosts of the leaves", diff (U, u, 1, 0));
  /* parents below dmin are not restricted, by either */
  double e = 0;
  for (int n = 0; n < nall; n++)
    if (!(nodes [n].node->flags & BAMR_NODE_LEAF) && nodes [n].p.level >= dmin) {
      const Real * h = block_data (nodes [n].node->s, u), * d = U + n * cells;
      block_cells (shape, 0)
        e = fmax (e, fabs ((double) d [c] - h [c]));
    }
  check ("restriction", e);
  dev_stencil ();
  cpu_stencil ();
#pragma omp target update from (V [0:NU])
  check ("stencil", diff (V, v, 1, 1));
  /* one prolongation from W = u on both */
  dev_prolong ();
  cpu_prolong ();
#pragma omp target update from (W [0:NU])
  check ("prolongation", diff (W, w, -1, 1));
}

/* ---------------------------------------------------------------------- */

int main (int argc, char ** argv)
{
  mode = argc > 1 ? argv [1] : "cpu";
  int dev = !strcmp (mode, "dev");
  int base = argc > 2 ? atoi (argv [2]) : dimension == 2 ? 4 : 2;
  int dmax = argc > 3 ? atoi (argv [3]) : dimension == 2 ? 8 : 5;
  int bd = argc > 4 ? atoi (argv [4]) : dimension == 2 ? 5 : 4;
  build (base, dmax, bd);
  store = tree_block_shape ();
  BlockShape shape = store;
  cells = shape.cells;
  depth = tree_get_depth (), dmin = tree_get_depth_min ();
  const TreeCaches * tc = tree_caches ();
  long m_d = 1;
  for (int a = 0; a < dimension; a++)
    m_d *= shape.m;
  int npar = 0;                       /* the parents restricted: >= dmin */
  for (int l = dmin; l < depth; l++)
    npar += tc->parent_first [l + 1] - tc->parent_first [l];
  printf ("T %d %s nodes %d leaves %d parents %d (from level %d) levels %d-%d m %d g %d"
          " cells/block %ld MB/scalar %.0f\n",
          dimension, prec, tc->nall, tc->nleaves, npar, dmin, dmin, depth, shape.m,
          shape.g, cells, tc->nall * cells * sizeof (Real) / 1e6);
  fill ();
  dev_setup ();
  printf ("T %d %s ghost cells: %ld copies, %ld level jumps, %ld walls (%ld of them"
          " ending in a jump)\n", dimension, prec, nentries [K_COPY], nentries [K_JUMP],
          nentries [K_WALLC] + nentries [K_WALLP], nentries [K_WALLP]);
  long ncopy = nentries [K_COPY], njump = nentries [K_JUMP];
  long nwall = nentries [K_WALLC] + nentries [K_WALLP];
  long nghost = ncopy + njump + nwall, nres = npar * m_d, ninner = tc->nleaves * m_d;

  if (dev) {
    dev_check ();
    report ("upload", timeit (dev_upload), NU);
    report ("download", timeit (dev_download), NU);
    report ("restriction", timeit (dev_restrict), nres);
    report ("ghosts:copies", timeit (dev_fill_copy), ncopy);
    report ("ghosts:jumps", timeit (dev_fill_jump), njump);
    report ("ghosts:walls", timeit (dev_fill_wall), nwall);
    report ("ghosts", timeit (dev_fill), nghost);
    report ("ghosts_update", timeit (dev_update), nres + nghost);
    report ("prolongation", timeit (dev_prolong), nres * NC);
    report ("stencil", timeit (dev_stencil), ninner);
    report ("step", timeit (dev_step), nres + nghost + ninner);
  }
  else {
    report ("restriction", timeit (cpu_restrict), nres);
    /* the copies alone, then with the jumps, then everything */
    ghost_ops (u, prolongation (prolong_none), boundary.all (bc_skip));
    double t_copy = timeit (cpu_fill);
    scalar_ops_reset (u);
    ghost_ops (u, boundary.all (bc_skip));
    double t_jump = timeit (cpu_fill);
    scalar_ops_reset (u);
    double t_all = timeit (cpu_fill);
    report ("ghosts:copies", t_copy, ncopy);
    report ("ghosts:jumps", t_jump - t_copy, njump);
    report ("ghosts:walls", t_all - t_jump, nwall);
    report ("ghosts", t_all, nghost);
    report ("ghosts_update", timeit (cpu_update), nres + nghost);
    report ("prolongation", timeit (cpu_prolong), nres * NC);
    report ("stencil", timeit (cpu_stencil), ninner);
    report ("step", timeit (cpu_step), nres + nghost + ninner);
  }
  tree_destroy ();
  return fails != 0;
}
