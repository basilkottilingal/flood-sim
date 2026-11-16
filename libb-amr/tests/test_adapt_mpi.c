/* adapt_wavelet () on the distributed tree (tree-mpi.h, milestone 4),
 * against the same on one process. Each process first builds the whole
 * tree (tree_replicated ()), refined about a sphere, and adapts it three
 * times to a front that moves, as one process does: the reference, the
 * counts of each step and, at the end, every leaf with every cell of two
 * cell scalars after a ghost update; then twice more on a constant, which
 * merges every family down to minlevel. Then the same tree is made again,
 * distributed, and adapted the same way: the counts of each step the
 * reference's; every leaf of the whole tree owned by one process; every
 * cell of the two scalars on every leaf of a process (interior and
 * ghosts) the reference's, bit for bit. The front refines and coarsens,
 * with 2:1 balance across processes (rounds of requests) and families
 * never split. Build with -D_MPI. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include "tree-mpi.h"
#include "block-ghosts.h"
#include "adapt.h"

scalar f, g;

static int checks, fails, rank, npe;
#define CHECK(c) do { checks++; if (!(c)) { fails++; \
  printf ("rank %d: FAIL line %d: %s\n", rank, __LINE__, #c); } } while (0)

enum { STEPS = 5 };
static const double max = 0.01;
#define MAXLEVEL (dimension == 3 ? 4 : 5)

static int index_coord (Index p, int ax)
{
  return ax == 0 ? p.i
#if dimension > 1
       : ax == 1 ? p.j
#endif
#if dimension > 2
       : ax == 2 ? p.k
#endif
       : 0;
}

/* the leaves near the sphere of radius 0.3 about the box's centre, refined */
static void refine_sphere (int rounds)
{
  for (int r = 0; r < rounds; r++) {
    const TreeCaches * tc = tree_caches ();
    int n = tc->nleaves, k = 0;
    Index * near = malloc ((size_t) n * sizeof (Index));
    for (int q = 0; q < n; q++) {
      Index p = tc->leaves [q].p;
      double h = 1. / (1 << p.level), d2 = 0;
      for (int ax = 0; ax < dimension; ax++) {
        double x = (index_coord (p, ax) + 0.5) * h - 0.5;
        d2 += x * x;
      }
      if (fabs (sqrt (d2) - 0.3) < h)
        near [k++] = p;
    }
    for (int q = 0; q < k; q++)
      tree_refine (near [q]);
    free (near);
  }
}

/* a front of radius 0.25 about a centre that moves with the step, then
 * (steps 3 and 4) f constant: every family merges down to minlevel, and
 * the cut must gather the families the merges make across processes; g
 * linear, carried along by the prolongations and restrictions */
static void field (int step)
{
  double cx = 0.38 + 0.09 * step, front = step < 3;
  foreach () {
    double r = sqrt ((x - cx) * (x - cx)
                     + (dimension > 1 ? (y - 0.5) * (y - 0.5) : 0)
                     + (dimension > 2 ? (z - 0.5) * (z - 0.5) : 0));
    f[] = front ? tanh ((r - 0.25) / 0.01) : 1;
    g[] = 1 + x - 2 * y + 0.5 * z;
  }
}

static void build (void)
{
  tree_replicated (MAXLEVEL, 4);
  refine_sphere (dimension == 3 ? 1 : 2);
}

/* a leaf of the reference, its blocks of f and g */
typedef struct { Index p; Real * b; } Leaf;
typedef struct { Leaf * a; int n; } Snapshot;
static long cells;

static int leaf_cmp (const void * a, const void * b)
{
  Index p = ((const Leaf *) a)->p, q = ((const Leaf *) b)->p;
  uint64_t x = tree_morton_key (p), y = tree_morton_key (q);
  if (x != y)
    return x < y ? -1 : 1;
  return p.level - q.level;
}

/* every leaf of the whole tree, f and g after a ghost update */
static Snapshot snapshot (void)
{
  ghosts_update (f, g);
  const TreeCaches * tc = tree_caches ();
  Snapshot v = { malloc ((size_t) tc->nleaves * sizeof (Leaf)), tc->nleaves };
  for (int q = 0; q < v.n; q++) {
    v.a [q].p = tc->leaves [q].p;
    v.a [q].b = malloc (2 * (size_t) cells * sizeof (Real));
    memcpy (v.a [q].b, block_data (tc->leaves [q].node->s, f),
            (size_t) cells * sizeof (Real));
    memcpy (v.a [q].b + cells, block_data (tc->leaves [q].node->s, g),
            (size_t) cells * sizeof (Real));
  }
  qsort (v.a, (size_t) v.n, sizeof (Leaf), leaf_cmp);
  return v;
}

/* the distributed tree against a snapshot: every leaf of a process in it,
 * every cell of f and g (after a ghost update) bit for bit, every leaf of
 * the snapshot owned once */
static void compare (const Snapshot * v, const char * what)
{
  ghosts_update (f, g);
  const TreeCaches * tc = tree_caches ();
  long found = 0, wrong = 0, compared = 0;
  for (int q = 0; q < tc->nowned; q++) {
    Leaf key = { tc->owned [q].p, NULL };
    const Leaf * r = bsearch (&key, v->a, (size_t) v->n, sizeof (Leaf),
                              leaf_cmp);
    if (!r)
      continue;
    found++;
    const Real * bf = block_data (tc->owned [q].node->s, f);
    const Real * bg = block_data (tc->owned [q].node->s, g);
    for (long c = 0; c < cells; c++) {
      wrong += memcmp (&bf [c], &r->b [c], sizeof (Real)) != 0;
      wrong += memcmp (&bg [c], &r->b [cells + c], sizeof (Real)) != 0;
      compared += 2;
    }
  }
  long mine [4] = { tc->nowned, found, wrong, compared }, all [4];
  MPI_Allreduce (mine, all, 4, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
  CHECK (all [0] == v->n && all [1] == v->n && all [2] == 0 && all [3] > 0);
  if (rank == 0)
    printf ("  %s: %ld leaves of %d found, %ld of %ld cells differ\n", what,
            all [1], v->n, all [2], all [3]);
}

static void snapshot_free (Snapshot * v)
{
  for (int q = 0; q < v->n; q++)
    free (v->a [q].b);
  free (v->a);
}

int main (void)
{
  build ();
  MPI_Comm_rank (MPI_COMM_WORLD, &rank);
  MPI_Comm_size (MPI_COMM_WORLD, &npe);
  cells = tree_block_shape ().cells;

  /* the reference: one process's adaptations of the whole tree, the
   * front's tree (after step 2) and the last */
  AdaptStats want [STEPS];
  Snapshot front, last;
  for (int k = 0; k < STEPS; k++) {
    field (k);
    want [k] = adapt_wavelet (1, (int []) { f }, &max, 1, MAXLEVEL);
    if (k == 2)
      front = snapshot ();
  }
  last = snapshot ();
  tree_destroy ();

  /* the same, distributed */
  build ();
  CHECK (tree_mpi_partition (NULL) == 0);
  int counts = 1, balance = 1;
  for (int k = 0; k < STEPS; k++) {
    field (k);
    AdaptStats got = adapt_wavelet (1, (int []) { f }, &max, 1, MAXLEVEL);
    counts &= got.refined == want [k].refined
              && got.coarsened == want [k].coarsened;
    if (rank == 0)
      printf ("  step %d: %d divided, %d merged (one process: %d, %d)\n", k,
              got.refined, got.coarsened, want [k].refined,
              want [k].coarsened);
    /* balanced again: each process within one family of the mean */
    long n = tree_caches ()->nowned, most, least, total;
    MPI_Allreduce (&n, &most, 1, MPI_LONG, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce (&n, &least, 1, MPI_LONG, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce (&n, &total, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
    balance &= most <= (double) total / npe + TREE_NCHILDREN
               && least >= (double) total / npe - TREE_NCHILDREN;
    if (k == 2)
      compare (&front, "the front");
  }
  CHECK (counts && balance);
  compare (&last, "the constant");
  snapshot_free (&front), snapshot_free (&last);
  tree_destroy ();
  int fl = fails, ftot;
  MPI_Allreduce (&fl, &ftot, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  if (rank == 0)
    printf ("dimension %d, %d processes: %d checks, %d failed\n", dimension,
            npe, checks, ftot);
  return ftot != 0;
}
