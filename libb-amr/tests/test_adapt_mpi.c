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
 * never split. The distributed run is made three times: with
 * tree_mpi_imbalance_max 1 (repartition whenever a process has more
 * leaves than a cut gives), 1.5 and 1e9 (only for a family of leaves split
 * by the merges): the same counts and cells each time; after every step
 * the partition table gives every node's pid, and no process has more
 * leaves than the threshold allows; the large thresholds repartition
 * less often. Last, a family split by a merge (split_case ()): only a
 * repartition lets it merge as on one process. Build with -D_MPI. */
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

/* the partition made since the tree changed, on every process: what
 * a repartition leaves (on one process alone it may only mean that its
 * tree did not change) */
static int repartitioned (void)
{
  int v = tree_partition_valid (), all;
  MPI_Allreduce (&v, &all, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
  return all;
}

/* the partition table gives the pid of every node held, its own leaves
 * this process's; after a repartition, its leaves and loads (unit
 * weights) those of every process, its first key that of each first leaf */
static int owners_by_table (void)
{
  const TreeCaches * tc = tree_caches ();
  int ok = 1;
  if (repartitioned ()) {
    const TreePartition * tp = tree_get_partition ();
    long n = tc->nowned, * all = malloc ((size_t) npe * sizeof (long));
    MPI_Allgather (&n, 1, MPI_LONG, all, 1, MPI_LONG, MPI_COMM_WORLD);
    ok &= tp->leaf_first [0] == 0;
    for (int r = 0; r < npe; r++)
      ok &= tp->leaf_first [r + 1] - tp->leaf_first [r] == all [r]
            && tp->load [r] == (double) all [r];
    free (all);
    if (n && rank > 0)
      ok &= tp->key_first [rank] == tree_morton_key (tc->owned [0].p);
  }
  for (int q = 0; q < tc->nall; q++)
    ok &= tree_partition_owner (tc->all [q].p) == tc->all [q].node->pid;
  for (int q = 0; q < tc->nowned; q++)
    ok &= tc->owned [q].node->pid == rank;
  return ok;
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

/*
.. A family split by a merge. Level 2 uniform, its first leaf P divided:
.. the leaves P's children, P's siblings s1 .. s_{2^D - 1}, then the other
.. families of level 2. Weights put P's children alone on rank 0: the cut
.. runs between them and s1 (P not a leaf, so no family rule). f constant:
.. the first adaptation merges P's children (and the other families of
.. leaves), making P, s1 ... a family of leaves across ranks 0 and 1; the
.. second merges it, which only its owner can, the family on one process:
.. adapt_wavelet () must repartition after the first even with
.. tree_mpi_imbalance_max 1e9. Against one process: the counts, every
.. leaf and cell.
*/
static void split_build (void)
{
  tree_replicated (MAXLEVEL, 4);
  int n = tree_caches ()->nleaves;
  Index * l1 = malloc ((size_t) n * sizeof (Index));
  for (int q = 0; q < n; q++)
    l1 [q] = tree_caches ()->leaves [q].p;
  for (int q = 0; q < n; q++)
    tree_refine (l1 [q]);
  free (l1);
  tree_refine (tree_caches ()->leaves [0].p);
}

static void split_case (void)
{
  if (npe < 2)
    return;
  AdaptStats want [2];
  Snapshot after [2];
  split_build ();
  for (int k = 0; k < 2; k++) {
    field (3);
    want [k] = adapt_wavelet (1, (int []) { f }, &max, 1, MAXLEVEL);
    after [k] = snapshot ();
  }
  tree_destroy ();

  tree_mpi_imbalance_max = 1e9;
  split_build ();
  const TreeCaches * tc = tree_caches ();
  int n = tc->nleaves, other = 0;
  double * w = malloc ((size_t) n * sizeof (double));
  for (int q = 0; q < n; q++)
    other += tc->leaves [q].p.level < 3;
  for (int q = 0; q < n; q++)
    w [q] = tc->leaves [q].p.level == 3
            ? (double) other / ((npe - 1) * TREE_NCHILDREN) : 1;
  CHECK (tree_mpi_partition (w) == 0);
  free (w);
  /* the cut where it was meant: rank 0 has P's children only */
  long mine = rank == 0 && tree_caches ()->nowned == TREE_NCHILDREN, ok;
  MPI_Allreduce (&mine, &ok, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
  CHECK (ok == 1);
  int counts = 1, parts = 0;
  for (int k = 0; k < 2; k++) {
    field (3);
    AdaptStats got = adapt_wavelet (1, (int []) { f }, &max, 1, MAXLEVEL);
    counts &= got.refined == want [k].refined
              && got.coarsened == want [k].coarsened;
    int part = repartitioned ();
    parts += part << k;
    if (rank == 0)
      printf ("  split family, step %d: %d merged (one process: %d); %s\n",
              k, got.coarsened, want [k].coarsened,
              part ? "repartitioned" : "kept");
    compare (&after [k], "  split family");
  }
  CHECK (counts && (parts & 1));
  tree_destroy ();
  snapshot_free (&after [0]), snapshot_free (&after [1]);
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

  /* the same, distributed, for three thresholds */
  const double eta [3] = { 1, 1.5, 1e9 };
  int parts [3];
  for (int e = 0; e < 3; e++) {
    tree_mpi_imbalance_max = eta [e];
    if (rank == 0)
      printf (" tree_mpi_imbalance_max %g\n", eta [e]);
    build ();
    CHECK (tree_mpi_partition (NULL) == 0);
    int counts = 1, balance = 1, owners = 1;
    parts [e] = 0;
    for (int k = 0; k < STEPS; k++) {
      field (k);
      AdaptStats got = adapt_wavelet (1, (int []) { f }, &max, 1, MAXLEVEL);
      counts &= got.refined == want [k].refined
                && got.coarsened == want [k].coarsened;
      /* repartitioned: the partition made since the tree changed */
      int part = repartitioned ();
      parts [e] += part;
      /* the leaves of each process: at most the threshold or a cut's */
      long n = tree_caches ()->nowned, most, total;
      MPI_Allreduce (&n, &most, 1, MPI_LONG, MPI_MAX, MPI_COMM_WORLD);
      MPI_Allreduce (&n, &total, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
      double mean = (double) total / npe;
      double cut = ceil (mean) + TREE_NCHILDREN - 1;
      balance &= most <= (eta [e] * mean > cut ? eta [e] * mean : cut);
      int ow = owners_by_table (), all;
      MPI_Allreduce (&ow, &all, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
      owners &= all;
      if (rank == 0)
        printf ("  step %d: %d divided, %d merged (one process: %d, %d); "
                "%s, most %ld, mean %.1f\n", k, got.refined, got.coarsened,
                want [k].refined, want [k].coarsened,
                part ? "repartitioned" : "kept", most, mean);
      if (k == 2)
        compare (&front, "the front");
    }
    CHECK (counts && balance && owners);
    compare (&last, "the constant");
    tree_destroy ();
  }
  /* fewer repartitions with a larger threshold */
  CHECK (parts [0] >= parts [1] && parts [1] >= parts [2]);
  CHECK (npe > 1 || parts [0] == 0);
  /* on 3 processes in 2D and 3D, means of 40 leaves and more: 1.5 keeps
   * cuts that 1 makes again */
  CHECK (npe != 3 || dimension == 1 || parts [1] < parts [0]);
  if (rank == 0)
    printf (" repartitions: %d, %d, %d of %d steps\n", parts [0], parts [1],
            parts [2], STEPS);
  snapshot_free (&front), snapshot_free (&last);
  split_case ();
  int fl = fails, ftot;
  MPI_Allreduce (&fl, &ftot, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  if (rank == 0)
    printf ("dimension %d, %d processes: %d checks, %d failed\n", dimension,
            npe, checks, ftot);
  return ftot != 0;
}
