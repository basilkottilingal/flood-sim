/* tree-mpi.h with MPI: the ghosts of a cell scalar, a face vector and a
 * vertex scalar, and the restriction of the parents, made on npe processes
 * by the plan's exchanges on the pruned trees, against the same made by
 * one process on the whole tree, bit for bit. Every block a process keeps
 * but does not own is poisoned (NaN) first, and its own leaves' ghosts
 * too: a block missing from the plan, or sent too early, leaves a NaN or a
 * different value; one pruned that is read, a crash (ASan). Leaves and
 * every node (parents too). The partition, made by rank 0 and broadcast,
 * against the same made by each process. Build with -D_MPI, run with
 * mpiexec -n 1, 2, 3, 4 ... */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include "tree-mpi.h"
#include "block-ghosts.h"

scalar a;
face_vector F;
vertex_scalar v;

static int checks, fails, rank, npe;
#define CHECK(c) do { checks++; if (!(c)) { fails++; \
  printf ("rank %d: FAIL line %d: %s\n", rank, __LINE__, #c); } } while (0)

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

static int list [1 + 3 + 1], nlist;
static TreeItem * nodes;         /* every node of the whole tree */
static int nnodes;
static char * alive;             /* still there: not pruned away */
static long cells;

/* every cell of every block of the scalars, set to value */
static void prefill (Real value)
{
  for (int q = 0; q < nnodes; q++)
    if (alive [q])
    for (int k = 0; k < nlist; k++) {
      Real * b = block_data (nodes [q].node->s, list [k]);
      for (long c = 0; c < cells; c++)
        b [c] = value;
    }
}

/* the fields, at their places: the cells, faces and vertices of leaves */
static void fill (void)
{
  foreach ()
    a[] = sin (3 * x + 1) * cos (2 * y - 0.5) + z * z;
  foreach_face ()
    F.x[] = cos (2 * x) + sin (3 * y + 0.7 * z);
  foreach_vertex ()
    v[] = x * x - 2 * y + sin (z + x);
}

/* every block of every node, the scalars one after another */
static Real * snapshot (void)
{
  Real * s = malloc ((size_t) nnodes * nlist * cells * sizeof (Real));
  for (int q = 0; q < nnodes; q++)
    if (alive [q])
    for (int k = 0; k < nlist; k++)
      memcpy (s + ((size_t) q * nlist + k) * cells,
              block_data (nodes [q].node->s, list [k]), cells * sizeof (Real));
  return s;
}

/* is cell x of a block, its interior coordinates, data the fill made (not
 * a ghost), for a scalar of location loc: the interior, and for a face the
 * face m along its axis, for a vertex those at m too */
static int owned (BlockShape sh, int loc, long c)
{
  /* the coordinates of storage cell c */
  int x [3] = { 0, 0, 0 };
  long r = c;
  x [0] = (int) (r % sh.nx) - sh.ox;
#if dimension > 1
  r /= sh.nx;
  x [1] = (int) (r % (sh.m + 2 * sh.g)) - sh.g;
#endif
#if dimension > 2
  r /= (sh.m + 2 * sh.g);
  x [2] = (int) r - sh.g;
#endif
  int in = 1;
  for (int ax = 0; ax < dimension; ax++) {
    int top = sh.m - 1;
    if (loc == CWRAP_VERTEX || loc == CWRAP_FACE_X + ax)
      top = sh.m;
    in &= x [ax] >= 0 && x [ax] <= top;
  }
  return in;
}

/* the references, on the whole tree, one process's work (every node is
 * this process's before the partition): which cells the update writes
 * (two prefills), and what */
static Real * r1 [2], * r2 [2];

static void reference (int parents)
{
  prefill (R (1e30)), fill (), ghosts_update_list (nlist, list, parents);
  r1 [parents] = snapshot ();
  prefill (R (-1e30)), fill (), ghosts_update_list (nlist, list, parents);
  r2 [parents] = snapshot ();
}

/* is a node of owner pid needed by this process: its own, or received
 * (by its address only: it may be freed) */
static int needed (const TreeMpiPlan * plan, BAMRNode * node, int pid)
{
  int need = pid == rank;
  for (int pa = 0; pa < 2 && !need; pa++)
    for (int kind = 0; kind < TREE_MPI_KINDS && !need; kind++)
      for (int l = 0; l < BAMR_MAX_DEPTH + 2 && !need; l++) {
        const TreeMpiExchange * ex = &plan->ex [pa][kind][l];
        for (int k = 0; k < ex->nrecv && !need; k++)
          for (int j = 0; j < ex->recv [k].n && !need; j++)
            need = ex->recv [k].item [j].node == node;
      }
  return need;
}

/* one update on npe processes against the same on one: parents 0 or 1 */
static void compare (int parents)
{
  BlockShape sh = tree_block_shape ();
  /* every block kept but not owned poisoned, and the ghosts of owned
   * leaves; owned parents wholly (their restriction) */
  prefill (R (-1e30)), fill ();
  for (int q = 0; q < nnodes; q++) {
    if (!alive [q])
      continue;
    BAMRNode * node = nodes [q].node;
    int mine = node->pid == rank, leaf = node->flags & BAMR_NODE_LEAF;
    for (int k = 0; k < nlist; k++) {
      Real * b = block_data (node->s, list [k]);
      int loc = scalar_location (list [k]);
      for (long c = 0; c < cells; c++)
        if (!mine || !leaf || !owned (sh, loc, c))
          b [c] = (Real) NAN;
    }
  }
  ghosts_update_list (nlist, list, parents);

  /* each owned node: every cell the update writes, as on one process */
  long compared = 0, wrong = 0, nmine = 0;
  for (int q = 0; q < nnodes; q++) {
    BAMRNode * node = nodes [q].node;
    if (!alive [q] || node->pid != rank)
      continue;
    nmine++;
    for (int k = 0; k < nlist; k++) {
      const Real * b = block_data (node->s, list [k]);
      size_t o = ((size_t) q * nlist + k) * cells;
      for (long c = 0; c < cells; c++) {
        if (memcmp (&r1 [parents][o + c], &r2 [parents][o + c], sizeof (Real)))
          continue;                 /* not written by the update */
        compared++;
        wrong += memcmp (&b [c], &r1 [parents][o + c], sizeof (Real)) != 0;
      }
    }
  }
  long all [3] = { compared, wrong, nmine }, sum [3];
  MPI_Allreduce (all, sum, 3, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
  CHECK (sum [1] == 0 && sum [0] > 0 && sum [2] == nnodes);

  /* what this process received, in blocks of the scalars */
  long blocks = 0, total;
  const TreeMpiPlan * plan = tree_mpi_plan ();
  for (int kind = 0; kind < TREE_MPI_KINDS; kind++)
    for (int l = 0; l < BAMR_MAX_DEPTH + 2; l++)
      for (int k = 0; k < plan->ex [parents][kind][l].nrecv; k++)
        blocks += plan->ex [parents][kind][l].recv [k].n;
  MPI_Allreduce (&blocks, &total, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
  if (rank == 0)
    printf ("  parents %d: %ld cells compared, %ld differ; %ld blocks "
            "received (of %d nodes)\n", parents, sum [0], sum [1], total,
            nnodes);
  free (r1 [parents]), free (r2 [parents]);
}

int main (int argc, char ** argv)
{
  MPI_Init (&argc, &argv);
  MPI_Comm_rank (MPI_COMM_WORLD, &rank);
  MPI_Comm_size (MPI_COMM_WORLD, &npe);
  tree (5, 4);
  refine_sphere (dimension == 3 ? 2 : 3);
  list [nlist++] = a;
  list [nlist++] = F.x;
#if dimension > 1
  list [nlist++] = F.y;
#endif
#if dimension > 2
  list [nlist++] = F.z;
#endif
  list [nlist++] = v;
  const TreeCaches * tc = tree_caches ();
  nnodes = tc->nall;
  nodes = malloc ((size_t) nnodes * sizeof (TreeItem));
  memcpy (nodes, tc->all, (size_t) nnodes * sizeof (TreeItem));
  cells = tree_block_shape ().cells;

  alive = malloc ((size_t) nnodes);
  memset (alive, 1, (size_t) nnodes);
  CHECK (tree_npe () == npe && tree_rank () == rank);
  int own = 1;
  for (int q = 0; q < nnodes; q++)
    own &= nodes [q].node->pid == rank;
  CHECK (own);
  reference (0);
  reference (1);

  /* the partition each process would make, to compare with rank 0's */
  int * pid = malloc ((size_t) nnodes * sizeof (int));
  int * first = malloc ((size_t) (npe + 1) * sizeof (int));
  CHECK (tree_partition (npe, NULL) == 0);
  for (int q = 0; q < nnodes; q++)
    pid [q] = nodes [q].node->pid;
  memcpy (first, tree_get_partition ()->leaf_first,
          (size_t) (npe + 1) * sizeof (int));
  double imbalance = tree_get_partition ()->imbalance;
  tree_partition_free ();
  for (int q = 0; q < nnodes; q++)
    nodes [q].node->pid = rank;   /* this process's again: rank 0 tells */

  CHECK (tree_mpi_partition (NULL) == 0);
  const TreeMpiPlan * plan = tree_mpi_plan ();
  CHECK (plan && plan->rank == rank && plan->npe == npe);
  CHECK (npe > 1 ? tree_pruned () > 0 : !tree_pruned ());
  {
    /* pruned: every node needed still there, the same node; the pids and
     * the partition rank 0's, as made here */
    const TreePartition * tp = tree_get_partition ();
    int there = 1, same = 1, nalive = 0;
    for (int q = 0; q < nnodes; q++) {
      alive [q] = index_is_allocated (nodes [q].p)
                  && index_node (nodes [q].p) == nodes [q].node;
      int need = needed (plan, nodes [q].node, pid [q]);
      there &= !need || alive [q];
      same &= !alive [q] || nodes [q].node->pid == pid [q];
      nalive += alive [q];
    }
    CHECK (there && same);
    CHECK (tp && tp->npe == npe && tp->imbalance == imbalance
           && !memcmp (tp->leaf_first, first, (size_t) (npe + 1)
                                                * sizeof (int)));
    CHECK (tree_caches ()->nall == nalive && tree_get_nblocks () == nalive);
    long kept = nalive, sum;
    MPI_Allreduce (&kept, &sum, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
    if (rank == 0)
      printf ("  %.2f of the %d nodes kept by a process\n",
              (double) sum / npe / nnodes, nnodes);
  }
  free (pid), free (first);

  compare (0);
  compare (1);

  /* off: the tree stays pruned; an update makes its own nodes, no crash */
  int npruned = tree_pruned ();
  tree_mpi_end ();
  CHECK (!tree_mpi_plan () && !tree_get_partition ()
         && tree_pruned () == npruned);
  ghosts_update_list (nlist, list, 1);

  free (nodes), free (alive);
  tree_destroy ();
  int f = fails, ftot;
  MPI_Allreduce (&f, &ftot, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  if (rank == 0)
    printf ("dimension %d, %d processes: %d checks, %d failed\n", dimension,
            npe, checks, ftot);
  MPI_Finalize ();
  return ftot != 0;
}
