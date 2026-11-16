/* tree-mpi.h with MPI: the ghosts of a cell scalar, a face vector and a
 * vertex scalar, and the restriction of the parents, made on npe processes
 * by the plan's exchanges on the pruned trees, against the same made by
 * one process on the whole tree, bit for bit. Every block a process keeps
 * but does not own is poisoned (NaN) first, and its own leaves' ghosts
 * too: a block missing from the plan, or sent too early, leaves a NaN or a
 * different value; one pruned that is read, a crash (ASan). Leaves and
 * every node (parents too). The partition, made by rank 0 and broadcast,
 * against the same made by each process; every node's owner from the
 * first keys alone; each process's leaves, every leaf once over them;
 * reductions over the processes, the whole tree once. MPI started by
 * tree (). Build with -D_MPI, run with mpiexec -n 1, 2, 3, 4 ... */
#include <stdint.h>
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

/* a plan's lists as entries, in its order: send, then receive, by
 * parents, kind, level, peer */
typedef struct { int parents, kind, level, send, peer; Index p; } Entry;

static int entries_of (const TreeMpiPlan * plan, Entry ** e)
{
  int n = 0, cap = 0;
  *e = NULL;
  for (int pa = 0; pa < 2; pa++)
    for (int kind = 0; kind < TREE_MPI_KINDS; kind++)
      for (int l = 0; l < BAMR_MAX_DEPTH + 2; l++)
        for (int send = 0; send < 2; send++) {
          const TreeMpiExchange * ex = &plan->ex [pa][kind][l];
          int nl = send ? ex->nsend : ex->nrecv;
          const TreeMpiList * list = send ? ex->send : ex->recv;
          for (int k = 0; k < nl; k++)
            for (int j = 0; j < list [k].n; j++) {
              if (n == cap)
                cap = cap ? 2 * cap : 256,
                *e = realloc (*e, (size_t) cap * sizeof (Entry));
              (*e) [n++] = (Entry) { pa, kind, l, send, list [k].peer,
                                     list [k].item [j].p };
            }
        }
  return n;
}

static int same_index (Index a, Index b)
{
  int same = a.level == b.level && a.i == b.i;
#if dimension > 1
  same &= a.j == b.j;
#endif
#if dimension > 2
  same &= a.k == b.k;
#endif
  return same;
}

/* the node q of nodes [] at p, by a search; -1 if none */
static int position (Index p)
{
  for (int q = 0; q < nnodes; q++)
    if (same_index (nodes [q].p, p))
      return q;
  return -1;
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

/* stencil loops with reduction (): every operator, a loop by axis, a loop
 * over a level; sums start at START, which must count once */
enum { START = 7 };
struct sums { long cells, faces, levels, level2, prod; int lmin, lmax;
              double xmin, xmax; int rmin, rmax; };
static struct sums whole;
static long whole_leaves;

static struct sums reduced (void)
{
  long cells = START, faces = START, levels = 0, level2 = 0, prod = 1;
  int lmin = 1 << 20, lmax = -1;
  double xmin = 2, xmax = -1;
  int rmin = 1 << 20, rmax = -1;      /* the ranks with cells: min and max
                                       * differ from process to process */
  foreach (reduction (+:cells, levels) reduction (*:prod)
           reduction (min:lmin, xmin, rmin) reduction (max:lmax, xmax, rmax)) {
    int l = (int) lround (log2 (1. / (Delta * tree_get_block_m ())));
    cells++;
    levels += l;
    prod *= 1 + 2 * (x < Delta && y < Delta && z < Delta);  /* one cell */
    lmin = l < lmin ? l : lmin;
    lmax = l > lmax ? l : lmax;
    xmin = x < xmin ? x : xmin;
    xmax = x > xmax ? x : xmax;
    int r = tree_rank ();
    rmin = r < rmin ? r : rmin;
    rmax = r > rmax ? r : rmax;
  }
  foreach_face (reduction (+:faces))
    faces++;
  foreach_level (2, reduction (+:level2))
    level2++;
  return (struct sums) { cells, faces, levels, level2, prod, lmin, lmax,
                         xmin, xmax, rmin, rmax };
}

int main (void)
{
  /* MPI started by tree_replicated () (as by tree ()), and finalised at
   * exit; the whole tree on every process, refined, then distributed */
  int on = 0;
  MPI_Initialized (&on);
  tree_replicated (5, 4);
  int now = 0;
  MPI_Initialized (&now);
  MPI_Comm_rank (MPI_COMM_WORLD, &rank);
  MPI_Comm_size (MPI_COMM_WORLD, &npe);
  CHECK (!on && now);
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
  whole = reduced ();
  whole_leaves = tree_caches ()->nleaves;

  /* the partition each process would make, to compare with rank 0's */
  int * pid = malloc ((size_t) nnodes * sizeof (int));
  int * first = malloc ((size_t) (npe + 1) * sizeof (int));
  CHECK (tree_partition (npe, NULL) == 0);
  for (int q = 0; q < nnodes; q++)
    pid [q] = nodes [q].node->pid;
  memcpy (first, tree_get_partition ()->leaf_first,
          (size_t) (npe + 1) * sizeof (int));
  double imbalance = tree_get_partition ()->imbalance;
  uint64_t * keys = malloc ((size_t) (npe + 1) * sizeof (uint64_t));
  memcpy (keys, tree_get_partition ()->key_first,
          (size_t) (npe + 1) * sizeof (uint64_t));

  /*
  .. the oracle, from the whole tree: the plan (tree_mpi_plan_new ()),
  .. what is needed (its own and what it receives), and so what a process
  .. keeps: the root, and every node whose parent has something needed
  .. strictly below; a needed node has children iff something needed is
  .. below it
  */
  TreeMpiPlan * oracle = tree_mpi_plan_new (rank);
  int oracle_depth [2] = { oracle->depth, oracle->depth_min };
  Entry * want;
  int nwant = entries_of (oracle, &want);
  char * need = calloc ((size_t) nnodes, 1);
  char * below = calloc ((size_t) nnodes, 1);
  char * keep = calloc ((size_t) nnodes, 1);
  int * oflags = malloc ((size_t) nnodes * sizeof (int)), nkeep = 0;
  for (int q = 0; q < nnodes; q++)
    need [q] = pid [q] == rank, oflags [q] = nodes [q].node->flags;
  for (int k = 0; k < nwant; k++)
    if (!want [k].send)
      need [position (want [k].p)] = 1;
  for (int q = 0; q < nnodes; q++)
    for (Index p = nodes [q].p; need [q] && p.level > 0; ) {
      p = index_parent (p);
      below [position (p)] = 1;
    }
  for (int q = 0; q < nnodes; q++) {
    keep [q] = !nodes [q].p.level
              || below [position (index_parent (nodes [q].p))];
    nkeep += keep [q];
  }
  tree_mpi_plan_free (oracle);
  tree_partition_free ();
  for (int q = 0; q < nnodes; q++)
    nodes [q].node->pid = rank;   /* this process's again: rank 0 tells */
  tree_caches_changed ();

  CHECK (tree_mpi_partition (NULL) == 0);
  const TreeMpiPlan * plan = tree_mpi_plan ();
  CHECK (plan && plan->rank == rank && plan->npe == npe);
  CHECK (npe > 1 ? tree_pruned () > 0 : !tree_pruned ());
  {
    /*
    .. made from the queries, not the whole tree: the plan the oracle's,
    .. entry for entry; the nodes kept, by index, the oracle's; their
    .. pids; a needed node a leaf, a parent or pruned as there. The nodes
    .. are found again by index: those freed and made again are others
    */
    Entry * got;
    int ngot = entries_of (plan, &got), equal = ngot == nwant;
    for (int k = 0; k < ngot && equal; k++)
      equal = got [k].parents == want [k].parents && got [k].kind == want [k].kind
              && got [k].level == want [k].level && got [k].send == want [k].send
              && got [k].peer == want [k].peer
              && same_index (got [k].p, want [k].p);
    CHECK (equal && plan->depth == oracle_depth [0]
           && plan->depth_min == oracle_depth [1]);
    const TreePartition * tp = tree_get_partition ();
    int kept = 1, pids = 1, flags = 1, nalive = 0;
    for (int q = 0; q < nnodes; q++) {
      BAMRNode * here = index_is_allocated (nodes [q].p)
                        ? index_node (nodes [q].p) : NULL;
      kept &= (here != NULL) == keep [q];
      pids &= !here || here->pid == pid [q];
      if (here && need [q]) {
        int f = here->flags & (BAMR_NODE_LEAF | BAMR_NODE_HAS_CHILDREN
                               | BAMR_NODE_IS_PRUNED);
        int w = (oflags [q] & BAMR_NODE_LEAF) ? BAMR_NODE_LEAF
                : below [q] ? BAMR_NODE_HAS_CHILDREN : BAMR_NODE_IS_PRUNED;
        flags &= f == w;
      }
      nodes [q].node = here;
      alive [q] = here != NULL;
      nalive += alive [q];
    }
    CHECK (kept && pids && flags);
    /* the counts the tree keeps, against its flags */
    const TreeCaches * tcc = tree_caches ();
    int npr = 0, nlf = 0;
    for (int q = 0; q < tcc->nall; q++) {
      npr += (tcc->all [q].node->flags & BAMR_NODE_IS_PRUNED) != 0;
      nlf += (tcc->all [q].node->flags & BAMR_NODE_LEAF) != 0;
    }
    CHECK (npr == tree_pruned () && nlf == tree_get_nleaves ()
           && nlf == tcc->nleaves);
    CHECK (tp && tp->npe == npe && tp->imbalance == imbalance
           && !memcmp (tp->leaf_first, first, (size_t) (npe + 1)
                                                * sizeof (int)));
    CHECK (tree_caches ()->nall == nkeep && tree_get_nblocks () == nkeep);
    free (got);
    /* the owner of every node of the whole tree, freed here or not, from
     * the first keys alone; the keys rank 0's, as made here */
    int owners = !memcmp (tp->key_first, keys, (size_t) (npe + 1)
                                                * sizeof (uint64_t));
    for (int q = 0; q < nnodes; q++)
      owners &= tree_partition_owner (nodes [q].p) == pid [q];
    CHECK (owners);
    /* the leaves of this process: those kept with its rank, in order;
     * over the processes every leaf once */
    const TreeCaches * tc = tree_caches ();
    int listed = 1, k = 0;
    for (int q = 0; q < nnodes; q++)
      if (pid [q] == rank && alive [q]
          && (nodes [q].node->flags & BAMR_NODE_LEAF)) {
        listed &= k < tc->nowned && tc->owned [k].node == nodes [q].node;
        k++;
      }
    listed &= k == tc->nowned;
    long counts [2] = { tc->nowned, nalive }, sums [2];
    MPI_Allreduce (counts, sums, 2, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
    CHECK (listed && sums [0] == whole_leaves);
    if (rank == 0)
      printf ("  %.2f of the %d nodes kept by a process\n",
              (double) sums [1] / npe / nnodes, nnodes);
  }
  free (pid), free (first), free (keys);
  free (want), free (need), free (below), free (keep), free (oflags);

  /* reductions over the processes, each its own leaves now: the whole
   * tree once. Before the partition each process owned all of it: npe
   * times (the value before a loop counted once, rank 0's) */
  {
    struct sums got = reduced ();
    /* the ranks with cells: from 0 to the last one with a leaf */
    const TreePartition * tp = tree_get_partition ();
    int last = npe - 1;
    while (last > 0 && tp->leaf_first [last] == tp->leaf_first [last + 1])
      last--;
    CHECK (got.rmin == 0 && got.rmax == last && whole.rmin == 0
           && whole.rmax == npe - 1);
    long p = 1;
    for (int r = 0; r < npe; r++)
      p *= got.prod;
    CHECK (got.cells - START == (whole.cells - START) / npe
           && (whole.cells - START) % npe == 0
           && got.faces - START == (whole.faces - START) / npe
           && got.level2 * npe == whole.level2 && got.levels * npe == whole.levels
           && got.prod == 3 && whole.prod == p
           && got.lmin == whole.lmin && got.lmax == whole.lmax
           && got.lmin < got.lmax
           && got.xmin == whole.xmin && got.xmax == whole.xmax
           && got.xmin < 0.05 && got.xmax > 0.95);
    if (rank == 0)
      printf ("  reductions: %ld cells, %ld faces, levels %d to %d, "
              "corner %ld, level-2 cells %ld: the whole tree once\n",
              got.cells - START, got.faces - START, got.lmin, got.lmax,
              got.prod, got.level2);
  }

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
  return ftot != 0;                   /* MPI_Finalize () at exit: tree ()'s */
}
