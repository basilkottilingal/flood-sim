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
#include "grid-scalar.h"

scalar a;
face_vector F;
vertex_scalar v;
scalar d;                        /* moved with its blocks: milestone 3 */

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

/* the block of scalar s in node n, NULL if none: an inactive family's */
static Real * blk (BAMRNode * n, int s)
{
  void ** g = n->s ? n->s [s / SCALAR_ARRAY_WIDTH] : NULL;
  return g ? g [s % SCALAR_ARRAY_WIDTH] : NULL;
}

/* every cell of every block of the scalars, set to value */
static void prefill (Real value)
{
  for (int q = 0; q < nnodes; q++)
    if (alive [q] && blk (nodes [q].node, list [0]))
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
    if (!alive [q] || !blk (nodes [q].node, list [0]))
      continue;                     /* no blocks: inactive */
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

/* a weight for leaf p, from its key: 7 or 8 on the second half of the
 * curve, 1 or 2 on the first, so that the cut moves far; integers, so that
 * sums are exact in any order */
static double weight_of (Index p)
{
  uint64_t key = tree_morton_key (p);
  int second = (key >> (dimension * BAMR_MAX_DEPTH - 1)) & 1;
  return 1 + 6 * second + (double) ((key * 2654435761u >> 7) % 2);
}

/*
.. an oracle, from the whole tree, for weights w (of the whole tree's
.. leaves, NULL: 1 each): the pids, the partition (tree_partition ()),
.. the plan (tree_mpi_plan_new ()), what is needed (its own and what it
.. receives: active), the skeleton (its own and N () of them: the
.. neighbours, the coarse leaf covering an absent one, a leaf's touching
.. children of a neighbour with children, the parent, the children), and
.. so what a process keeps: the root, and every node whose parent has
.. something of the skeleton strictly below; a node has children iff
.. something of the skeleton is below it
*/
typedef struct {
  int * pid, * first, * oflags, nwant, nkeep, depth [2];
  uint64_t * keys;
  double imbalance, * load;
  Entry * want;
  char * need, * below, * keep, * sk;
} Oracle;

static int geometry = 1;           /* the skeleton holds what is read */

/* N (p) on the whole tree, marked in sk */
static void skeleton_of (Index p, char * sk)
{
  BAMRNode * node = index_node (p);
  int nb = 1 << p.level;
  for (int dir = 0; dir < TREE_NDIRS; dir++) {
    Index b = p;
    int in = 1, t = dir;
    for (int a = 0; a < dimension; a++, t /= 3) {
      int d = t % 3 - 1;
      if (a == 0) b.i += d;
#if dimension > 1
      if (a == 1) b.j += d;
#endif
#if dimension > 2
      if (a == 2) b.k += d;
#endif
      int x = index_coord (b, a);
      in &= x >= 0 && x < nb;
    }
    if (dir == TREE_DIR_CENTRE || !in)
      continue;
    if (!index_is_allocated (b)) {
      sk [position (index_parent (b))] = 1;
      continue;
    }
    sk [position (b)] = 1;
    if ((node->flags & BAMR_NODE_LEAF)
        && (index_node (b)->flags & BAMR_NODE_HAS_CHILDREN))
      for (int c = 0; c < TREE_NCHILDREN; c++) {
        int touch = 1, u = dir;
        for (int a = 0; a < dimension; a++, u /= 3) {
          int d = u % 3 - 1;
          touch &= !d || ((c >> a) & 1) == (d < 0);
        }
        if (touch)
          sk [position (index_child (b, c))] = 1;
      }
  }
  if (p.level > 0)
    sk [position (index_parent (p))] = 1;
  if (node->flags & BAMR_NODE_HAS_CHILDREN)
    for (int c = 0; c < TREE_NCHILDREN; c++)
      sk [position (index_child (p, c))] = 1;
}

static Oracle oracle_of (const double * w)
{
  Oracle o;
  CHECK (tree_partition (npe, w) == 0);
  const TreePartition * tp = tree_get_partition ();
  o.pid = malloc ((size_t) nnodes * sizeof (int));
  o.first = malloc ((size_t) (npe + 1) * sizeof (int));
  o.keys = malloc ((size_t) (npe + 1) * sizeof (uint64_t));
  o.load = malloc ((size_t) npe * sizeof (double));
  for (int q = 0; q < nnodes; q++)
    o.pid [q] = nodes [q].node->pid;
  memcpy (o.first, tp->leaf_first, (size_t) (npe + 1) * sizeof (int));
  memcpy (o.keys, tp->key_first, (size_t) (npe + 1) * sizeof (uint64_t));
  memcpy (o.load, tp->load, (size_t) npe * sizeof (double));
  o.imbalance = tp->imbalance;
  TreeMpiPlan * plan = tree_mpi_plan_new (rank);
  o.depth [0] = plan->depth, o.depth [1] = plan->depth_min;
  o.nwant = entries_of (plan, &o.want);
  tree_mpi_plan_free (plan);
  o.need = calloc ((size_t) nnodes, 1);
  o.below = calloc ((size_t) nnodes, 1);
  o.keep = calloc ((size_t) nnodes, 1);
  o.oflags = malloc ((size_t) nnodes * sizeof (int));
  o.nkeep = 0;
  for (int q = 0; q < nnodes; q++)
    o.need [q] = o.pid [q] == rank, o.oflags [q] = nodes [q].node->flags;
  for (int k = 0; k < o.nwant; k++)
    if (!o.want [k].send)
      o.need [position (o.want [k].p)] = 1;
  char * sk = o.sk = calloc ((size_t) nnodes, 1);
  for (int q = 0; q < nnodes; q++)
    if (o.pid [q] == rank) {
      sk [q] = 1;
      skeleton_of (nodes [q].p, sk);
    }
  for (int q = 0; q < nnodes; q++) {
    geometry &= !o.need [q] || sk [q];   /* what is read is in the skeleton */
    for (Index p = nodes [q].p; sk [q] && p.level > 0; ) {
      p = index_parent (p);
      o.below [position (p)] = 1;
    }
  }

  for (int q = 0; q < nnodes; q++) {
    o.keep [q] = !nodes [q].p.level
                 || o.below [position (index_parent (nodes [q].p))];
    o.nkeep += o.keep [q];
  }
  tree_partition_free ();
  return o;
}

static void oracle_free (Oracle * o)
{
  free (o->pid), free (o->first), free (o->keys), free (o->load);
  free (o->want), free (o->need), free (o->below), free (o->keep);
  free (o->oflags), free (o->sk);
}

/*
.. the distributed tree against an oracle, made from the queries, not the
.. whole tree: the plan the oracle's, entry for entry; the nodes kept, by
.. index, the oracle's; their pids; a needed node a leaf, a parent or
.. pruned as there; the counts; the partition and its keys; the leaves of
.. this process. The nodes are found again by index (nodes [], alive []):
.. those freed and made again are others
*/
static void check_oracle (const Oracle * o, const char * what)
{
  const TreeMpiPlan * plan = tree_mpi_plan ();
  CHECK (npe == 1 || (plan && plan->rank == rank && plan->npe == npe));
  if (!plan)
    return;
  Entry * got;
  int ngot = entries_of (plan, &got), equal = ngot == o->nwant;
  for (int k = 0; k < ngot && equal; k++)
    equal = got [k].parents == o->want [k].parents
            && got [k].kind == o->want [k].kind
            && got [k].level == o->want [k].level
            && got [k].send == o->want [k].send
            && got [k].peer == o->want [k].peer
            && same_index (got [k].p, o->want [k].p);
  CHECK (equal && plan->depth == o->depth [0]
         && plan->depth_min == o->depth [1]);
  const TreePartition * tp = tree_get_partition ();
  int kept = 1, pids = 1, flags = 1, blocks = 1, nalive = 0;
  for (int q = 0; q < nnodes; q++) {
    BAMRNode * here = index_is_allocated (nodes [q].p)
                      ? index_node (nodes [q].p) : NULL;
    kept &= (here != NULL) == o->keep [q];
    pids &= !here || here->pid == o->pid [q];
    if (here) {
      int f = here->flags & (BAMR_NODE_LEAF | BAMR_NODE_HAS_CHILDREN
                             | BAMR_NODE_IS_PRUNED);
      int w = (o->oflags [q] & BAMR_NODE_LEAF) ? BAMR_NODE_LEAF
              : o->below [q] ? BAMR_NODE_HAS_CHILDREN : BAMR_NODE_IS_PRUNED;
      /* every node of the skeleton as on the whole tree (a leaf, a
       * parent, or pruned); active iff needed; blocks iff its family
       * has an active member */
      /* a node of the skeleton as on the whole tree; another (a sibling
       * a family brings) may be of a structure unknown here: pruned */
      flags &= f == w || (!o->sk [q] && f == BAMR_NODE_IS_PRUNED);
      flags &= (!(here->flags & BAMR_NODE_IS_INACTIVE)) == o->need [q];
      int fam = 0;
      if (!nodes [q].p.level)
        fam = o->need [q];
      else
        for (int c = 0; c < TREE_NCHILDREN; c++)
          fam |= o->need [position (index_child (index_parent (nodes [q].p), c))];
      blocks &= (blk (here, list [0]) != NULL) == fam;
    }
    nodes [q].node = here;
    alive [q] = here != NULL;
    nalive += alive [q];
  }
  CHECK (kept && pids && flags && blocks && geometry);
  CHECK (grid_scalar_check ());
  {
    /* a scalar made now: blocks in the active families only */
    scalar t;
    (void) t;
    CHECK (grid_scalar_check ());
  }
  CHECK (grid_scalar_check ());
  /* the peers both ways: a 0/1 row of npe sent to all (a test only) */
  int * row = calloc ((size_t) npe, sizeof (int));
  int * col = calloc ((size_t) npe, sizeof (int));
  for (int k = 0; k < plan->npeers; k++)
    row [plan->peers [k]] = 1;
  MPI_Alltoall (row, 1, MPI_INT, col, 1, MPI_INT, MPI_COMM_WORLD);
  int symmetric = !row [rank];
  for (int r = 0; r < npe; r++)
    symmetric &= row [r] == col [r];
  /* and every rank of the plan's lists a peer */
  for (int k = 0; k < ngot; k++)
    symmetric &= row [got [k].peer];
  CHECK (symmetric);
  free (row), free (col), free (got);
  /* the counts the tree keeps, against its flags */
  const TreeCaches * tc = tree_caches ();
  int npr = 0, nlf = 0;
  for (int q = 0; q < tc->nall; q++) {
    npr += (tc->all [q].node->flags & BAMR_NODE_IS_PRUNED) != 0;
    nlf += (tc->all [q].node->flags & BAMR_NODE_LEAF) != 0;
  }
  CHECK (npr == tree_pruned () && nlf == tree_get_nleaves ()
         && nlf == tc->nleaves);
  CHECK (tc->nall == o->nkeep && tree_get_nblocks () == o->nkeep);
  /* the partition rank 0's (or the prefix sums'), as the oracle's; the
   * owner of every node of the whole tree, freed here or not, from the
   * first keys alone */
  int loads = 1;
  for (int r = 0; r < npe; r++)
    loads &= fabs (tp->load [r] - o->load [r]) <= 1e-12 * o->load [r];
  CHECK (tp && tp->npe == npe && loads
         && fabs (tp->imbalance - o->imbalance) <= 1e-12 * o->imbalance
         && !memcmp (tp->leaf_first, o->first, (size_t) (npe + 1)
                                              * sizeof (int))
         && !memcmp (tp->key_first, o->keys, (size_t) (npe + 1)
                                             * sizeof (uint64_t)));
  int owners = 1;
  for (int q = 0; q < nnodes; q++)
    owners &= tree_partition_owner (nodes [q].p) == o->pid [q];
  CHECK (owners);
  /* the leaves of this process: those kept with its rank, in order;
   * over the processes every leaf once */
  int listed = 1, k = 0;
  for (int q = 0; q < nnodes; q++)
    if (o->pid [q] == rank && alive [q]
        && (nodes [q].node->flags & BAMR_NODE_LEAF)) {
      listed &= k < tc->nowned && tc->owned [k].node == nodes [q].node;
      k++;
    }
  listed &= k == tc->nowned;
  long nblk = 0;
  for (int q = 0; q < nnodes; q++)
    nblk += alive [q] && blk (nodes [q].node, list [0]) != NULL;
  long counts [3] = { tc->nowned, nalive, nblk }, sums [3];
  MPI_Allreduce (counts, sums, 3, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
  CHECK (listed && sums [0] == whole_leaves);
  if (rank == 0)
    printf ("  %s: a process keeps %.2f of the %d nodes, %.2f with blocks; "
            "imbalance %.3f\n", what, (double) sums [1] / npe / nnodes,
            nnodes, (double) sums [2] / npe / nnodes, tp->imbalance);
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

  /* the oracles, from the whole tree: unit weights, then weights by
   * index, the whole tree's leaves in order */
  double * w = malloc ((size_t) whole_leaves * sizeof (double));
  for (int q = 0; q < whole_leaves; q++)
    w [q] = weight_of (tree_caches ()->leaves [q].p);
  Oracle unit = oracle_of (NULL), weighted = oracle_of (w);
  free (w);
  for (int q = 0; q < nnodes; q++)
    nodes [q].node->pid = rank;   /* this process's again: rank 0 tells */
  tree_caches_changed ();

  CHECK (tree_mpi_partition (NULL) == 0);
  check_oracle (&unit, "partition");

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

  /*
  .. repartition on the distributed tree (milestone 3). A field d, filled
  .. on this process's leaves, every other block of it poisoned (NaN):
  .. after a repartition, a leaf that came here without its data shows.
  .. Unit weights again: the same cut, nothing sent, the same plan. Then
  .. weights by index: the oracle's cut, plan and nodes; d exact on every
  .. leaf of this process; the ghosts again bit for bit
  */
  for (int l = 0; l <= BAMR_MAX_DEPTH; l++)     /* every node: parents too */
    foreach_level (l)
      d[] = sin (5 * x + 0.3) + cos (3 * y) * z;
  for (int q = 0; q < nnodes; q++)
    if (alive [q] && blk (nodes [q].node, d)) {
      BAMRNode * node = nodes [q].node;
      int mine = node->pid == rank;
      Real * b = block_data (node->s, d);
      for (long c = 0; c < cells; c++)
        if (!mine || !owned (tree_block_shape (), CWRAP_CENTRE, c))
          b [c] = (Real) NAN;
    }
  long moved = tree_mpi_repartition (NULL), allmoved;
  MPI_Allreduce (&moved, &allmoved, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
  CHECK (allmoved == 0);
  check_oracle (&unit, "unit weights again");
  {
    const TreeCaches * tc = tree_caches ();
    int nw = tc->nowned;
    double * wo = malloc ((size_t) (nw ? nw : 1) * sizeof (double));
    for (int q = 0; q < nw; q++)
      wo [q] = weight_of (tc->owned [q].p);
    moved = tree_mpi_repartition (wo);
    free (wo);
  }
  MPI_Allreduce (&moved, &allmoved, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
  CHECK (npe == 1 ? allmoved == 0 : allmoved > 0);
  check_oracle (&weighted, "weights");
  {
    /* d on every node of this process, parents too, as filled */
    long bad = 0, seen = 0, all [2], mine [2];
    for (int l = 0; l <= BAMR_MAX_DEPTH; l++)
      foreach_level (l, reduction (+:bad, seen)) {
        bad += !(d[] == (Real) (sin (5 * x + 0.3) + cos (3 * y) * z));
        seen++;
      }
    /* every other scalar too, on the leaves: the values the last update
     * left, which fill () and the same update make again (a refluxed face
     * is the mean of the fine faces beside it) */
    BlockShape sh = tree_block_shape ();
    const TreeCaches * tc = tree_caches ();
    Real * kept = malloc ((size_t) (tc->nowned ? tc->nowned : 1) * nlist
                          * cells * sizeof (Real));
    for (int q = 0; q < tc->nowned; q++)
      for (int k = 0; k < nlist; k++)
        memcpy (kept + ((size_t) q * nlist + k) * cells,
                block_data (tc->owned [q].node->s, list [k]),
                cells * sizeof (Real));
    fill ();
    ghosts_update_list (nlist, list, 1);
    tc = tree_caches ();
    for (int q = 0; q < tc->nowned; q++)
      for (int k = 0; k < nlist; k++) {
        const Real * b = block_data (tc->owned [q].node->s, list [k]);
        const Real * o = kept + ((size_t) q * nlist + k) * cells;
        for (long c = 0; c < cells; c++)
          if (owned (sh, scalar_location (list [k]), c))
            bad += memcmp (&b [c], &o [c], sizeof (Real)) != 0, seen++;
      }
    free (kept);
    mine [0] = bad, mine [1] = seen;
    MPI_Allreduce (mine, all, 2, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
    CHECK (all [0] == 0 && all [1] > 0);
    if (rank == 0)
      printf ("  %ld nodes moved; every scalar on every node of its owner "
              "as before: %ld of %ld cells wrong\n", allmoved, all [0],
              all [1]);
  }
  compare (0);
  compare (1);
  oracle_free (&unit), oracle_free (&weighted);
  for (int p = 0; p < 2; p++)
    free (r1 [p]), free (r2 [p]);

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
