#include <assert.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef _MPI
# include <mpi.h>
#endif

#include "block-ghosts.h"
#include "tree-mpi.h"

/* a node needed by the process being pruned for, and the root of a
.. subtree with nothing needed: set only during tree_mpi_prune (), bits no
.. other code uses */
#define TREE_MPI_NEEDED (1 << 30)
#define TREE_MPI_EMPTY  (1 << 29)

/* ---------------------------------------------------------------------- */
/* 1. the partition                                                       */
/* ---------------------------------------------------------------------- */

static TreePartition part;
static int made = 0;

void tree_partition_free (void)
{
  free (part.leaf_first);
  free (part.load);
  free (part.key_first);
  part = (TreePartition) { 0 };
  made = 0;
}

int tree_partition_owner (Index p)
{
  if (!made)
    return tree_rank ();              /* no partition: the tree is ours */
  /* the last process whose first key is at most p's: an empty one has the
   * next one's key, so it is passed over */
  uint64_t key = tree_morton_key (p);
  int lo = 0, hi = part.npe;          /* key_first [lo] <= key < [hi] */
  while (hi - lo > 1)
  {
    int mid = (lo + hi) / 2;
    if (part.key_first [mid] <= key)
      lo = mid;
    else
      hi = mid;
  }
  return lo;
}

const TreePartition * tree_get_partition (void)
{
  return made ? &part : NULL;
}

int tree_partition_valid (void)
{
  return made && part.generation == tree_generation;
}

/*
.. the partition recorded from the nodes' pids, set already, and each
.. process's load (taken). The owners changed: the lists again (the
.. leaves of this process)
*/
static void partition_record (int npe, double * load)
{
  tree_caches_changed ();
  const TreeCaches * tc = tree_caches ();
  int n = tc->nleaves;
  int * first = malloc ((size_t) (npe + 1) * sizeof (int));
  uint64_t * key = malloc ((size_t) (npe + 1) * sizeof (uint64_t));
  if (!first || !key)
    abort ();
  /* the first leaf of each process: pid does not decrease along the curve;
   * its key, or the next process's for one with no leaf */
  for (int r = 0, q = 0; r <= npe; r++)
  {
    while (q < n && tc->leaves [q].node->pid < r)
      q++;
    first [r] = q;
  }
  for (int r = npe; r >= 0; r--)
    key [r] = r == npe ? UINT64_MAX
            : first [r] < first [r + 1]
            ? tree_morton_key (tc->leaves [first [r]].p)
            : key [r + 1];
  key [0] = 0;                        /* the curve starts with process 0 */
  double total = 0, most = 0;
  for (int r = 0; r < npe; r++)
  {
    total += load [r];
    most = load [r] > most ? load [r] : most;
  }
  tree_partition_free ();
  part = (TreePartition) {
    .npe = npe, .nleaves = n, .leaf_first = first, .load = load,
    .key_first = key,
    .imbalance = n ? most / (total / npe) : 1, .generation = tree_generation,
  };
  made = 1;
}

int tree_partition (int npe, const double * weight)
{
  if (npe < 1 || tree_pruned ())
    return -1;
  const TreeCaches * tc = tree_caches ();
  int n = tc->nleaves;
  double total = 0;
  if (weight)
    for (int q = 0; q < n; q++)
    {
      if (!(weight [q] > 0))        /* NaN too */
        return -1;
      total += weight [q];
    }
  double * load = calloc ((size_t) npe, sizeof (double));
  if (!load)
    abort ();

  /*
  .. the leaves: leaf q to the process whose share of the curve holds its
  .. middle; in integers for unit weights, so the pieces are exact
  */
  double before = 0;
  for (int q = 0; q < n; q++)
  {
    int pid;
    double w = weight ? weight [q] : 1;
    if (weight)
    {
      pid = (int) floor (npe * (before + w / 2) / total);
      before += w;
    }
    else
      pid = (int) ((long long) npe * (2 * (long long) q + 1)
                   / (2 * (long long) n));
    pid = pid < 0 ? 0 : pid >= npe ? npe - 1 : pid;
    tc->leaves [q].node->pid = pid;
    load [pid] += w;
  }
  /* the parents, finest first: the smallest pid of their children */
  for (int l = tc->depth - 1; l >= 0; l--)
    for (int b = tc->parent_first [l]; b < tc->parent_first [l + 1]; b++)
    {
      Index p = tc->parents_by_level [b].p;
      int pid = npe;
      for (int c = 0; c < TREE_NCHILDREN; c++)
      {
        int cp = index_node (index_child (p, c))->pid;
        pid = cp < pid ? cp : pid;
      }
      tc->parents_by_level [b].node->pid = pid;
    }

  partition_record (npe, load);
  return 0;
}

/* ---------------------------------------------------------------------- */
/* 2. the neighbours                                                      */
/* ---------------------------------------------------------------------- */

static int index_axis (Index p, int a)
{
  return a == 0 ? p.i
  #if dimension > 1
       : a == 1 ? p.j
  #endif
  #if dimension > 2
       : a == 2 ? p.k
  #endif
       : 0;
}

static Index index_shift (Index p, int a, int d)
{
  if (a == 0)
    p.i += d;
  #if dimension > 1
  if (a == 1)
    p.j += d;
  #endif
  #if dimension > 2
  if (a == 2)
    p.k += d;
  #endif
  return p;
}

/* p + d for direction dir; 0 if outside the domain */
static int index_dir (Index p, int dir, Index * b)
{
  int nb = 1 << p.level;
  *b = p;
  for (int a = 0; a < dimension; a++, dir /= 3)
  {
    *b = index_shift (*b, a, dir % 3 - 1);
    int x = index_axis (*b, a);
    if (x < 0 || x >= nb)
      return 0;
  }
  return 1;
}

TreeNeighbours tree_neighbours (Index p)
{
  TreeNeighbours m = { 0, 0, 0, { 0 } };
  for (int dir = 0; dir < TREE_NDIRS; dir++)
  {
    Index b;
    m.pid [dir] = -1;
    if (dir == TREE_DIR_CENTRE || !index_dir (p, dir, &b))
      continue;
    uint32_t bit = (uint32_t) 1 << dir;
    if (index_is_allocated (b))
    {
      BAMRNode * node = index_node (b);
      if (node->flags & BAMR_NODE_HAS_CHILDREN)
        m.fine |= bit;
      else
        m.same |= bit;
      m.pid [dir] = node->pid;
    }
    else
    {
      Index C = index_parent (b);
      assert (p.level > 0 && index_is_allocated (C));
      m.coarse |= bit;
      m.pid [dir] = index_node (C)->pid;
    }
  }
  return m;
}

/* ---------------------------------------------------------------------- */
/* 3. the plan                                                            */
/* ---------------------------------------------------------------------- */

typedef struct
{
  int parents, kind, level, send, peer;
  TreeItem it;
} Need;

typedef struct
{
  Need * a;
  long n, cap;
  int local;              /* receives only: the sends are asked for */
} Needs;

/* by level, then k, j, i: the order of both sides */
static int index_cmp (Index a, Index b)
{
  if (a.level != b.level)
    return a.level < b.level ? -1 : 1;
  for (int ax = dimension - 1; ax >= 0; ax--)
  {
    int x = index_axis (a, ax), y = index_axis (b, ax);
    if (x != y)
      return x < y ? -1 : 1;
  }
  return 0;
}

static int need_cmp (const void * pa, const void * pb)
{
  const Need * a = pa, * b = pb;
  int ka [5] = { a->parents, a->kind, a->level, a->send, a->peer };
  int kb [5] = { b->parents, b->kind, b->level, b->send, b->peer };
  for (int q = 0; q < 5; q++)
    if (ka [q] != kb [q])
      return ka [q] < kb [q] ? -1 : 1;
  return index_cmp (a->it.p, b->it.p);
}

/* the node X (at x) read by a node of process r, for rank's plan */
static void need (Needs * v, int rank, int parents, int kind, int level,
                  int r, Index x)
{
  BAMRNode * X = index_node (x);
  assert (X);
  int owner = X->pid, send;
  if (owner == r)
    return;
  if (r == rank)
    send = 0;
  else if (owner == rank && !v->local)
    send = 1;
  else
    return;
  if (v->n == v->cap)
  {
    v->cap = v->cap ? 2 * v->cap : 1024;
    v->a = realloc (v->a, (size_t) v->cap * sizeof (Need));
    if (!v->a)
      abort ();
  }
  v->a [v->n++] = (Need) { parents, kind, level, send, send ? r : owner,
                           { x, X } };
}

/*
.. what rank receives (and sends, unless v->local) by the tree here: the
.. ghosts of its nodes, the reflux of its leaves, the restriction of its
.. parents down to level depth_min (every level for the parents' plan)
*/
static void plan_needs (Needs * v, int rank, int depth_min)
{
  const TreeCaches * tc = tree_caches ();
  for (int parents = 0; parents < 2; parents++)
  {
    /* the ghosts of every node updated: the leaves, or every node */
    for (int q = 0; q < tc->nall; q++)
    {
      Index p = tc->all [q].p;
      BAMRNode * node = tc->all [q].node;
      int leaf = node->flags & BAMR_NODE_LEAF, r = node->pid;
      if ((!parents && !leaf) || (v->local && r != rank))
        continue;
      TreeNeighbours m = tree_neighbours (p);
      for (int dir = 0; dir < TREE_NDIRS; dir++)
      {
        Index b;
        uint32_t bit = (uint32_t) 1 << dir;
        if (!((m.same | m.fine | m.coarse) & bit))
          continue;
        index_dir (p, dir, &b);
        need (v, rank, parents, TREE_MPI_GHOSTS, p.level, r,
              m.coarse & bit ? index_parent (b) : b);
      }
      /* the reflux of a leaf: the upper children along a of the parent
       * below it along a */
      if (leaf)
        for (int a = 0; a < dimension; a++)
        {
          int w = 1;
          for (int t = 0; t < a; t++)
            w *= 3;
          int dir = TREE_DIR_CENTRE - w;    /* d = -1 along a */
          if (!(m.fine & ((uint32_t) 1 << dir)))
            continue;
          Index lo = index_shift (p, a, -1);
          for (int c = 0; c < TREE_NCHILDREN; c++)
            if ((c >> a) & 1)
              need (v, rank, parents, TREE_MPI_FACES, p.level, r,
                    index_child (lo, c));
        }
    }
    /* the restriction: each parent, from its children */
    int lowest = parents ? 0 : depth_min;
    for (int l = lowest; l < tc->depth; l++)
      for (int b = tc->parent_first [l]; b < tc->parent_first [l + 1]; b++)
        if (!v->local || tc->parents_by_level [b].node->pid == rank)
          for (int c = 0; c < TREE_NCHILDREN; c++)
            need (v, rank, parents, TREE_MPI_RESTRICT, l,
                  tc->parents_by_level [b].node->pid,
                  index_child (tc->parents_by_level [b].p, c));
  }
}

/* v sorted, each once, then cut into plan's lists: an exchange, a peer
.. each, in the order both sides share; v freed */
static void plan_lists (TreeMpiPlan * plan, Needs * v)
{
  if (v->n)
    qsort (v->a, (size_t) v->n, sizeof (Need), need_cmp);
  long u = 0;
  for (long q = 0; q < v->n; q++)
    if (!u || need_cmp (&v->a [u - 1], &v->a [q]))
      v->a [u++] = v->a [q];
  for (long q = 0; q < u; )
  {
    long e = q;
    const Need * f = &v->a [q];
    while (e < u && v->a [e].parents == f->parents
           && v->a [e].kind == f->kind && v->a [e].level == f->level
           && v->a [e].send == f->send && v->a [e].peer == f->peer)
      e++;
    TreeMpiExchange * ex = &plan->ex [f->parents][f->kind][f->level];
    int * n = f->send ? &ex->nsend : &ex->nrecv;
    TreeMpiList ** list = f->send ? &ex->send : &ex->recv;
    *list = realloc (*list, (size_t) (*n + 1) * sizeof (TreeMpiList));
    TreeItem * item = malloc ((size_t) (e - q) * sizeof (TreeItem));
    if (!*list || !item)
      abort ();
    for (long k = q; k < e; k++)
      item [k - q] = v->a [k].it;
    (*list) [(*n)++] = (TreeMpiList) { f->peer, (int) (e - q), item };
    q = e;
  }
  free (v->a);
  *v = (Needs) { NULL, 0, 0, 0 };
}

TreeMpiPlan * tree_mpi_plan_new (int rank)
{
  assert (tree_partition_valid () && !tree_pruned ());
  assert (rank >= 0 && rank < part.npe);
  const TreeCaches * tc = tree_caches ();
  TreeMpiPlan * plan = calloc (1, sizeof (TreeMpiPlan));
  if (!plan)
    abort ();
  plan->rank = rank, plan->npe = part.npe;
  plan->depth = tc->depth, plan->depth_min = tree_get_depth_min ();
  plan->generation = tree_generation;
  Needs v = { NULL, 0, 0, 0 };
  plan_needs (&v, rank, plan->depth_min);
  plan_lists (plan, &v);
  return plan;
}

void tree_mpi_plan_free (TreeMpiPlan * plan)
{
  if (!plan)
    return;
  for (int parents = 0; parents < 2; parents++)
    for (int kind = 0; kind < TREE_MPI_KINDS; kind++)
      for (int l = 0; l < BAMR_MAX_DEPTH + 2; l++)
      {
        TreeMpiExchange * ex = &plan->ex [parents][kind][l];
        for (int k = 0; k < ex->nsend; k++)
          free (ex->send [k].item);
        for (int k = 0; k < ex->nrecv; k++)
          free (ex->recv [k].item);
        free (ex->send), free (ex->recv);
      }
  free (plan);
}

/* ---------------------------------------------------------------------- */
/* 4. pruning                                                             */
/* ---------------------------------------------------------------------- */

/* every node of process rank marked needed */
static void need_own (int rank)
{
  const TreeCaches * tc = tree_caches ();
  for (int q = 0; q < tc->nall; q++)
    if (tc->all [q].node->pid == rank)
      tc->all [q].node->flags |= TREE_MPI_NEEDED;
}

static int prune_unneeded (void);

int tree_mpi_prune (const TreeMpiPlan * plan)
{
  assert (plan->generation == tree_generation);
  /* needed: the process's nodes and what it receives */
  need_own (plan->rank);
  for (int parents = 0; parents < 2; parents++)
    for (int kind = 0; kind < TREE_MPI_KINDS; kind++)
      for (int l = 0; l < BAMR_MAX_DEPTH + 2; l++)
      {
        const TreeMpiExchange * ex = &plan->ex [parents][kind][l];
        for (int k = 0; k < ex->nrecv; k++)
          for (int q = 0; q < ex->recv [k].n; q++)
            ex->recv [k].item [q].node->flags |= TREE_MPI_NEEDED;
      }
  return prune_unneeded ();
}

/* the largest subtrees with no node marked needed below their root,
.. pruned; the marks cleared. The subtrees pruned */
static int prune_unneeded (void)
{
  const TreeCaches * tc = tree_caches ();
  /*
  .. in Morton order a parent comes right after its children's subtrees:
  .. any [l] says whether a node of level l since the last parent of level
  .. l - 1 has something needed in its subtree. A parent with nothing
  .. needed below is EMPTY
  */
  int any [BAMR_MAX_DEPTH + 2] = { 0 };
  for (int q = 0; q < tc->nall; q++)
  {
    BAMRNode * node = tc->all [q].node;
    int l = tc->all [q].p.level, below = 0;
    if (node->flags & BAMR_NODE_HAS_CHILDREN)
    {
      below = any [l + 1];
      any [l + 1] = 0;
      if (!below)
        node->flags |= TREE_MPI_EMPTY;
    }
    any [l] |= below || (node->flags & TREE_MPI_NEEDED);
  }
  /* the largest: EMPTY, its parent not */
  Index * root = malloc ((size_t) (tc->nall ? tc->nall : 1) * sizeof (Index));
  if (!root)
    abort ();
  int n = 0;
  for (int q = 0; q < tc->nall; q++)
  {
    Index p = tc->all [q].p;
    if ((tc->all [q].node->flags & TREE_MPI_EMPTY)
        && (!p.level || !(index_node (index_parent (p))->flags
                          & TREE_MPI_EMPTY)))
      root [n++] = p;
  }
  for (int q = 0; q < tc->nall; q++)
    tc->all [q].node->flags &= ~(TREE_MPI_NEEDED | TREE_MPI_EMPTY);
  for (int k = 0; k < n; k++)
    tree_prune (root [k]);
  free (root);
  return n;
}

/* ---------------------------------------------------------------------- */
/* 5. MPI                                                                 */
/* ---------------------------------------------------------------------- */

#ifdef _MPI

static TreeMpiPlan * applied = NULL;

const TreeMpiPlan * tree_mpi_plan (void)
{
  return applied;
}

void tree_mpi_end (void)
{
  tree_mpi_plan_free (applied);
  applied = NULL;
  tree_partition_free ();
}

/*
.. ints exchanged with every process: out [r][0 .. nout [r]) to r, into
.. in [r][0 .. nin [r]) from r (malloc'd, NULL for none); the counts
.. first (MPI_Alltoall). One's own is copied
*/
static void exchange_ints (int ** out, const int * nout, int ** in, int * nin,
                           int tag)
{
  int npe = tree_npe (), me = tree_rank (), nr = 0;
  MPI_Alltoall ((void *) nout, 1, MPI_INT, nin, 1, MPI_INT, MPI_COMM_WORLD);
  MPI_Request * req = malloc ((size_t) (2 * npe) * sizeof (MPI_Request));
  if (!req)
    abort ();
  for (int r = 0; r < npe; r++)
  {
    in [r] = NULL;
    if (!nin [r])
      continue;
    in [r] = malloc ((size_t) nin [r] * sizeof (int));
    if (!in [r])
      abort ();
    if (r == me)
      memcpy (in [r], out [r], (size_t) nin [r] * sizeof (int));
    else
      MPI_Irecv (in [r], nin [r] * (int) sizeof (int), MPI_BYTE, r, tag,
                 MPI_COMM_WORLD, &req [nr++]);
  }
  for (int r = 0; r < npe; r++)
    if (nout [r] && r != me)
      MPI_Isend (out [r], nout [r] * (int) sizeof (int), MPI_BYTE, r, tag,
                 MPI_COMM_WORLD, &req [nr++]);
  MPI_Waitall (nr, req, MPI_STATUSES_IGNORE);
  free (req);
}

/* an Index as 4 ints, and back */
static void put_index (int * a, Index p)
{
  a [0] = p.level, a [1] = p.i, a [2] = 0, a [3] = 0;
  #if dimension > 1
  a [2] = p.j;
  #endif
  #if dimension > 2
  a [3] = p.k;
  #endif
}

static Index get_index (const int * a)
{
  Index p;
  p.level = a [0], p.i = a [1];
  #if dimension > 1
  p.j = a [2];
  #endif
  #if dimension > 2
  p.k = a [3];
  #endif
  return p;
}

/* a growing array of ints, for each process */
typedef struct { int * a; int n, cap; } Ints;

static void ints_push (Ints * v, const int * x, int n)
{
  if (v->n + n > v->cap)
  {
    v->cap = v->n + n > 2 * v->cap ? v->n + n : 2 * v->cap;
    v->a = realloc (v->a, (size_t) v->cap * sizeof (int));
    if (!v->a)
      abort ();
  }
  memcpy (v->a + v->n, x, (size_t) n * sizeof (int));
  v->n += n;
}

/* what is at index b, as its owner holds it: 0 nothing (a coarser leaf
.. covers it), 1 a leaf, 2 a node with children */
enum { HALO_NONE, HALO_LEAF, HALO_PARENT };

static int answer (Index b)
{
  if (!index_is_allocated (b))
    return HALO_NONE;
  BAMRNode * node = index_node (b);
  /* the owner of b's first leaf holds it and its ancestors whole */
  assert (!(node->flags & BAMR_NODE_IS_PRUNED));
  return node->flags & BAMR_NODE_HAS_CHILDREN ? HALO_PARENT : HALO_LEAF;
}

typedef struct
{
  Index p;
  int owner, answer;
} Query;

static int query_cmp (const void * a, const void * b)
{
  return index_cmp (((const Query *) a)->p, ((const Query *) b)->p);
}

/*
.. The halo: every node this process reads of others, made here, with its
.. ancestors and their families. For each of its nodes, each index beside
.. it is asked of its owner (tree_partition_owner ()), which answers what
.. it holds there (answer ()); then made: a leaf or a node with children
.. beside it, or the coarser leaf covering the index; for a leaf, the
.. upper children of a node with children below it along an axis (the
.. reflux). The children of its own parents it holds already. The nodes
.. made get their owners from the keys.
*/
static void halo_by_queries (void)
{
  int me = tree_rank (), npe = tree_npe ();
  const TreeCaches * tc = tree_caches ();

  /* 1. every index beside this process's nodes, each once */
  Query * qs = NULL;
  long nq = 0, cap = 0;
  for (int q = 0; q < tc->nall; q++)
  {
    if (tc->all [q].node->pid != me)
      continue;
    for (int dir = 0; dir < TREE_NDIRS; dir++)
    {
      Index b;
      if (dir == TREE_DIR_CENTRE || !index_dir (tc->all [q].p, dir, &b))
        continue;
      if (nq == cap)
      {
        cap = cap ? 2 * cap : 1024;
        qs = realloc (qs, (size_t) cap * sizeof (Query));
        if (!qs)
          abort ();
      }
      qs [nq++] = (Query) { b, tree_partition_owner (b), -1 };
    }
  }
  if (nq)
    qsort (qs, (size_t) nq, sizeof (Query), query_cmp);
  long u = 0;
  for (long q = 0; q < nq; q++)
    if (!u || query_cmp (&qs [u - 1], &qs [q]))
      qs [u++] = qs [q];
  nq = u;

  /* 2. asked of their owners, in that order; this process's answered here */
  Ints * ask = calloc ((size_t) npe, sizeof (Ints));
  int ** out = malloc ((size_t) npe * sizeof (int *));
  int ** in = malloc ((size_t) npe * sizeof (int *));
  int * nout = calloc ((size_t) npe, sizeof (int));
  int * nin = calloc ((size_t) npe, sizeof (int));
  if (!ask || !out || !in || !nout || !nin)
    abort ();
  for (long q = 0; q < nq; q++)
    if (qs [q].owner == me)
      qs [q].answer = answer (qs [q].p);
    else
    {
      int a [4];
      put_index (a, qs [q].p);
      ints_push (&ask [qs [q].owner], a, 4);
    }
  for (int r = 0; r < npe; r++)
    out [r] = ask [r].a, nout [r] = ask [r].n;
  exchange_ints (out, nout, in, nin, 201);

  /* 3. what is asked of this process, answered */
  for (int r = 0; r < npe; r++)
  {
    free (ask [r].a);
    ask [r] = (Ints) { NULL, 0, 0 };
    for (int k = 0; k < nin [r]; k += 4)
    {
      int a = answer (get_index (in [r] + k));
      ints_push (&ask [r], &a, 1);
    }
    free (in [r]);
    out [r] = ask [r].a, nout [r] = ask [r].n;
  }
  exchange_ints (out, nout, in, nin, 202);
  for (int r = 0; r < npe; r++)
    nout [r] = 0;                     /* the answers read so far */
  for (long q = 0; q < nq; q++)
    if (qs [q].owner != me)
    {
      int r = qs [q].owner;
      assert (nout [r] < nin [r]);
      qs [q].answer = in [r][nout [r]++];
    }

  /* 4. the nodes read, made: the list of this process's nodes first, as
   * making changes the tree's lists */
  int nown = 0;
  Index * own = malloc ((size_t) (tc->nall ? tc->nall : 1) * sizeof (Index));
  char * leaf = malloc ((size_t) (tc->nall ? tc->nall : 1));
  if (!own || !leaf)
    abort ();
  for (int q = 0; q < tc->nall; q++)
    if (tc->all [q].node->pid == me)
    {
      leaf [nown] = (tc->all [q].node->flags & BAMR_NODE_LEAF) != 0;
      own [nown++] = tc->all [q].p;
    }
  for (int k = 0; k < nown; k++)
    for (int dir = 0; dir < TREE_NDIRS; dir++)
    {
      Query key;
      if (dir == TREE_DIR_CENTRE || !index_dir (own [k], dir, &key.p))
        continue;
      const Query * f = bsearch (&key, qs, (size_t) nq, sizeof (Query),
                                 query_cmp);
      assert (f && f->answer >= 0);
      if (f->owner == me)
        continue;
      if (f->answer == HALO_NONE)
        tree_insert (index_parent (f->p), 1);
      else
        tree_insert (f->p, f->answer == HALO_LEAF);
      /* the reflux: d = -1 along an axis a, a node with children there */
      int a = -1;
      for (int t = 0, w = 1; t < dimension; t++, w *= 3)
        if (dir == TREE_DIR_CENTRE - w)
          a = t;
      if (leaf [k] && a >= 0 && f->answer == HALO_PARENT)
        for (int c = 0; c < TREE_NCHILDREN; c++)
          if ((c >> a) & 1)
          {
            Index x = index_child (f->p, c);
            if (tree_partition_owner (x) != me)
              tree_insert (x, 1);
          }
    }

  /* 5. the owners of the nodes made */
  tc = tree_caches ();
  for (int q = 0; q < tc->nall; q++)
    if (tc->all [q].node->pid < 0)
      tc->all [q].node->pid = tree_partition_owner (tc->all [q].p);

  for (int r = 0; r < npe; r++)
    free (ask [r].a), free (in [r]);
  free (ask), free (out), free (in), free (nout), free (nin);
  free (qs), free (own), free (leaf);
}

/*
.. This process's plan from what it holds: what it receives, from the
.. tree here (plan_needs ()); what it sends, as the others ask, each its
.. receive lists (7 ints a node: parents, kind, level, the index). The
.. levels, depth and depth_min, those of the whole tree.
*/
static TreeMpiPlan * plan_local (int depth, int depth_min)
{
  int me = tree_rank (), npe = tree_npe ();
  TreeMpiPlan * plan = calloc (1, sizeof (TreeMpiPlan));
  Ints * ask = calloc ((size_t) npe, sizeof (Ints));
  int ** out = malloc ((size_t) npe * sizeof (int *));
  int ** in = malloc ((size_t) npe * sizeof (int *));
  int * nout = calloc ((size_t) npe, sizeof (int));
  int * nin = calloc ((size_t) npe, sizeof (int));
  if (!plan || !ask || !out || !in || !nout || !nin)
    abort ();
  plan->rank = me, plan->npe = npe;
  plan->depth = depth, plan->depth_min = depth_min;

  Needs v = { NULL, 0, 0, 1 };
  plan_needs (&v, me, depth_min);
  for (long q = 0; q < v.n; q++)
  {
    const Need * f = &v.a [q];
    int a [7] = { f->parents, f->kind, f->level };
    put_index (a + 3, f->it.p);
    ints_push (&ask [f->peer], a, 7);
  }
  for (int r = 0; r < npe; r++)
    out [r] = ask [r].a, nout [r] = ask [r].n;
  exchange_ints (out, nout, in, nin, 203);
  for (int r = 0; r < npe; r++)
    for (int k = 0; k < nin [r]; k += 7)
    {
      Index p = get_index (in [r] + k + 3);
      BAMRNode * X = index_is_allocated (p) ? index_node (p) : NULL;
      if (!X || X->pid != me)
      {
        fprintf (stderr, "tree-mpi: rank %d asked by %d for a node it does "
                 "not own\n", me, r);
        abort ();
      }
      if (v.n == v.cap)
      {
        v.cap = v.cap ? 2 * v.cap : 1024;
        v.a = realloc (v.a, (size_t) v.cap * sizeof (Need));
        if (!v.a)
          abort ();
      }
      v.a [v.n++] = (Need) { in [r][k], in [r][k + 1], in [r][k + 2], 1, r,
                             { p, X } };
    }
  plan_lists (plan, &v);
  plan->generation = tree_generation;
  for (int r = 0; r < npe; r++)
    free (ask [r].a), free (in [r]);
  free (ask), free (out), free (in), free (nout), free (nin);
  return plan;
}

/*
.. bytes exchanged with every process, as exchange_ints (): out [r]
.. [0 .. nout [r]) to r, in [r][0 .. nin [r]) from r (malloc'd)
*/
static void exchange_bytes (char ** out, const int * nout, char ** in,
                            int * nin, int tag)
{
  int npe = tree_npe (), me = tree_rank (), nr = 0;
  MPI_Alltoall ((void *) nout, 1, MPI_INT, nin, 1, MPI_INT, MPI_COMM_WORLD);
  MPI_Request * req = malloc ((size_t) (2 * npe) * sizeof (MPI_Request));
  if (!req)
    abort ();
  for (int r = 0; r < npe; r++)
  {
    in [r] = NULL;
    if (!nin [r])
      continue;
    in [r] = malloc ((size_t) nin [r]);
    if (!in [r])
      abort ();
    if (r == me)
      memcpy (in [r], out [r], (size_t) nin [r]);
    else
      MPI_Irecv (in [r], nin [r], MPI_BYTE, r, tag, MPI_COMM_WORLD,
                 &req [nr++]);
  }
  for (int r = 0; r < npe; r++)
    if (nout [r] && r != me)
      MPI_Isend (out [r], nout [r], MPI_BYTE, r, tag, MPI_COMM_WORLD,
                 &req [nr++]);
  MPI_Waitall (nr, req, MPI_STATUSES_IGNORE);
  free (req);
}

/* a growing array of bytes */
typedef struct { char * a; size_t n, cap; } Bytes;

static void bytes_push (Bytes * v, const void * x, size_t n)
{
  if (v->n + n > v->cap)
  {
    v->cap = v->n + n > 2 * v->cap ? v->n + n : 2 * v->cap;
    v->a = realloc (v->a, v->cap);
    if (!v->a)
      abort ();
  }
  memcpy (v->a + v->n, x, n);
  v->n += n;
}

/* the halo and the plan made again for this process, from what it holds:
.. its own nodes and ancestors kept, the halo asked, the plan built */
static void distribute_again (void)
{
  tree_mpi_plan_free (applied);
  applied = NULL;
  need_own (tree_rank ());
  prune_unneeded ();
  halo_by_queries ();
  const TreeCaches * tc = tree_caches ();
  int lv [2] = { 0, -(BAMR_MAX_DEPTH + 1) }, all [2];
  for (int q = 0; q < tc->nowned; q++)
  {
    int l = tc->owned [q].p.level;
    lv [0] = l > lv [0] ? l : lv [0];
    lv [1] = -l > lv [1] ? -l : lv [1];
  }
  MPI_Allreduce (lv, all, 2, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
  applied = plan_local (all [0], -all [1]);
  ghosts_mark_all_dirty ();
}

int tree_mpi_repartition (const double * weight)
{
  int me = tree_rank (), npe = tree_npe ();
  const TreeCaches * tc = tree_caches ();
  int n = tc->nowned, bad = 0, any;
  double wsum = 0;
  for (int q = 0; q < n; q++)
  {
    double w = weight ? weight [q] : 1;
    bad |= !(w > 0);
    wsum += w;
  }
  MPI_Allreduce (&bad, &any, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
  if (any || !made)
    return -1;

  /*
  .. 1. the new cut: the leaves before this process's, and their weight,
  .. by prefix sums; then each leaf's process as tree_partition () gives
  .. it (in integers for unit weights)
  */
  long cnt = n, before = 0, total;
  double wbefore = 0, wtotal;
  MPI_Exscan (&cnt, &before, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
  MPI_Exscan (&wsum, &wbefore, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  if (me == 0)
    before = 0, wbefore = 0;          /* Exscan leaves rank 0's undefined */
  MPI_Allreduce (&cnt, &total, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce (&wsum, &wtotal, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  uint64_t * kmin = malloc ((size_t) npe * sizeof (uint64_t));
  uint64_t * kall = malloc ((size_t) npe * sizeof (uint64_t));
  double * load = calloc ((size_t) npe, sizeof (double));
  double * lall = malloc ((size_t) npe * sizeof (double));
  long * count = calloc ((size_t) npe, sizeof (long));
  long * call = malloc ((size_t) npe * sizeof (long));
  if (!kmin || !kall || !load || !lall || !count || !call)
    abort ();
  for (int r = 0; r < npe; r++)
    kmin [r] = UINT64_MAX;
  for (int q = 0; q < n; q++)
  {
    double w = weight ? weight [q] : 1;
    long g = before + q;
    int pid = weight
      ? (int) floor (npe * (wbefore + w / 2) / wtotal)
      : (int) ((long long) npe * (2 * (long long) g + 1) / (2 * (long long) total));
    wbefore += w;
    pid = pid < 0 ? 0 : pid >= npe ? npe - 1 : pid;
    uint64_t key = tree_morton_key (tc->owned [q].p);
    kmin [pid] = key < kmin [pid] ? key : kmin [pid];
    load [pid] += w;
    count [pid]++;
  }
  /* every process's first key, load and leaves, everywhere */
  MPI_Allreduce (kmin, kall, npe, MPI_UNSIGNED_LONG_LONG, MPI_MIN,
                 MPI_COMM_WORLD);
  MPI_Allreduce (load, lall, npe, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce (count, call, npe, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
  int * first = malloc ((size_t) (npe + 1) * sizeof (int));
  uint64_t * key = malloc ((size_t) (npe + 1) * sizeof (uint64_t));
  if (!first || !key)
    abort ();
  first [0] = 0;
  for (int r = 0; r < npe; r++)
    first [r + 1] = first [r] + (int) call [r];
  key [npe] = UINT64_MAX;
  for (int r = npe - 1; r >= 0; r--)
    key [r] = call [r] ? kall [r] : key [r + 1];
  key [0] = 0;
  double most = 0;
  for (int r = 0; r < npe; r++)
    most = lall [r] > most ? lall [r] : most;

  /*
  .. 2. every node of this process that changes hands, sent to its new
  .. owner, found from the new keys: its index, whether a leaf, the blocks
  .. of every scalar. The old partition's owners are the pids still
  */
  free (part.leaf_first), free (part.load), free (part.key_first);
  part = (TreePartition) {
    .npe = npe, .nleaves = (int) total, .leaf_first = first, .load = lall,
    .key_first = key, .imbalance = total ? most / (wtotal / npe) : 1,
    .generation = tree_generation,
  };
  int ns = 0, ids [NSCALAR_MAX];
  for (int s = 0; s < NSCALAR_MAX; s++)
    if (scalar_in_use (s))
      ids [ns++] = s;
  size_t bytes = (size_t) tree_block_shape ().cells * sizeof (Real);
  Bytes * box = calloc ((size_t) npe, sizeof (Bytes));
  char ** out = malloc ((size_t) npe * sizeof (char *));
  char ** in = malloc ((size_t) npe * sizeof (char *));
  int * nout = calloc ((size_t) npe, sizeof (int));
  int * nin = calloc ((size_t) npe, sizeof (int));
  if (!box || !out || !in || !nout || !nin)
    abort ();
  int sent = 0;
  for (int q = 0; q < tc->nall; q++)
  {
    BAMRNode * node = tc->all [q].node;
    int r = tree_partition_owner (tc->all [q].p);
    if (node->pid != me || r == me)
      continue;
    int head [5] = { (node->flags & BAMR_NODE_LEAF) != 0 };
    put_index (head + 1, tc->all [q].p);
    bytes_push (&box [r], head, sizeof (head));
    for (int k = 0; k < ns; k++)
      bytes_push (&box [r], block_data (node->s, ids [k]), bytes);
    sent++;
  }
  for (int r = 0; r < npe; r++)
  {
    assert (box [r].n <= INT_MAX);
    out [r] = box [r].a, nout [r] = (int) box [r].n;
  }
  exchange_bytes (out, nout, in, nin, 204);
  /* the nodes received, made here, their blocks filled; the leaves
   * first, so that a parent is made as their ancestor */
  size_t record = 5 * sizeof (int) + (size_t) ns * bytes;
  for (int pass = 1; pass >= 0; pass--)
    for (int r = 0; r < npe; r++)
      for (int k = 0; k + (int) record <= nin [r]; k += (int) record)
      {
        int head [5];
        memcpy (head, in [r] + k, sizeof (head));
        if (head [0] != pass)
          continue;
        Index p = get_index (head + 1);
        tree_insert (p, head [0]);
        BAMRNode * node = index_node (p);
        const char * b = in [r] + k + sizeof (head);
        for (int j = 0; j < ns; j++, b += bytes)
          memcpy (block_data (node->s, ids [j]), b, bytes);
      }

  /* 3. every node's owner by the new keys; then the halo and the plan */
  tc = tree_caches ();
  for (int q = 0; q < tc->nall; q++)
    tc->all [q].node->pid = tree_partition_owner (tc->all [q].p);
  tree_caches_changed ();
  distribute_again ();
  part.generation = tree_generation;  /* the partition of the tree now */

  for (int r = 0; r < npe; r++)
    free (box [r].a), free (in [r]);
  free (box), free (out), free (in), free (nout), free (nin);
  free (kmin), free (kall), free (load), free (count), free (call);
  return sent;
}

int tree_mpi_partition (const double * weight)
{
  int rank = tree_rank (), npe = tree_npe ();
  tree_mpi_end ();
  /* rank 0 partitions: whether it did, and its tree's size */
  int head [2] = { 0, tree_caches ()->nall };
  if (rank == 0)
    head [0] = tree_partition (npe, weight);
  MPI_Bcast (head, 2, MPI_INT, 0, MPI_COMM_WORLD);
  if (head [0])
    return -1;
  const TreeCaches * tc = tree_caches ();
  if (head [1] != tc->nall || tree_pruned ())
  {
    fprintf (stderr, "tree_mpi_partition (): rank %d: the tree is not rank "
             "0's (%d nodes, rank 0 %d; %d pruned)\n", rank, tc->nall,
             head [1], tree_pruned ());
    abort ();
  }
  /* the others: every node's pid, in Morton order, and the loads */
  int * pid = malloc ((size_t) (tc->nall ? tc->nall : 1) * sizeof (int));
  double * load = malloc ((size_t) npe * sizeof (double));
  if (!pid || !load)
    abort ();
  if (rank == 0)
  {
    for (int q = 0; q < tc->nall; q++)
      pid [q] = tc->all [q].node->pid;
    memcpy (load, part.load, (size_t) npe * sizeof (double));
  }
  MPI_Bcast (pid, tc->nall, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast (load, npe, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  if (rank)
  {
    for (int q = 0; q < tc->nall; q++)
      tc->all [q].node->pid = pid [q];
    partition_record (npe, load);
  }
  else
    free (load);
  free (pid);

  /*
  .. each process keeps its own nodes and their ancestors only, then
  .. learns its halo from the owners and makes its plan from what it holds:
  .. no process uses the whole tree from here
  */
  distribute_again ();
  return 0;
}

/* the blocks of s [0..n) of each node of a list, to or from buf */
static void pack (const TreeMpiList * list, int n, const int * s, char * buf,
                  size_t bytes, int out)
{
  for (int q = 0; q < list->n; q++)
    for (int k = 0; k < n; k++, buf += bytes)
    {
      Real * b = block_data (list->item [q].node->s, s [k]);
      if (out)
        memcpy (buf, b, bytes);
      else
        memcpy (b, buf, bytes);
    }
}

void tree_mpi_exchange (int parents, int kind, int level, int n, const int * s)
{
  if (!applied || n <= 0 || level < 0 || level > BAMR_MAX_DEPTH + 1)
    return;
  const TreeMpiExchange * ex = &applied->ex [parents != 0][kind][level];
  int nm = ex->nsend + ex->nrecv;
  if (!nm)
    return;
  size_t bytes = (size_t) tree_block_shape ().cells * sizeof (Real);
  int tag = ((parents != 0) * TREE_MPI_KINDS + kind) * (BAMR_MAX_DEPTH + 2)
            + level;
  MPI_Request * req = malloc ((size_t) nm * sizeof (MPI_Request));
  char ** buf = malloc ((size_t) nm * sizeof (char *));
  if (!req || !buf)
    abort ();
  for (int k = 0; k < nm; k++)
  {
    const TreeMpiList * list = k < ex->nrecv ? &ex->recv [k]
                                             : &ex->send [k - ex->nrecv];
    size_t len = (size_t) list->n * (size_t) n * bytes;
    assert (len <= INT_MAX);
    buf [k] = malloc (len);
    if (!buf [k])
      abort ();
    if (k < ex->nrecv)
      MPI_Irecv (buf [k], (int) len, MPI_BYTE, list->peer, tag, MPI_COMM_WORLD,
                 &req [k]);
    else
    {
      pack (list, n, s, buf [k], bytes, 1);
      MPI_Isend (buf [k], (int) len, MPI_BYTE, list->peer, tag, MPI_COMM_WORLD,
                 &req [k]);
    }
  }
  MPI_Waitall (nm, req, MPI_STATUSES_IGNORE);
  for (int k = 0; k < ex->nrecv; k++)
    pack (&ex->recv [k], n, s, buf [k], bytes, 0);
  for (int k = 0; k < nm; k++)
    free (buf [k]);
  free (buf), free (req);
}

#endif /* _MPI */
