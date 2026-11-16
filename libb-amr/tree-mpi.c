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
  free (plan->peers);
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

/* ---------------------------------------------------------------------- */
/* the skeleton, the peers, point-to-point exchanges                      */
/* ---------------------------------------------------------------------- */

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

/*
.. bytes with the given processes only, point to point: out [k] to to [k],
.. in [k] from from [k] (malloc'd, NULL if none); the sizes first. Each
.. side must name the other: a size, 0 or not, goes to every one named
*/
static void exchange_with (int nto, const int * to, char ** out,
                           const int * nout, int nfrom, const int * from,
                           char ** in, int * nin, int tag)
{
  MPI_Request * req = malloc ((size_t) (nto + nfrom + 1)
                              * sizeof (MPI_Request));
  if (!req)
    abort ();
  int nr = 0;
  for (int k = 0; k < nfrom; k++)
    MPI_Irecv (&nin [k], (int) sizeof (int), MPI_BYTE, from [k], tag,
               MPI_COMM_WORLD, &req [nr++]);
  for (int k = 0; k < nto; k++)
    MPI_Isend ((void *) &nout [k], (int) sizeof (int), MPI_BYTE, to [k], tag,
               MPI_COMM_WORLD, &req [nr++]);
  MPI_Waitall (nr, req, MPI_STATUSES_IGNORE);
  nr = 0;
  for (int k = 0; k < nfrom; k++)
  {
    in [k] = NULL;
    if (!nin [k])
      continue;
    in [k] = malloc ((size_t) nin [k]);
    if (!in [k])
      abort ();
    MPI_Irecv (in [k], nin [k], MPI_BYTE, from [k], tag + 1, MPI_COMM_WORLD,
               &req [nr++]);
  }
  for (int k = 0; k < nto; k++)
    if (nout [k])
      MPI_Isend (out [k], nout [k], MPI_BYTE, to [k], tag + 1,
                 MPI_COMM_WORLD, &req [nr++]);
  MPI_Waitall (nr, req, MPI_STATUSES_IGNORE);
  free (req);
}

/*
.. N (p), the indices related to a node p by what the plans read, either
.. way, each given to f: p's same-level neighbours (ghosts, same and fine),
.. the coarser leaf covering an absent one (coarse); for a leaf, the
.. children of a neighbour with children that touch it (the reflux, and
.. the other side of a coarse neighbour); p's parent and children (the
.. restriction). Read from the tree here, which must hold them: the
.. skeleton (own nodes, ancestors and N of those)
*/
static void related (Index p, BAMRNode * node, void (* f) (Index, void *),
                     void * arg)
{
  for (int dir = 0; dir < TREE_NDIRS; dir++)
  {
    Index b;
    if (dir == TREE_DIR_CENTRE || !index_dir (p, dir, &b))
      continue;
    if (!index_is_allocated (b))
    {
      f (index_parent (b), arg);
      continue;
    }
    f (b, arg);
    if (!(node->flags & BAMR_NODE_LEAF)
        || !(index_node (b)->flags & BAMR_NODE_HAS_CHILDREN))
      continue;
    /* the children of b on p's side: along each axis with d != 0, the
     * half towards p */
    for (int c = 0; c < TREE_NCHILDREN; c++)
    {
      int touch = 1, t = dir;
      for (int a = 0; a < dimension; a++, t /= 3)
      {
        int d = t % 3 - 1;
        touch &= !d || ((c >> a) & 1) == (d < 0);
      }
      if (touch)
        f (index_child (b, c), arg);
    }
  }
  if (p.level > 0)
    f (index_parent (p), arg);
  if (node->flags & BAMR_NODE_HAS_CHILDREN)
    for (int c = 0; c < TREE_NCHILDREN; c++)
      f (index_child (p, c), arg);
}

static void mark_needed (Index x, void * arg)
{
  (void) arg;
  assert (index_is_allocated (x));
  index_node (x)->flags |= TREE_MPI_NEEDED;
}

/* the ranks of the related nodes, other than this one's, gathered */
typedef struct { int * a; int n, cap, me; } Ranks;

static void add_rank (Index x, void * arg)
{
  Ranks * v = arg;
  int r = index_node (x)->pid;
  if (r == v->me)
    return;
  if (v->n == v->cap)
  {
    v->cap = v->cap ? 2 * v->cap : 64;
    v->a = realloc (v->a, (size_t) v->cap * sizeof (int));
    if (!v->a)
      abort ();
  }
  v->a [v->n++] = r;
}

static int int_cmp (const void * a, const void * b)
{
  int x = *(const int *) a, y = *(const int *) b;
  return (x > y) - (x < y);
}

/* the place of rank r among the sorted peers, -1 if not one */
static int peer_slot (const TreeMpiPlan * plan, int r)
{
  const int * f = bsearch (&r, plan->peers, (size_t) plan->npeers,
                           sizeof (int), int_cmp);
  return f ? (int) (f - plan->peers) : -1;
}

/*
.. every node's blocks as it is used: a node of this process's or one it
.. receives is active, every other inactive (BAMR_NODE_IS_INACTIVE); a
.. family (the root) has blocks iff it has an active member
*/
static void activate (const TreeMpiPlan * plan)
{
  const TreeCaches * tc = tree_caches ();
  int me = tree_rank ();
  for (int q = 0; q < tc->nall; q++)
  {
    BAMRNode * node = tc->all [q].node;
    if (node->pid == me)
      node->flags &= ~BAMR_NODE_IS_INACTIVE;
    else
      node->flags |= BAMR_NODE_IS_INACTIVE;
  }
  for (int parents = 0; parents < 2; parents++)
    for (int kind = 0; kind < TREE_MPI_KINDS; kind++)
      for (int l = 0; l < BAMR_MAX_DEPTH + 2; l++)
      {
        const TreeMpiExchange * ex = &plan->ex [parents][kind][l];
        for (int k = 0; k < ex->nrecv; k++)
          for (int j = 0; j < ex->recv [k].n; j++)
            ex->recv [k].item [j].node->flags &= ~BAMR_NODE_IS_INACTIVE;
      }
  BAMRNode * root = index_node (TREE_ROOT_INDEX);
  if (root->flags & BAMR_NODE_IS_INACTIVE)
    grid_scalar_root_free (root->s);
  else
    grid_scalar_root_blocks (root->s);
  for (int q = 0; q < tc->nall; q++)
  {
    if (!(tc->all [q].node->flags & BAMR_NODE_HAS_CHILDREN))
      continue;
    void *** t [TREE_NCHILDREN];
    int on = 0;
    for (int c = 0; c < TREE_NCHILDREN; c++)
    {
      BAMRNode * child = index_node (index_child (tc->all [q].p, c));
      t [c] = child->s;
      on |= !(child->flags & BAMR_NODE_IS_INACTIVE);
    }
    grid_scalar_family_blocks (t, on);
  }
}

/* one node made active now, its family's blocks given if missing */
static void activate_one (Index p)
{
  BAMRNode * node = index_node (p);
  node->flags &= ~BAMR_NODE_IS_INACTIVE;
  if (!p.level)
  {
    grid_scalar_root_blocks (node->s);
    return;
  }
  void *** t [TREE_NCHILDREN];
  for (int c = 0; c < TREE_NCHILDREN; c++)
    t [c] = index_node (index_child (index_parent (p), c))->s;
  grid_scalar_family_blocks (t, 1);
}

/*
.. The skeleton and the plan, from the tree here, every node's pid right:
.. 1. the skeleton: this process's nodes, their ancestors, and N () of
..    each (related ()); the rest pruned;
.. 2. the peers: the ranks of N () of its nodes, other than its own; the
..    relation goes both ways (N (x) holds p when N (p) holds x, by 2:1
..    balance), so each peer names this process too;
.. 3. the plan: what it receives, from the tree here (plan_needs ()); what
..    it sends, as its peers ask it, each its receive lists (7 ints a
..    node: parents, kind, level, the index), point to point;
.. 4. the blocks: active nodes and families (activate ()).
.. The levels, depth and depth_min, the whole tree's (MPI_Allreduce).
*/
static TreeMpiPlan * skeleton_and_plan (void)
{
  int me = tree_rank (), npe = tree_npe ();

  /* 1. */
  const TreeCaches * tc = tree_caches ();
  int nown = 0;
  TreeItem * own = malloc ((size_t) (tc->nall ? tc->nall : 1)
                           * sizeof (TreeItem));
  if (!own)
    abort ();
  for (int q = 0; q < tc->nall; q++)
    if (tc->all [q].node->pid == me)
      own [nown++] = tc->all [q];
  for (int k = 0; k < nown; k++)
  {
    own [k].node->flags |= TREE_MPI_NEEDED;
    related (own [k].p, own [k].node, mark_needed, NULL);
  }
  prune_unneeded ();

  /* 2. */
  Ranks rs = { NULL, 0, 0, me };
  for (int k = 0; k < nown; k++)
    related (own [k].p, index_node (own [k].p), add_rank, &rs);
  if (rs.n)
    qsort (rs.a, (size_t) rs.n, sizeof (int), int_cmp);
  int u = 0;
  for (int k = 0; k < rs.n; k++)
    if (!u || rs.a [u - 1] != rs.a [k])
      rs.a [u++] = rs.a [k];
  free (own);

  /* 3. */
  TreeMpiPlan * plan = calloc (1, sizeof (TreeMpiPlan));
  if (!plan)
    abort ();
  plan->rank = me, plan->npe = npe;
  plan->npeers = u, plan->peers = rs.a;
  tc = tree_caches ();
  int lv [2] = { 0, -(BAMR_MAX_DEPTH + 1) }, all [2];
  for (int q = 0; q < tc->nowned; q++)
  {
    int l = tc->owned [q].p.level;
    lv [0] = l > lv [0] ? l : lv [0];
    lv [1] = -l > lv [1] ? -l : lv [1];
  }
  MPI_Allreduce (lv, all, 2, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
  plan->depth = all [0], plan->depth_min = -all [1];

  Needs v = { NULL, 0, 0, 1 };
  plan_needs (&v, me, plan->depth_min);
  int m = u ? u : 1;
  Bytes * ask = calloc ((size_t) m, sizeof (Bytes));
  char ** out = malloc ((size_t) m * sizeof (char *));
  char ** in = malloc ((size_t) m * sizeof (char *));
  int * nout = calloc ((size_t) m, sizeof (int));
  int * nin = calloc ((size_t) m, sizeof (int));
  if (!ask || !out || !in || !nout || !nin)
    abort ();
  for (long q = 0; q < v.n; q++)
  {
    const Need * f = &v.a [q];
    int k = peer_slot (plan, f->peer);
    if (k < 0)
    {
      fprintf (stderr, "tree-mpi: rank %d receives from %d, not a peer\n",
               me, f->peer);
      abort ();
    }
    int a [7] = { f->parents, f->kind, f->level };
    put_index (a + 3, f->it.p);
    bytes_push (&ask [k], a, sizeof (a));
  }
  for (int k = 0; k < u; k++)
    out [k] = ask [k].a, nout [k] = (int) ask [k].n;
  exchange_with (u, plan->peers, out, nout, u, plan->peers, in, nin, 203);
  for (int k = 0; k < u; k++)
    for (int b = 0; b + 7 * (int) sizeof (int) <= nin [k];
         b += 7 * (int) sizeof (int))
    {
      int a [7];
      memcpy (a, in [k] + b, sizeof (a));
      Index p = get_index (a + 3);
      BAMRNode * X = index_is_allocated (p) ? index_node (p) : NULL;
      if (!X || X->pid != me)
      {
        fprintf (stderr, "tree-mpi: rank %d asked by %d for a node it does "
                 "not own\n", me, plan->peers [k]);
        abort ();
      }
      if (v.n == v.cap)
      {
        v.cap = v.cap ? 2 * v.cap : 1024;
        v.a = realloc (v.a, (size_t) v.cap * sizeof (Need));
        if (!v.a)
          abort ();
      }
      v.a [v.n++] = (Need) { a [0], a [1], a [2], 1, plan->peers [k],
                             { p, X } };
    }
  plan_lists (plan, &v);
  for (int k = 0; k < u; k++)
    free (ask [k].a), free (in [k]);
  free (ask), free (out), free (in), free (nout), free (nin);

  /* 4. */
  activate (plan);
  plan->generation = tree_generation;
  return plan;
}

/* the plan made again for this process, from the tree here */
static void distribute_again (void)
{
  tree_mpi_plan_free (applied);
  applied = skeleton_and_plan ();
  ghosts_mark_all_dirty ();
}

/* the owner of key in a partition's first keys k [0 .. npe] */
static int owner_of_key (const uint64_t * k, int npe, uint64_t key)
{
  int lo = 0, hi = npe;
  while (hi - lo > 1)
  {
    int mid = (lo + hi) / 2;
    if (k [mid] <= key)
      lo = mid;
    else
      hi = mid;
  }
  return lo;
}

/* a node's structure into a message, no data: a leaf (1) or a node with
.. children here or elsewhere (2); its index */
static void pack_structure (Index x, void * arg)
{
  int head [6] = { index_node (x)->flags & BAMR_NODE_LEAF ? 1 : 2, 0 };
  put_index (head + 2, x);
  bytes_push (arg, head, sizeof (head));
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
      : (int) ((long long) npe * (2 * (long long) g + 1)
               / (2 * (long long) total));
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
  uint64_t * old = part.key_first;    /* the old partition's, kept here */
  free (part.leaf_first), free (part.load);
  part = (TreePartition) {
    .npe = npe, .nleaves = (int) total, .leaf_first = first, .load = lall,
    .key_first = key, .imbalance = total ? most / (wtotal / npe) : 1,
    .generation = tree_generation,
  };

  /*
  .. 2. every node of this process that changes hands, sent to its new
  .. owner (the pids are still the old ones): its index, whether a leaf,
  .. the block of every scalar; with N () of it, structure only, so that
  .. the new owner holds the skeleton around it. The partners from the
  .. two partitions: to the new owners of this process's nodes, from the
  .. old owners of its new piece of the curve
  */
  int ns = 0, ids [NSCALAR_MAX];
  for (int s = 0; s < NSCALAR_MAX; s++)
    if (scalar_in_use (s))
      ids [ns++] = s;
  size_t bytes = (size_t) tree_block_shape ().cells * sizeof (Real);
  int * to = malloc ((size_t) (tc->nall + 1) * sizeof (int)), nto = 0;
  int * from = malloc ((size_t) npe * sizeof (int)), nfrom = 0;
  if (!to || !from)
    abort ();
  for (int q = 0; q < tc->nall; q++)
    if (tc->all [q].node->pid == me)
    {
      int r = owner_of_key (key, npe, tree_morton_key (tc->all [q].p));
      if (r != me)
        to [nto++] = r;
    }
  if (nto)
    qsort (to, (size_t) nto, sizeof (int), int_cmp);
  int u = 0;
  for (int k = 0; k < nto; k++)
    if (!u || to [u - 1] != to [k])
      to [u++] = to [k];
  nto = u;
  if (key [me] < key [me + 1])
    for (int o = owner_of_key (old, npe, key [me]);
         o <= owner_of_key (old, npe, key [me + 1] - 1); o++)
      if (o != me && old [o] < old [o + 1])
        from [nfrom++] = o;
  int mt = nto ? nto : 1, mf = nfrom ? nfrom : 1;
  Bytes * box = calloc ((size_t) mt, sizeof (Bytes));
  char ** out = malloc ((size_t) mt * sizeof (char *));
  int * nout = calloc ((size_t) mt, sizeof (int));
  char ** in = malloc ((size_t) mf * sizeof (char *));
  int * nin = calloc ((size_t) mf, sizeof (int));
  if (!box || !out || !nout || !in || !nin)
    abort ();
  int sent = 0;
  for (int q = 0; q < tc->nall; q++)
  {
    BAMRNode * node = tc->all [q].node;
    Index p = tc->all [q].p;
    if (node->pid != me)
      continue;
    int r = owner_of_key (key, npe, tree_morton_key (p));
    if (r == me)
      continue;
    int k = 0;
    while (to [k] != r)
      k++;
    int head [6] = { node->flags & BAMR_NODE_LEAF ? 1 : 2, 1 };
    put_index (head + 2, p);
    bytes_push (&box [k], head, sizeof (head));
    for (int j = 0; j < ns; j++)
      bytes_push (&box [k], block_data (node->s, ids [j]), bytes);
    related (p, node, pack_structure, &box [k]);
    sent++;
  }
  for (int k = 0; k < nto; k++)
  {
    assert (box [k].n <= INT_MAX);
    out [k] = box [k].a, nout [k] = (int) box [k].n;
  }
  exchange_with (nto, to, out, nout, nfrom, from, in, nin, 204);
  /* the nodes received made here; those with data active, and filled */
  for (int k = 0; k < nfrom; k++)
    for (size_t b = 0; b + sizeof (int [6]) <= (size_t) nin [k]; )
    {
      int head [6];
      memcpy (head, in [k] + b, sizeof (head));
      b += sizeof (head);
      Index p = get_index (head + 2);
      tree_insert (p, head [0] == 1);
      if (!head [1])
        continue;
      activate_one (p);
      BAMRNode * node = index_node (p);
      for (int j = 0; j < ns; j++, b += bytes)
        memcpy (block_data (node->s, ids [j]), in [k] + b, bytes);
    }

  /* 3. every node's owner by the new keys; then the skeleton and plan */
  tc = tree_caches ();
  for (int q = 0; q < tc->nall; q++)
    tc->all [q].node->pid = tree_partition_owner (tc->all [q].p);
  tree_caches_changed ();
  distribute_again ();
  part.generation = tree_generation;  /* the partition of the tree now */

  for (int k = 0; k < nto; k++)
    free (box [k].a);
  for (int k = 0; k < nfrom; k++)
    free (in [k]);
  free (box), free (out), free (nout), free (in), free (nin);
  free (to), free (from), free (old);
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
  /* from here no process uses the whole tree */
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
