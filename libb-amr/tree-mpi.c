#include <assert.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef _MPI
# include <mpi.h>
#endif

#include <cwrap/cache.h>

#include "bamr-error.h"
#include "block-ghosts.h"
#include "tree-mpi.h"
#include "grid-scalar.h"

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
#ifdef _MPI
double tree_mpi_imbalance_max = 1;
#endif

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
    BAMR_ABORT ("out of memory");
  /* the first leaf of each process: pid does not decrease along the curve;
   * its key, or the next process's for one with no leaf */
  for (int r = 0, q = 0; r <= npe; r++)
  {
    while (q < n && tc->leaves [q].node->pid < r)
      q++;
    first [r] = q;
  }
  key [npe] = UINT64_MAX;
  for (int r = npe - 1; r > 0; r--)
    key [r] = first [r] < first [r + 1]
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

/* are p and q the same index? */
static int same_index (Index p, Index q)
{
  int same = p.level == q.level && p.i == q.i;
  #if dimension > 1
  same &= p.j == q.j;
  #endif
  #if dimension > 2
  same &= p.k == q.k;
  #endif
  return same;
}

/*
.. The units of the cut: a family of sibling leaves, never split, or a
.. leaf alone. In a list of leaves in Morton order (the tree's, or a
.. process's own), the leaves of the unit starting at q: TREE_NCHILDREN if
.. leaf q is child 0 of a family whose children are all leaves, the
.. family's leaves following it in the list; else 1
*/
static int unit_at (const TreeItem * leaves, int n, int q)
{
  Index p = leaves [q].p;
  if (!p.level || q + TREE_NCHILDREN > n)
    return 1;
  Index up = index_parent (p);
  for (int c = 0; c < TREE_NCHILDREN; c++)
  {
    Index x = index_child (up, c);
    if (!same_index (leaves [q + c].p, x)
        || !(leaves [q + c].node->flags & BAMR_NODE_LEAF))
      return 1;
  }
  return TREE_NCHILDREN;
}

/* the process of a leaf, the g-th of n (weight w, after weight before, of
.. total), by the middle of the leaf: in integers for unit weights */
static int leaf_pid (int npe, long g, long n, double before, double w,
                     double total, int weighted)
{
  int pid = weighted
    ? (int) floor (npe * (before + w / 2) / total)
    : (int) ((long long) npe * (2 * (long long) g + 1) / (2 * (long long) n));
  return pid < 0 ? 0 : pid >= npe ? npe - 1 : pid;
}

/*
.. Cut the whole tree (one process holds it) into npe pieces of the curve:
.. 1. each leaf to the process whose share of the total weight holds the
..    middle of the leaf (leaf_pid ());
.. 2. each family of leaves whole to the process of its child 0
..    (unit_at ()), so the next coarsening stays on one process;
.. 3. each parent, finest first, to the smallest pid of its children: the
..    process of its first leaf, as tree_partition_owner () finds it;
.. 4. the first keys and the loads recorded (partition_record ()).
.. Nothing changes if npe < 1, a weight is not > 0 or the tree is pruned.
*/
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
    BAMR_ABORT ("out of memory");

  /*
  .. leaf q to the process whose share of the curve holds its middle (in
  .. integers for unit weights, so the pieces are exact); then a family of
  .. leaves (unit_at ()) whole to the process of its child 0
  */
  double before = 0;
  for (int q = 0; q < n; q++)
  {
    double w = weight ? weight [q] : 1;
    tc->leaves [q].node->pid = leaf_pid (npe, q, n, before, w, total,
                                         weight != NULL);
    before += w;
  }
  for (int q = 0; q < n; )
  {
    int m = unit_at (tc->leaves, n, q);
    for (int k = 1; k < m; k++)
      tc->leaves [q + k].node->pid = tc->leaves [q].node->pid;
    q += m;
  }
  for (int q = 0; q < n; q++)
    load [tc->leaves [q].node->pid] += weight ? weight [q] : 1;
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

/*
.. The 3^D - 1 neighbours of p, each sorted by what is there: a node with
.. children (fine), a leaf or pruned node of p's level (same), or no node,
.. so the coarser leaf covering it by 2:1 balance (coarse); with the pid of
.. that node. A direction outside the domain is in none, its pid -1.
*/
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
  cache (Need) l;         /* the needs, l.a [0 .. l.n) */
  int local;              /* receives only: the sends are asked for */
} Needs;

/* by level, then k, j, i: the order of both sides */
static int index_cmp (Index a, Index b)
{
  if (a.level != b.level)
    return a.level < b.level ? -1 : 1;
  #if dimension > 2
  if (a.k != b.k)
    return a.k < b.k ? -1 : 1;
  #endif
  #if dimension > 1
  if (a.j != b.j)
    return a.j < b.j ? -1 : 1;
  #endif
  return (a.i > b.i) - (a.i < b.i);
}

static int need_cmp (const void * pa, const void * pb)
{
  const Need * a = pa, * b = pb;
  /* written out: gcc -O2 keeps a loop over arrays (and a stack check) */
  #define NEED_CMP(f) if (a->f != b->f) return a->f < b->f ? -1 : 1
  NEED_CMP (parents);
  NEED_CMP (kind);
  NEED_CMP (level);
  NEED_CMP (send);
  NEED_CMP (peer);
  #undef NEED_CMP
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
  cache_push (v->l, (Need) { parents, kind, level, send, send ? r : owner,
                             { x, X } });
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
      /* every node's parent (the parents' plan): what the wavelet of
       * adapt_wavelet () reads, with its ghosts, made at the level above
       * just before */
      if (parents && p.level > 0)
        need (v, rank, parents, TREE_MPI_GHOSTS, p.level, r,
              index_parent (p));
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
  if (v->l.n)
    qsort (v->l.a, (size_t) v->l.n, sizeof (Need), need_cmp);
  long u = 0;
  for (long q = 0; q < v->l.n; q++)
    if (!u || need_cmp (&v->l.a [u - 1], &v->l.a [q]))
      v->l.a [u++] = v->l.a [q];
  for (long q = 0; q < u; )
  {
    long e = q;
    const Need * f = &v->l.a [q];
    while (e < u && v->l.a [e].parents == f->parents
           && v->l.a [e].kind == f->kind && v->l.a [e].level == f->level
           && v->l.a [e].send == f->send && v->l.a [e].peer == f->peer)
      e++;
    TreeMpiExchange * ex = &plan->ex [f->parents][f->kind][f->level];
    int * n = f->send ? &ex->nsend : &ex->nrecv;
    TreeMpiList ** list = f->send ? &ex->send : &ex->recv;
    *list = realloc (*list, (size_t) (*n + 1) * sizeof (TreeMpiList));
    TreeItem * item = malloc ((size_t) (e - q) * sizeof (TreeItem));
    if (!*list || !item)
      BAMR_ABORT ("out of memory");
    for (long k = q; k < e; k++)
      item [k - q] = v->l.a [k].it;
    (*list) [(*n)++] = (TreeMpiList) { f->peer, (int) (e - q), item };
    q = e;
  }
  cache_free (v->l);
}

TreeMpiPlan * tree_mpi_plan_new (int rank)
{
  assert (tree_partition_valid () && !tree_pruned ());
  assert (rank >= 0 && rank < part.npe);
  const TreeCaches * tc = tree_caches ();
  TreeMpiPlan * plan = calloc (1, sizeof (TreeMpiPlan));
  if (!plan)
    BAMR_ABORT ("out of memory");
  plan->rank = rank, plan->npe = part.npe;
  plan->depth = tc->depth, plan->depth_min = tree_get_depth_min ();
  plan->generation = tree_generation;
  Needs v = { .local = 0 };
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
  Index * root = malloc ((tc->nall > 0 ? (size_t) tc->nall : 1)
                        * sizeof (Index));
  if (!root)
    BAMR_ABORT ("out of memory");
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

/* n bytes of x appended to the message v (a cache aborts beyond INT_MAX
.. bytes, which MPI's int counts could not send anyway) */
static void bytes_push (cache (char) * v, const void * x, size_t n)
{
  cache_reserve (*v, v->n + (long) n);
  memcpy (v->a + v->n, x, n);
  v->n += (int) n;
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
    BAMR_ABORT ("out of memory");
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
      BAMR_ABORT ("out of memory");
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
typedef struct { cache (int) l; int me; } Ranks;

static void add_rank (Index x, void * arg)
{
  Ranks * v = arg;
  int r = index_node (x)->pid;
  if (r != v->me)
    cache_push (v->l, r);
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
  TreeItem * own = malloc ((tc->nall > 0 ? (size_t) tc->nall : 1)
                           * sizeof (TreeItem));
  if (!own)
    BAMR_ABORT ("out of memory");
  for (int q = 0; q < tc->nall; q++)
    if (tc->all [q].node->pid == me)
      own [nown++] = tc->all [q];
  for (int k = 0; k < nown; k++)
  {
    own [k].node->flags |= TREE_MPI_NEEDED;
    related (own [k].p, own [k].node, mark_needed, NULL);
  }
  /* a leaf out of the skeleton (a sibling a family brings) of a structure
   * no one keeps up to date here: pruned, unknown */
  for (int q = 0; q < tc->nall; q++)
    if (!(tc->all [q].node->flags & TREE_MPI_NEEDED)
        && (tc->all [q].node->flags & BAMR_NODE_LEAF))
      tree_forget (tc->all [q].p);
  prune_unneeded ();

  /* 2. */
  Ranks rs = { .me = me };
  for (int k = 0; k < nown; k++)
    related (own [k].p, index_node (own [k].p), add_rank, &rs);
  if (rs.l.n)
    qsort (rs.l.a, (size_t) rs.l.n, sizeof (int), int_cmp);
  int u = 0;
  for (int k = 0; k < rs.l.n; k++)
    if (!u || rs.l.a [u - 1] != rs.l.a [k])
      rs.l.a [u++] = rs.l.a [k];
  free (own);

  /* 3. */
  TreeMpiPlan * plan = calloc (1, sizeof (TreeMpiPlan));
  if (!plan)
    BAMR_ABORT ("out of memory");
  plan->rank = me, plan->npe = npe;
  plan->npeers = u, plan->peers = rs.l.a;   /* the plan's now */
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

  Needs v = { .local = 1 };
  plan_needs (&v, me, plan->depth_min);
  int m = u ? u : 1;
  cache (char) * ask = calloc ((size_t) m, sizeof * ask);
  char ** out = malloc ((size_t) m * sizeof (char *));
  char ** in = malloc ((size_t) m * sizeof (char *));
  int * nout = calloc ((size_t) m, sizeof (int));
  int * nin = calloc ((size_t) m, sizeof (int));
  if (!ask || !out || !in || !nout || !nin)
    BAMR_ABORT ("out of memory");
  for (int q = 0; q < v.l.n; q++)
  {
    const Need * f = &v.l.a [q];
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
    out [k] = ask [k].a, nout [k] = ask [k].n;
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
      cache_push (v.l, (Need) { a [0], a [1], a [2], 1, plan->peers [k],
                                { p, X } });
    }
  plan_lists (plan, &v);
  for (int k = 0; k < u; k++)
    cache_free (ask [k]), free (in [k]);
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

void tree_mpi_redistribute (void)
{
  distribute_again ();
}

/* ---------------------------------------------------------------------- */
/* adapting the distributed tree                                          */
/* ---------------------------------------------------------------------- */

/* x and y swap their arrays: cwrap refuses to assign a cache */
static void indices_swap (cache (Index) * x, cache (Index) * y)
{
  Index * a = x->a;
  int n = x->n, cap = x->cap;
  x->a = y->a, x->n = y->n, x->cap = y->cap;
  y->a = a, y->n = n, y->cap = cap;
}

static int index_order (const void * a, const void * b)
{
  return index_cmp (* (const Index *) a, * (const Index *) b);
}

/* v sorted, coarsest first, each once */
static void indices_unique (cache (Index) * v)
{
  if (v->n)
    qsort (v->a, (size_t) v->n, sizeof (Index), index_order);
  int u = 0;
  for (int q = 0; q < v->n; q++)
    if (!u || index_cmp (v->a [u - 1], v->a [q]))
      v->a [u++] = v->a [q];
  v->n = u;
}

/*
.. One round's messages with every peer, both ways: to peer k, the
.. indices ask [k] (a request) and every index of made (what this process
.. changed); from each, the same, into asked and done. 4 ints an index,
.. after the two counts
*/
static void peer_round (const TreeMpiPlan * plan, cache (Index) * ask,
                        const cache (Index) * made, cache (Index) * asked,
                        cache (Index) * done, int tag)
{
  int np = plan->npeers, m = np ? np : 1;
  char ** out = malloc ((size_t) m * sizeof (char *));
  char ** in = malloc ((size_t) m * sizeof (char *));
  int * nout = calloc ((size_t) m, sizeof (int));
  int * nin = calloc ((size_t) m, sizeof (int));
  cache (char) * box = calloc ((size_t) m, sizeof * box);
  if (!out || !in || !nout || !nin || !box)
    BAMR_ABORT ("out of memory");
  for (int k = 0; k < np; k++)
  {
    int head [2] = { ask [k].n, made->n }, a [4];
    bytes_push (&box [k], head, sizeof (head));
    for (int q = 0; q < ask [k].n; q++)
      put_index (a, ask [k].a [q]), bytes_push (&box [k], a, sizeof (a));
    for (int q = 0; q < made->n; q++)
      put_index (a, made->a [q]), bytes_push (&box [k], a, sizeof (a));
    out [k] = box [k].a, nout [k] = box [k].n;
  }
  exchange_with (np, plan->peers, out, nout, np, plan->peers, in, nin, tag);
  for (int k = 0; k < np; k++)
  {
    int head [2];
    memcpy (head, in [k], sizeof (head));
    const int * a = (const int *) (in [k] + sizeof (head));
    for (int q = 0; q < head [0]; q++, a += 4)
      cache_push (*asked, get_index (a));
    for (int q = 0; q < head [1]; q++, a += 4)
      cache_push (*done, get_index (a));
    free (in [k]), cache_free (box [k]);
  }
  free (out), free (in), free (nout), free (nin), free (box);
}

/* is anything left to do on any process? */
static int any_left (int left)
{
  int any;
  MPI_Allreduce (&left, &any, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
  return any;
}

/*
.. Divide the leaves of list on the distributed tree, with 2:1 balance
.. kept across processes. Rounds, until no process has work left:
.. 1. here, sweeps until nothing changes: a leaf of this process is
..    divided once every same-level index around it holds a node; else
..    the coarser leaf covering a missing one is divided first, by this
..    process if it owns it, else asked of its owner, and the leaf waits;
.. 2. one round with every peer (peer_round ()): the requests to the
..    owners, and every leaf divided here, which the peers graft
..    (structure only: tree_graft ());
.. 3. the requests received join the waiting leaves; any_left () decides.
.. A cascade crosses a process per round, and the coarsest waiting leaf is
.. divided each round: at most a few rounds a level (the cap aborts).
.. Returns the leaves this process divided. The skeleton and plan are left
.. to the caller (tree_mpi_redistribute ()).
*/
int tree_mpi_refine (int n, const Index * list)
{
  assert (applied);
  const TreeMpiPlan * plan = applied;
  int me = tree_rank (), np = plan->npeers, divided = 0;
  cache (Index) pending;
  for (int q = 0; q < n; q++)
    cache_push (pending, list [q]);
  cache (Index) * ask = calloc ((size_t) (np ? np : 1), sizeof * ask);
  if (!ask)
    BAMR_ABORT ("out of memory");
  for (int round = 0; ; round++)
  {
    /* each round divides at least the coarsest leaf waiting: a few rounds
     * a level at most */
    if (round > 4 * (BAMR_MAX_DEPTH + 1))
    {
      fprintf (stderr, "tree_mpi_refine (): rank %d: no end after %d "
               "rounds\n", me, round);
      abort ();
    }
    /*
    .. the cascade here, sweeps coarsest first until nothing changes: a
    .. leaf is divided once every index of its level beside it holds a
    .. node; a coarser leaf covering one is divided first, here if this
    .. process's, else asked of its owner, and the leaf waits
    */
    cache (Index) made;               /* freed at the end of the round */
    for (int k = 0; k < np; k++)
      cache_clear (ask [k]);
    for (;;)
    {
      indices_unique (&pending);
      cache (Index) next;
      int progress = 0, before = pending.n;
      for (int q = 0; q < pending.n; q++)
      {
        Index p = pending.a [q];
        BAMRNode * node = index_is_allocated (p) ? index_node (p) : NULL;
        if (!node || !(node->flags & BAMR_NODE_LEAF) || node->pid != me
            || p.level >= tree_get_depth_max ())
          continue;                   /* done, or not this process's */
        int blocked = 0;
        for (int dir = 0; dir < TREE_NDIRS; dir++)
        {
          Index b;
          if (dir == TREE_DIR_CENTRE || !index_dir (p, dir, &b)
              || index_is_allocated (b))
            continue;
          Index c = index_parent (b);   /* the coarser leaf, by 2:1 */
          BAMRNode * cover = index_node (c);
          assert (cover && (cover->flags & BAMR_NODE_LEAF));
          blocked = 1;
          if (cover->pid == me)
            cache_push (next, c);
          else
            cache_push (ask [peer_slot (plan, cover->pid)], c);
        }
        if (blocked)
          cache_push (next, p);
        else if (tree_divide (p))
        {
          divided++, progress = 1;
          cache_push (made, p);
        }
      }
      indices_swap (&pending, &next);   /* the old list goes with next */
      indices_unique (&pending);
      if (!progress && pending.n <= before)
        break;
    }

    /* the requests to the owners, what was divided to every peer */
    cache (Index) asked;
    cache (Index) done;
    int sent = made.n;
    for (int k = 0; k < np; k++)
    {
      indices_unique (&ask [k]);
      sent += ask [k].n;
    }
    peer_round (plan, ask, &made, &asked, &done, 205);
    for (int q = 0; q < done.n; q++)
    {
      Index p = done.a [q];
      if (index_is_allocated (p) && (index_node (p)->flags & BAMR_NODE_LEAF)
          && index_node (p)->pid != me)
        tree_graft (p);
    }
    for (int q = 0; q < asked.n; q++)
      cache_push (pending, asked.a [q]);
    if (!any_left (pending.n > 0 || sent > 0))
      break;
  }
  for (int k = 0; k < np; k++)
    cache_free (ask [k]);
  free (ask);
  return divided;
}

/* may family p be merged by 2:1 balance? a leaf of a neighbour's touching
.. children, each: what the skeleton holds (they are in N () of p's
.. children) */
static int may_merge (Index p)
{
  for (int dir = 0; dir < TREE_NDIRS; dir++)
  {
    Index b;
    if (dir == TREE_DIR_CENTRE || !index_dir (p, dir, &b)
        || !index_is_allocated (b))
      continue;
    for (int c = 0; c < TREE_NCHILDREN; c++)
    {
      Index x = index_child (b, c);
      if (tree_child_touches (p, b, c) && index_is_allocated (x)
          && !(index_node (x)->flags & BAMR_NODE_LEAF))
        return 0;
    }
  }
  return 1;
}

/*
.. Merge the families of list (their parents) on the distributed tree,
.. level by level from the finest any process has (one MPI_Allreduce), as
.. one process would: a merge at one level may let a coarser one happen.
.. At each level a family of this process is merged if no neighbour has a
.. touching child with children (may_merge (), the exact 2:1 rule); then
.. one round with the peers, which ungraft what was merged here. No
.. request goes to other processes: a family is never split (the
.. family-aware cut), so its owner decides alone. Returns the families
.. merged here; the skeleton and plan are left to the caller.
*/
int tree_mpi_coarsen (int n, const Index * list)
{
  assert (applied);
  const TreeMpiPlan * plan = applied;
  int me = tree_rank (), np = plan->npeers, merged = 0, top = 0, all;
  for (int q = 0; q < n; q++)
    top = list [q].level > top ? list [q].level : top;
  MPI_Allreduce (&top, &all, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
  cache (Index) * ask = calloc ((size_t) (np ? np : 1), sizeof * ask);
  if (!ask)
    BAMR_ABORT ("out of memory");
  /* level by level, finest first, as one process does: what is merged
   * at a level may let a coarser family merge */
  for (int l = all; l >= 0; l--)
  {
    cache (Index) made;               /* each freed at the scope's end */
    cache (Index) asked;
    cache (Index) done;
    for (int q = 0; q < n; q++)
      if (list [q].level == l && index_is_allocated (list [q])
          && index_node (list [q])->pid == me && may_merge (list [q])
          && tree_merge (list [q]))
      {
        merged++;
        cache_push (made, list [q]);
      }
    peer_round (plan, ask, &made, &asked, &done, 207);
    for (int q = 0; q < done.n; q++)
    {
      Index p = done.a [q];
      if (index_is_allocated (p)
          && (index_node (p)->flags & BAMR_NODE_HAS_CHILDREN))
        tree_ungraft (p);
    }
  }
  free (ask);
  return merged;
}

/* process s's entry of the new partition, from the one process holding
.. its first leaf: that leaf's key, the weight and the leaves before it */
typedef struct { uint64_t key; double before; long nbefore; } PartEntry;

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

/*
.. Cut the distributed tree again and move every node that changes hands,
.. in one pass; no process holds the whole tree. weight [q] is the cost of
.. own leaf q (tree_caches ()->owned order), or NULL for 1 each.
.. 1. the cut: the leaves and the weight before this process's
..    (MPI_Exscan) and the totals; each own leaf by its middle, as
..    tree_partition () would give it; a family of leaves whole to the
..    process of its child 0, told to the peers if it is split; the new
..    first key of each process from the one holding its first leaf,
..    with the weight and leaves before it (MPI_Allgather of each
..    process's last pid, then MPI_Allgatherv of the entries): no
..    reduction, the loads the differences of the weights before;
.. 2. the move: each node of this process whose new owner is another,
..    sent there with the blocks of every scalar in use and N () of it,
..    structure only (related ()); senders and receivers known from the
..    new and old first keys, one round (exchange_with ());
.. 3. every node's pid from the new keys, then the skeleton, peers, plan
..    and blocks again (distribute_again ()).
.. Returns the nodes this process sent, or -1 on every process (nothing
.. changed) if a weight is not > 0 or no partition was made.
*/
int tree_mpi_repartition (const double * weight)
{
  int me = tree_rank (), npe = tree_npe ();
  const TreeCaches * tc = tree_caches ();
  int n = tc->nowned, bad = 0;
  double wsum = 0;
  for (int q = 0; q < n; q++)
  {
    double w = weight ? weight [q] : 1;
    bad |= !(w > 0);
    wsum += w;
  }

  /*
  .. 1. the new cut: the leaves of the whole tree, their weight and the
  .. processes with a weight not > 0 (one MPI_Allreduce), the leaves and
  .. weight before this process's (one MPI_Exscan): doubles, exact for
  .. counts below 2^53; then each leaf's process as tree_partition ()
  .. gives it (in integers for unit weights)
  */
  double sums [3] = { n, wsum, bad }, tot [3], pre [2] = { 0, 0 };
  MPI_Allreduce (sums, tot, 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  if (tot [2] > 0 || !made)
    return -1;
  MPI_Exscan (sums, pre, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  if (me == 0)
    pre [0] = pre [1] = 0;            /* Exscan leaves rank 0's undefined */
  long before = (long) pre [0], total = (long) tot [0];
  double wbefore = pre [1], wtotal = tot [1], wfirst = wbefore;
  /* each leaf by its middle, as tree_partition () */
  int * pid_of = malloc ((size_t) (n ? n : 1) * sizeof (int));
  if (!pid_of)
    BAMR_ABORT ("out of memory");
  for (int q = 0; q < n; q++)
  {
    double w = weight ? weight [q] : 1;
    pid_of [q] = leaf_pid (npe, before + q, total, wbefore, w, wtotal,
                       weight != NULL);
    wbefore += w;
  }
  /*
  .. a family of leaves whole to the process of its child 0: here if this
  .. process holds child 0 (the family's leaves follow it), else from the
  .. owner of child 0, a peer (the siblings are in N () of each other),
  .. which sends it to every peer: parent and process, 5 ints
  */
  const TreeMpiPlan * plan = applied;
  int np = plan ? plan->npeers : 0;
  cache (char) tell;                  /* freed as the function ends */
  for (int q = 0; q < n; )
  {
    Index p = tc->owned [q].p;
    int m = 1, family = p.level > 0, split = 0;
    Index up = family ? index_parent (p) : p;
    for (int c = 0; c < TREE_NCHILDREN && family; c++)
    {
      Index x = index_child (up, c);
      BAMRNode * sib = index_is_allocated (x) ? index_node (x) : NULL;
      family &= sib && (sib->flags & BAMR_NODE_LEAF);
      split |= sib && sib->pid != me;
    }
    if (family && same_index (p, index_child (up, 0)))
    {
      /* the family's leaves here follow child 0 along the curve */
      while (q + m < n && tc->owned [q + m].p.level == p.level
             && same_index (index_parent (tc->owned [q + m].p), up))
        pid_of [q + m] = pid_of [q], m++;
      if (split)
      {
        int a [5];
        put_index (a, up);
        a [4] = pid_of [q];
        bytes_push (&tell, a, sizeof (a));
      }
    }
    q += m;
  }
  {
    int mm = np ? np : 1;
    char ** out = malloc ((size_t) mm * sizeof (char *));
    char ** in = malloc ((size_t) mm * sizeof (char *));
    int * nout = calloc ((size_t) mm, sizeof (int));
    int * nin = calloc ((size_t) mm, sizeof (int));
    if (!out || !in || !nout || !nin)
      BAMR_ABORT ("out of memory");
    for (int k = 0; k < np; k++)
      out [k] = tell.a, nout [k] = tell.n;
    exchange_with (np, plan ? plan->peers : NULL, out, nout,
                   np, plan ? plan->peers : NULL, in, nin, 209);
    for (int k = 0; k < np; k++)
    {
      for (int b = 0; b + (int) sizeof (int [5]) <= nin [k];
           b += (int) sizeof (int [5]))
      {
        int a [5];
        memcpy (a, in [k] + b, sizeof (a));
        Index up = get_index (a);
        /* this process's leaves of that family: children of up */
        for (int q = 0; q < n; q++)
          if (tc->owned [q].p.level == up.level + 1
              && same_index (index_parent (tc->owned [q].p), up))
            pid_of [q] = a [4];
      }
      free (in [k]);
    }
    free (out), free (in), free (nout), free (nin);
  }
  /*
  .. the new table. pid_of does not decrease along the curve, so K_s, the
  .. first key of process s, is the key of the first leaf with pid_of >= s:
  .. one process holds it. This one gives the entries s = last + 1 .. the
  .. pid of its last leaf, last the largest pid of the leaves before its
  .. own: every process's last pid (-1 with no leaf; MPI_Allgather) tells
  .. each how many entries every process gives, in the order of s; the
  .. entries themselves by MPI_Allgatherv. s past every leaf's pid: an
  .. empty process at the end, the key UINT64_MAX
  */
  int lastpid = n ? pid_of [n - 1] : -1, last = -1, before_me = -1;
  int ngot = 0, nmine = 0;
  int * lasts = malloc ((size_t) npe * sizeof (int));
  int * count = malloc ((size_t) npe * sizeof (int));
  int * at = malloc ((size_t) npe * sizeof (int));
  if (!lasts || !count || !at)
    BAMR_ABORT ("out of memory");
  MPI_Allgather (&lastpid, 1, MPI_INT, lasts, 1, MPI_INT, MPI_COMM_WORLD);
  for (int r = 0; r < npe; r++)
  {
    int c = lasts [r] > last ? lasts [r] - last : 0;
    if (r == me)
      before_me = last, nmine = c;
    count [r] = c * (int) sizeof (PartEntry);
    at [r] = ngot * (int) sizeof (PartEntry);
    ngot += c;
    last = lasts [r] > last ? lasts [r] : last;
  }
  PartEntry * mine = malloc ((size_t) (nmine ? nmine : 1)
                             * sizeof (PartEntry));
  PartEntry * all = malloc ((size_t) (ngot ? ngot : 1) * sizeof (PartEntry));
  if (!mine || !all)
    BAMR_ABORT ("out of memory");
  {
    int s = before_me + 1, k = 0;     /* the first entry given here */
    double w0 = wfirst;
    for (int q = 0; q < n && k < nmine; q++)
    {
      for (; k < nmine && pid_of [q] >= s; s++, k++)
        mine [k] = (PartEntry) { tree_morton_key (tc->owned [q].p), w0,
                                 before + q };
      w0 += weight ? weight [q] : 1;
    }
    assert (k == nmine);
  }
  free (pid_of);
  MPI_Allgatherv (mine, count [me], MPI_BYTE, all, count, at, MPI_BYTE,
                  MPI_COMM_WORLD);
  int * first = malloc ((size_t) (npe + 1) * sizeof (int));
  uint64_t * key = malloc ((size_t) (npe + 1) * sizeof (uint64_t));
  double * lall = malloc ((size_t) npe * sizeof (double));
  double * wb = malloc ((size_t) (npe + 1) * sizeof (double));
  if (!first || !key || !lall || !wb)
    BAMR_ABORT ("out of memory");
  for (int r = 0; r <= npe; r++)
  {
    int got = r < ngot;
    key [r] = got ? all [r].key : UINT64_MAX;
    first [r] = got ? (int) all [r].nbefore : (int) total;
    wb [r] = got ? all [r].before : wtotal;
  }
  key [0] = 0, first [0] = 0, wb [0] = 0;   /* the curve starts with 0 */
  double most = 0;
  for (int r = 0; r < npe; r++)
  {
    lall [r] = first [r + 1] > first [r] ? wb [r + 1] - wb [r] : 0;
    most = lall [r] > most ? lall [r] : most;
  }
  free (mine), free (all), free (lasts), free (count), free (at), free (wb);
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
  foreach_scalar (s)                  /* ascending: both sides' order */
    ids [ns++] = s;
  size_t bytes = (size_t) tree_block_shape ().cells * sizeof (Real);
  int * to = malloc ((size_t) (tc->nall + 1) * sizeof (int)), nto = 0;
  /* the old owners of the new piece: a range of processes */
  int lo = 0, hi = -1;
  if (key [me] < key [me + 1])
    lo = owner_of_key (old, npe, key [me]),
    hi = owner_of_key (old, npe, key [me + 1] - 1);
  int * from = malloc ((size_t) (hi - lo + 2) * sizeof (int)), nfrom = 0;
  if (!to || !from)
    BAMR_ABORT ("out of memory");
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
  for (int o = lo; o <= hi; o++)
    if (o != me && old [o] < old [o + 1])
      from [nfrom++] = o;
  int mt = nto ? nto : 1, mf = nfrom ? nfrom : 1;
  cache (char) * box = calloc ((size_t) mt, sizeof * box);
  char ** out = malloc ((size_t) mt * sizeof (char *));
  int * nout = calloc ((size_t) mt, sizeof (int));
  char ** in = malloc ((size_t) mf * sizeof (char *));
  int * nin = calloc ((size_t) mf, sizeof (int));
  if (!box || !out || !nout || !in || !nin)
    BAMR_ABORT ("out of memory");
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
    out [k] = box [k].a, nout [k] = box [k].n;
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
    cache_free (box [k]);
  for (int k = 0; k < nfrom; k++)
    free (in [k]);
  free (box), free (out), free (nout), free (in), free (nin);
  free (to), free (from), free (old);
  return sent;
}

/*
.. The first partition, of a tree every process made alike (tree ()):
.. rank 0 cuts it (tree_partition ()) and broadcasts every node's pid, in
.. Morton order, and the loads; the others check they hold the same tree
.. (its size, nothing pruned) and record the same partition. Then each
.. process keeps its skeleton only (distribute_again ()). Returns -1 on
.. every process if rank 0 could not cut it.
*/
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
  int * pid = malloc ((tc->nall > 0 ? (size_t) tc->nall : 1) * sizeof (int));
  double * load = malloc ((size_t) npe * sizeof (double));
  if (!pid || !load)
    BAMR_ABORT ("out of memory");
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

/* ---------------------------------------------------------------------- */
/* the start: a uniform level, distributed with no process holding it all  */
/* ---------------------------------------------------------------------- */

/* the node of Morton index m at level l (tree_morton_key (p) >> dimension
.. (BAMR_MAX_DEPTH - l)): its bits dealt back to the axes */
static Index index_of_morton (uint64_t m, int l)
{
  Index p;
  memset (&p, 0, sizeof (p));
  p.level = l;
  for (int b = 0; b < l; b++)
  {
    p.i |= (int) ((m >> (dimension * b)) & 1) << b;
    #if dimension > 1
    p.j |= (int) ((m >> (dimension * b + 1)) & 1) << b;
    #endif
    #if dimension > 2
    p.k |= (int) ((m >> (dimension * b + 2)) & 1) << b;
    #endif
  }
  return p;
}

/*
.. The process of family f of a uniform level of n leaves, as
.. tree_partition () cuts it with unit weights: the process whose share of
.. the curve holds the middle of the family's child 0, leaf f 2^D. Exact
.. in integers: npe (2 g + 1) < 2^63 for n <= 2^30, npe < 2^31.
*/
static int start_pid (uint64_t f, uint64_t n, int npe)
{
  uint64_t g = f * TREE_NCHILDREN;
  return (int) ((uint64_t) npe * (2 * g + 1) / (2 * n));
}

/* p and its same-level neighbours held here, leaves at level L: N (p) on
.. a uniform tree, but for the parent and children, which tree_insert ()
.. makes with every family whole */
static void start_around (Index p, int L)
{
  tree_insert (p, p.level == L);
  for (int dir = 0; dir < TREE_NDIRS; dir++)
  {
    Index b;
    if (dir != TREE_DIR_CENTRE && index_dir (p, dir, &b))
      tree_insert (b, b.level == L);
  }
}

/*
.. The distributed tree of a uniform level L, made by each process alone:
.. tree () on several processes, after tree_replicated () made level 1.
.. 1. the cut, a formula of L and npe (start_pid ()): every process's first
..    leaf, key and load, here, no message;
.. 2. this process's leaves and N () of each, its parents (those whose
..    first leaf is its) and N () of each (tree_insert () makes the
..    ancestors, every family whole); every node's pid by the formula;
.. 3. the skeleton, peers, plan and blocks (distribute_again (): one
..    MPI_Allreduce of the levels, one round with the peers); every block
..    0, as the one-process tree's.
.. The same tree, pids, plan and blocks as tree_replicated (), a uniform
.. refinement to L and tree_mpi_partition (), in O(own + halo) memory and
.. work instead of O(2^(D L)) on every process. At most 2^30 leaves.
*/
int tree_mpi_start (int L)
{
  int me = tree_rank (), npe = tree_npe ();
  assert (L >= 1 && L <= tree_get_depth_max ());
  assert (tree_get_depth () == 1 && !tree_pruned ());
  if (dimension * L > 30)
    BAMR_ABORT ("tree_mpi_start (): more than 2^30 leaves at that level");
  tree_mpi_end ();

  /* 1. */
  uint64_t n = (uint64_t) 1 << (dimension * L), nf = n / TREE_NCHILDREN;
  int shift = dimension * (BAMR_MAX_DEPTH - L);   /* a key to a leaf's q */
  int * first = malloc ((size_t) (npe + 1) * sizeof (int));
  uint64_t * key = malloc ((size_t) (npe + 1) * sizeof (uint64_t));
  double * load = malloc ((size_t) npe * sizeof (double));
  if (!first || !key || !load)
    BAMR_ABORT ("out of memory");
  for (int r = 0; r < npe; r++)
  {
    uint64_t lo = 0, hi = nf;         /* the first family of r or after */
    while (lo < hi)
    {
      uint64_t mid = lo + (hi - lo) / 2;
      if (start_pid (mid, n, npe) >= r)
        hi = mid;
      else
        lo = mid + 1;
    }
    first [r] = (int) (lo * TREE_NCHILDREN);
  }
  first [npe] = (int) n;
  key [npe] = UINT64_MAX;
  for (int r = npe - 1; r >= 0; r--)
    key [r] = first [r] < first [r + 1] ? (uint64_t) first [r] << shift
                                        : key [r + 1];
  key [0] = 0;
  double most = 0;
  for (int r = 0; r < npe; r++)
  {
    load [r] = first [r + 1] - first [r];
    most = load [r] > most ? load [r] : most;
  }
  part = (TreePartition) {
    .npe = npe, .nleaves = (int) n, .leaf_first = first, .load = load,
    .key_first = key, .imbalance = most / ((double) n / npe),
  };
  made = 1;

  /* 2. the level 1 leaves tree_replicated () made: parents now, or unknown */
  if (L > 1)
    for (int c = 0; c < TREE_NCHILDREN; c++)
      tree_forget (index_child (TREE_ROOT_INDEX, c));
  for (int q = first [me]; q < first [me + 1]; q++)
  {
    start_around (index_of_morton ((uint64_t) q, L), L);
    /* the parents whose first leaf is q, coarser and coarser; their
     * children, the rest of N (), are there: families are whole, and at
     * level L they are this process's leaves (a family is never split) */
    uint64_t m = (uint64_t) q;
    for (int l = L - 1; l >= 0 && !(m & (TREE_NCHILDREN - 1)); l--)
    {
      m >>= dimension;
      start_around (index_of_morton (m, l), L);
    }
  }
  const TreeCaches * tc = tree_caches ();
  for (int q = 0; q < tc->nall; q++)
    tc->all [q].node->pid =
      start_pid ((tree_morton_key (tc->all [q].p) >> shift) / TREE_NCHILDREN,
                 n, npe);
  tree_caches_changed ();             /* the owners changed: owned */

  /* 3. */
  distribute_again ();
  size_t bytes = (size_t) tree_block_shape ().cells * sizeof (Real);
  tc = tree_caches ();
  for (int q = 0; q < tc->nall; q++)
  {
    void *** t = tc->all [q].node->s;
    foreach_scalar (s)
      if (t [s / SCALAR_ARRAY_WIDTH]
          && t [s / SCALAR_ARRAY_WIDTH][s % SCALAR_ARRAY_WIDTH])
        memset (block_data (t, s), 0, bytes);
  }
  part.generation = tree_generation;
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

/*
.. One exchange of the plan applied: the blocks of the n scalars s, for the
.. nodes of the lists of (parents, kind, level), with every peer that has
.. one. A buffer a list, the receives posted with the sends, one
.. MPI_Waitall, then unpacked. The tag is unique to (parents, kind,
.. level), so exchanges of different lists never match each other.
*/
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
    BAMR_ABORT ("out of memory");
  for (int k = 0; k < nm; k++)
  {
    const TreeMpiList * list = k < ex->nrecv ? &ex->recv [k]
                                             : &ex->send [k - ex->nrecv];
    size_t len = (size_t) list->n * (size_t) n * bytes;
    assert (len <= INT_MAX);
    buf [k] = malloc (len);
    if (!buf [k])
      BAMR_ABORT ("out of memory");
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
