#include <assert.h>
#include <limits.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#ifdef _MPI
# include <mpi.h>
#endif

#include "block-ghosts.h"
#include "tree-mpi.h"

/* a node needed by the process being pruned for, set only during
.. tree_prune (): a bit no other code uses */
#define TREE_MPI_NEEDED (1 << 30)

/* ---------------------------------------------------------------------- */
/* 1. the partition                                                       */
/* ---------------------------------------------------------------------- */

static TreePartition part;
static int made = 0, npruned = 0;

void tree_partition_free (void)
{
  free (part.leaf_first);
  free (part.load);
  part = (TreePartition) { 0 };
  made = 0;
}

const TreePartition * tree_get_partition (void)
{
  return made ? &part : NULL;
}

int tree_partition_valid (void)
{
  return made && part.generation == tree_generation;
}

int tree_partition (int npe, const double * weight)
{
  if (npe < 1)
    return -1;
  if (npruned)
    tree_unprune ();
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
  int * first = malloc ((size_t) (npe + 1) * sizeof (int));
  double * load = calloc ((size_t) npe, sizeof (double));
  if (!first || !load)
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
  /* the first leaf of each process: pid does not decrease along the curve */
  for (int r = 0, q = 0; r <= npe; r++)
  {
    while (q < n && tc->leaves [q].node->pid < r)
      q++;
    first [r] = q;
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

  double mean = (weight ? total : n) / npe, most = 0;
  for (int r = 0; r < npe; r++)
    most = load [r] > most ? load [r] : most;
  tree_partition_free ();
  part = (TreePartition) {
    .npe = npe, .nleaves = n, .leaf_first = first, .load = load,
    .imbalance = n ? most / mean : 1, .generation = tree_generation,
  };
  made = 1;
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

static Index index_parent (Index p)
{
  p.level--;
  p.i >>= 1;
  #if dimension > 1
  p.j >>= 1;
  #endif
  #if dimension > 2
  p.k >>= 1;
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
  TreeNeighbours m = { 0, 0, 0 };
  for (int dir = 0; dir < TREE_NDIRS; dir++)
  {
    Index b;
    if (dir == TREE_DIR_CENTRE || !index_dir (p, dir, &b))
      continue;
    uint32_t bit = (uint32_t) 1 << dir;
    if (index_is_allocated (b))
    {
      if (index_node (b)->flags & BAMR_NODE_HAS_CHILDREN)
        m.fine |= bit;
      else
        m.same |= bit;
    }
    else
    {
      assert (p.level > 0 && index_is_allocated (index_parent (b)));
      m.coarse |= bit;
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
  else if (owner == rank)
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

TreeMpiPlan * tree_mpi_plan_new (int rank)
{
  assert (tree_partition_valid () && !npruned);
  assert (rank >= 0 && rank < part.npe);
  const TreeCaches * tc = tree_caches ();
  TreeMpiPlan * plan = calloc (1, sizeof (TreeMpiPlan));
  if (!plan)
    abort ();
  plan->rank = rank, plan->npe = part.npe;
  plan->depth = tc->depth, plan->depth_min = tree_get_depth_min ();

  Needs v = { NULL, 0, 0 };
  for (int parents = 0; parents < 2; parents++)
  {
    /* the ghosts of every node updated: the leaves, or every node */
    for (int q = 0; q < tc->nall; q++)
    {
      Index p = tc->all [q].p;
      BAMRNode * node = tc->all [q].node;
      int leaf = node->flags & BAMR_NODE_LEAF, r = node->pid;
      if (!parents && !leaf)
        continue;
      TreeNeighbours m = tree_neighbours (p);
      for (int dir = 0; dir < TREE_NDIRS; dir++)
      {
        Index b;
        uint32_t bit = (uint32_t) 1 << dir;
        if (!((m.same | m.fine | m.coarse) & bit))
          continue;
        index_dir (p, dir, &b);
        need (&v, rank, parents, TREE_MPI_GHOSTS, p.level, r,
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
              need (&v, rank, parents, TREE_MPI_FACES, p.level, r,
                    index_child (lo, c));
        }
    }
    /* the restriction: each parent, from its children */
    int lowest = parents ? 0 : plan->depth_min;
    for (int l = lowest; l < tc->depth; l++)
      for (int b = tc->parent_first [l]; b < tc->parent_first [l + 1]; b++)
        for (int c = 0; c < TREE_NCHILDREN; c++)
          need (&v, rank, parents, TREE_MPI_RESTRICT, l,
                tc->parents_by_level [b].node->pid,
                index_child (tc->parents_by_level [b].p, c));
  }

  /* sorted, each once, then cut into lists: an exchange, a peer each */
  if (v.n)
    qsort (v.a, (size_t) v.n, sizeof (Need), need_cmp);
  long u = 0;
  for (long q = 0; q < v.n; q++)
    if (!u || need_cmp (&v.a [u - 1], &v.a [q]))
      v.a [u++] = v.a [q];
  for (long q = 0; q < u; )
  {
    long e = q;
    while (e < u && v.a [e].parents == v.a [q].parents
           && v.a [e].kind == v.a [q].kind && v.a [e].level == v.a [q].level
           && v.a [e].send == v.a [q].send && v.a [e].peer == v.a [q].peer)
      e++;
    const Need * f = &v.a [q];
    TreeMpiExchange * ex = &plan->ex [f->parents][f->kind][f->level];
    int * n = v.a [q].send ? &ex->nsend : &ex->nrecv;
    TreeMpiList ** list = v.a [q].send ? &ex->send : &ex->recv;
    *list = realloc (*list, (size_t) (*n + 1) * sizeof (TreeMpiList));
    TreeItem * item = malloc ((size_t) (e - q) * sizeof (TreeItem));
    if (!*list || !item)
      abort ();
    for (long k = q; k < e; k++)
      item [k - q] = v.a [k].it;
    (*list) [(*n)++] = (TreeMpiList) { v.a [q].peer, (int) (e - q), item };
    q = e;
  }
  free (v.a);
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

/* the lists made again; the partition, if it was valid, still is */
static void flags_changed (void)
{
  int valid = tree_partition_valid ();
  tree_caches_changed ();
  if (valid)
    part.generation = tree_generation;
}

int tree_pruned (void)
{
  return npruned;
}

int tree_prune (const TreeMpiPlan * plan)
{
  if (npruned)
    tree_unprune ();
  const TreeCaches * tc = tree_caches ();
  /* needed: the process's nodes and what it receives */
  for (int q = 0; q < tc->nall; q++)
    if (tc->all [q].node->pid == plan->rank)
      tc->all [q].node->flags |= TREE_MPI_NEEDED;
  for (int parents = 0; parents < 2; parents++)
    for (int kind = 0; kind < TREE_MPI_KINDS; kind++)
      for (int l = 0; l < BAMR_MAX_DEPTH + 2; l++)
      {
        const TreeMpiExchange * ex = &plan->ex [parents][kind][l];
        for (int k = 0; k < ex->nrecv; k++)
          for (int q = 0; q < ex->recv [k].n; q++)
            ex->recv [k].item [q].node->flags |= TREE_MPI_NEEDED;
      }
  /*
  .. in Morton order a parent comes right after its children's subtrees:
  .. any [l] says whether a node of level l since the last parent of level
  .. l - 1 has something needed in its subtree
  */
  int any [BAMR_MAX_DEPTH + 2] = { 0 }, n = 0;
  for (int q = 0; q < tc->nall; q++)
  {
    BAMRNode * node = tc->all [q].node;
    int l = tc->all [q].p.level, below = 0;
    if (node->flags & BAMR_NODE_HAS_CHILDREN)
    {
      below = any [l + 1];
      any [l + 1] = 0;
      if (!below)
        node->flags |= BAMR_NODE_IS_PRUNED, n++;
    }
    any [l] |= below || (node->flags & TREE_MPI_NEEDED);
    node->flags &= ~TREE_MPI_NEEDED;
  }
  npruned = n;
  flags_changed ();
  return n;
}

/* every node below p, pruned or not, unpruned: how many were */
static int unprune_below (Index p, BAMRNode * node)
{
  int n = (node->flags & BAMR_NODE_IS_PRUNED) != 0;
  node->flags &= ~BAMR_NODE_IS_PRUNED;
  if (node->flags & BAMR_NODE_HAS_CHILDREN)
    for (int c = 0; c < TREE_NCHILDREN; c++)
    {
      Index q = index_child (p, c);
      n += unprune_below (q, index_node (q));
    }
  return n;
}

int tree_unprune (void)
{
  if (!npruned)
    return 0;
  int n = unprune_below (TREE_ROOT_INDEX, index_node (TREE_ROOT_INDEX));
  assert (n == npruned);
  npruned = 0;
  flags_changed ();
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
  if (!applied)
    return;
  tree_unprune ();
  tree_mpi_plan_free (applied);
  applied = NULL;
}

int tree_mpi_partition (const double * weight)
{
  int rank, npe;
  MPI_Comm_rank (MPI_COMM_WORLD, &rank);
  MPI_Comm_size (MPI_COMM_WORLD, &npe);
  tree_mpi_end ();
  if (tree_partition (npe, weight))
    return -1;
  applied = tree_mpi_plan_new (rank);
  tree_prune (applied);
  ghosts_mark_all_dirty ();
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
