/* tree-mpi.h: the balanced partition of the tree, on a uniform tree and on
 * one refined about a sphere, against references computed apart: the
 * closed form of the unit-weight cut, the loads summed again from the pids,
 * the optimal contiguous cut (bisection with a greedy packing), and each
 * parent's pid from the leaves inside its box. Then the neighbours' maps
 * and ranks, the plans, and pruning, against the geometry of the leaves. */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "tree-mpi.h"
#include "grid-scalar.h"

scalar u;              /* blocks for every node, to be checked */

static int checks, fails;
#define CHECK(c) do { checks++; if (!(c)) { fails++; \
  printf ("FAIL line %d: %s\n", __LINE__, #c); } } while (0)

static int index_coord (Index p, int a)
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
      for (int a = 0; a < dimension; a++) {
        double x = (index_coord (p, a) + 0.5) * h - 0.5;
        d2 += x * x;
      }
      if (fabs (sqrt (d2) - 0.3) < h)
        near [k++] = p;
    }
    for (int q = 0; q < k; q++)
      tree_refine (near [q]);       /* a leaf still, or 0 */
    free (near);
  }
}

/* the smallest bottleneck of a cut of w [0 .. n) into at most npe
 * contiguous runs: bisection on it, each tried by a greedy packing */
static int packs (const double * w, int n, int npe, double cap)
{
  int runs = 1;
  double s = 0;
  for (int q = 0; q < n; q++) {
    if (w [q] > cap)
      return 0;
    if (s + w [q] > cap)
      runs++, s = 0;
    s += w [q];
  }
  return runs <= npe;
}

static double optimal (const double * w, int n, int npe)
{
  double lo = 0, hi = 0;
  for (int q = 0; q < n; q++)
    hi += w [q];
  for (int it = 0; it < 200; it++) {
    double mid = (lo + hi) / 2;
    if (packs (w, n, npe, mid))
      hi = mid;
    else
      lo = mid;
  }
  return hi;
}

/* is leaf L in the box of node p (p's level coarser)? */
static int inside (Index L, Index p)
{
  int s = L.level - p.level, in = s >= 0;
  for (int a = 0; a < dimension && in; a++)
    in = (index_coord (L, a) >> s) == index_coord (p, a);
  return in;
}

/* the leaves' pairs that touch by a face and lie on two processes, of all
 * such pairs: the surface of the pieces. Each pair counted once, from its
 * finer side (from below at one level) */
static double cut_fraction (int * pairs)
{
  const TreeCaches * tc = tree_caches ();
  int all = 0, cut = 0;
  for (int q = 0; q < tc->nleaves; q++) {
    Index p = tc->leaves [q].p;
    int nb = 1 << p.level, pid = tc->leaves [q].node->pid;
    for (int a = 0; a < dimension; a++)
      for (int side = -1; side <= 1; side += 2) {
        Index b = p;
        if (a == 0) b.i += side;
#if dimension > 1
        if (a == 1) b.j += side;
#endif
#if dimension > 2
        if (a == 2) b.k += side;
#endif
        if (index_coord (b, a) < 0 || index_coord (b, a) >= nb)
          continue;
        BAMRNode * o = NULL;
        if (index_is_allocated (b)) {
          o = index_node (b);
          if (!(o->flags & BAMR_NODE_LEAF) || side < 0)
            continue;               /* finer there, or counted from below */
        }
        else {                      /* the coarser leaf (2:1) */
          Index C = b;
          C.level--;
          C.i >>= 1;
#if dimension > 1
          C.j >>= 1;
#endif
#if dimension > 2
          C.k >>= 1;
#endif
          o = index_node (C);
        }
        all++;
        cut += o->pid != pid;
      }
  }
  *pairs = all;
  return all ? (double) cut / all : 0;
}

/* every check of the partition just made, for npe and weight */
static void check (int npe, const double * weight)
{
  const TreeCaches * tc = tree_caches ();
  const TreePartition * tp = tree_get_partition ();
  int n = tc->nleaves;
  CHECK (tp && tp->npe == npe && tp->nleaves == n && tree_partition_valid ());

  /* the pids: in range, never decreasing along the curve; with unit
   * weights the closed form */
  int range = 1, order = 1, closed = 1;
  for (int q = 0; q < n; q++) {
    int pid = tc->leaves [q].node->pid;
    range &= pid >= 0 && pid < npe;
    order &= !q || pid >= tc->leaves [q - 1].node->pid;
    closed &= weight || pid == (int) ((long long) npe * (2 * q + 1) / (2 * n));
  }
  CHECK (range && order && closed);

  /* leaf_first and the loads, again from the pids */
  double * load = calloc ((size_t) npe, sizeof (double)), total = 0, wmax = 0;
  int * count = calloc ((size_t) npe, sizeof (int));
  for (int q = 0; q < n; q++) {
    double w = weight ? weight [q] : 1;
    load [tc->leaves [q].node->pid] += w, count [tc->leaves [q].node->pid]++;
    total += w, wmax = w > wmax ? w : wmax;
  }
  int firsts = tp->leaf_first [0] == 0 && tp->leaf_first [npe] == n, loads = 1;
  int empty = 0, cmin = n, cmax = 0;
  double most = 0;
  for (int r = 0; r < npe; r++) {
    firsts &= tp->leaf_first [r + 1] - tp->leaf_first [r] == count [r];
    loads &= fabs (tp->load [r] - load [r]) <= 1e-12 * total;
    empty += !count [r];
    cmin = count [r] < cmin ? count [r] : cmin;
    cmax = count [r] > cmax ? count [r] : cmax;
    most = load [r] > most ? load [r] : most;
  }
  CHECK (firsts && loads);
  CHECK (fabs (tp->imbalance - most / (total / npe)) <= 1e-12 * tp->imbalance);

  /* balance: unit weights within one leaf, none empty when there are
   * enough; weights within the largest of the optimal cut */
  if (!weight)
    CHECK (cmax - cmin <= 1 && (n < npe || !empty));
  else {
    double opt = optimal (weight, n, npe);
    CHECK (most <= total / npe + wmax * (1 + 1e-12));
    CHECK (most <= opt + wmax * (1 + 1e-12) && most >= opt * (1 - 1e-9));
  }

  /* a parent: the smallest pid of the leaves in its box, by a search */
  int parents = 1;
  for (int b = 0; b < tc->parent_first [tc->depth]; b++) {
    Index p = tc->parents_by_level [b].p;
    int least = npe;
    for (int q = 0; q < n; q++)
      if (inside (tc->leaves [q].p, p) && tc->leaves [q].node->pid < least)
        least = tc->leaves [q].node->pid;
    parents &= tc->parents_by_level [b].node->pid == least;
  }
  CHECK (parents);

  /* the keys: increasing along the leaves, a node's its first leaf's; the
   * owner of every node from them alone, its pid; the leaves of this
   * process (rank 0): those of pid 0, in order */
  int keys = 1, owners = 1, owned = 1, k = 0;
  for (int q = 1; q < n; q++)
    keys &= tree_morton_key (tc->leaves [q].p)
            > tree_morton_key (tc->leaves [q - 1].p);
  for (int q = 0; q < tc->nall; q++) {
    Index p = tc->all [q].p, f = p;
    while (!(index_node (f)->flags & BAMR_NODE_LEAF))
      f = index_child (f, 0);
    keys &= tree_morton_key (p) == tree_morton_key (f);
    owners &= tree_partition_owner (p) == tc->all [q].node->pid;
  }
  for (int q = 0; q < n; q++)
    if (tc->leaves [q].node->pid == 0)
      owned &= k < tc->nowned && tc->owned [k++].node == tc->leaves [q].node;
  CHECK (keys && owners && owned && k == tc->nowned);
  free (load), free (count);
}


/* ---------------------------------------------------------------------- */
/* the bit maps, the plans and pruning, against geometry                  */
/* ---------------------------------------------------------------------- */

/* the nodes of the tree, by Index: a search apart from the level index */
static TreeItem * sorted;
static int nsorted;
static int geometry_ok = 1;        /* the searches below found what they want */

static int item_cmp (const void * pa, const void * pb)
{
  Index a = ((const TreeItem *) pa)->p, b = ((const TreeItem *) pb)->p;
  if (a.level != b.level)
    return a.level < b.level ? -1 : 1;
  for (int ax = dimension - 1; ax >= 0; ax--)
    if (index_coord (a, ax) != index_coord (b, ax))
      return index_coord (a, ax) < index_coord (b, ax) ? -1 : 1;
  return 0;
}

static void sort_nodes (void)
{
  const TreeCaches * tc = tree_caches ();
  free (sorted);
  nsorted = tc->nall;
  sorted = malloc ((size_t) nsorted * sizeof (TreeItem));
  for (int q = 0; q < nsorted; q++)
    sorted [q] = tc->all [q];
  qsort (sorted, (size_t) nsorted, sizeof (TreeItem), item_cmp);
}

static BAMRNode * find (Index p)
{
  TreeItem key = { p, NULL };
  TreeItem * f = bsearch (&key, sorted, (size_t) nsorted, sizeof (TreeItem),
                          item_cmp);
  return f ? f->node : NULL;
}

/* d of direction dir along axis a */
static int dir_d (int dir, int a)
{
  for (int t = 0; t < a; t++)
    dir /= 3;
  return dir % 3 - 1;
}

/*
.. what is at p + d by the leaves that overlap it, on the finest level's
.. integer lattice: 0 outside, 1 same (a leaf of p's level), 2 fine (finer
.. leaves), 3 coarse (a leaf of level l - 1 covering it); *X the node read
*/
static int classify (Index p, int dir, Index * X)
{
  const TreeCaches * tc = tree_caches ();
  int F = tree_get_depth (), l = p.level;
  Index b = p;
  long lo [3], hi [3];
  for (int a = 0; a < dimension; a++) {
    int x = index_coord (p, a) + dir_d (dir, a);
    if (x < 0 || x >= 1 << l)
      return 0;
    if (a == 0) b.i = x;
#if dimension > 1
    if (a == 1) b.j = x;
#endif
#if dimension > 2
    if (a == 2) b.k = x;
#endif
    lo [a] = (long) x << (F - l), hi [a] = (long) (x + 1) << (F - l);
  }
  int coarsest = BAMR_MAX_DEPTH + 1;
  Index cover = b;
  for (int q = 0; q < tc->nleaves; q++) {
    Index L = tc->leaves [q].p;
    int over = 1;
    for (int a = 0; a < dimension && over; a++) {
      long x0 = (long) index_coord (L, a) << (F - L.level);
      long x1 = (long) (index_coord (L, a) + 1) << (F - L.level);
      over = x0 < hi [a] && x1 > lo [a];
    }
    if (over && L.level < coarsest)
      coarsest = L.level, cover = L;
  }
  if (coarsest < l) {
    geometry_ok &= coarsest == l - 1;
    *X = cover;
    return 3;
  }
  *X = b;
  return coarsest == l ? 1 : 2;
}

/* the maps of every node against classify () */
static void check_neighbours (void)
{
  sort_nodes ();
  const TreeCaches * tc = tree_caches ();
  int ok = 1;
  for (int q = 0; q < tc->nall; q++) {
    Index p = tc->all [q].p;
    TreeNeighbours m = tree_neighbours (p);
    for (int dir = 0; dir < TREE_NDIRS; dir++) {
      Index X;
      uint32_t bit = (uint32_t) 1 << dir;
      int c = dir == TREE_DIR_CENTRE ? 0 : classify (p, dir, &X);
      int got = (m.same & bit ? 1 : 0) + (m.fine & bit ? 2 : 0)
                + (m.coarse & bit ? 3 : 0);
      ok &= got == c && !((m.same & m.fine) | (m.same & m.coarse)
                         | (m.fine & m.coarse));
      /* the rank of what is read there: the owner of X, by search */
      BAMRNode * o = c ? find (X) : NULL;
      ok &= c ? o && m.pid [dir] == o->pid : m.pid [dir] == -1;
    }
  }
  CHECK (ok && geometry_ok);
}

/* an entry of a plan: what a process receives, or sends */
typedef struct { int parents, kind, level, peer; Index p; } Entry;

static int entry_cmp (const void * pa, const void * pb)
{
  const Entry * a = pa, * b = pb;
  int ka [4] = { a->parents, a->kind, a->level, a->peer };
  int kb [4] = { b->parents, b->kind, b->level, b->peer };
  for (int q = 0; q < 4; q++)
    if (ka [q] != kb [q])
      return ka [q] < kb [q] ? -1 : 1;
  TreeItem x = { a->p, NULL }, y = { b->p, NULL };
  return item_cmp (&x, &y);
}

typedef struct { Entry * a; int n, cap; } Entries;

static void add (Entries * e, int parents, int kind, int level, Index x, int me)
{
  BAMRNode * X = find (x);
  geometry_ok &= X != NULL;
  if (!X || X->pid == me)
    return;
  if (e->n == e->cap)
    e->cap = e->cap ? 2 * e->cap : 256,
    e->a = realloc (e->a, (size_t) e->cap * sizeof (Entry));
  e->a [e->n++] = (Entry) { parents, kind, level, X->pid, x };
}

static void uniq (Entries * e)
{
  if (e->n)
    qsort (e->a, (size_t) e->n, sizeof (Entry), entry_cmp);
  int u = 0;
  for (int q = 0; q < e->n; q++)
    if (!u || entry_cmp (&e->a [u - 1], &e->a [q]))
      e->a [u++] = e->a [q];
  e->n = u;
}

/* what process me must receive, from the geometry of the leaves */
static Entries expected (int me)
{
  const TreeCaches * tc = tree_caches ();
  Entries e = { NULL, 0, 0 };
  for (int parents = 0; parents < 2; parents++) {
    for (int q = 0; q < tc->nall; q++) {
      Index p = tc->all [q].p, X;
      BAMRNode * node = tc->all [q].node;
      int leaf = node->flags & BAMR_NODE_LEAF;
      if ((!parents && !leaf) || node->pid != me)
        continue;
      for (int dir = 0; dir < TREE_NDIRS; dir++)
        if (dir != TREE_DIR_CENTRE && classify (p, dir, &X))
          add (&e, parents, TREE_MPI_GHOSTS, p.level, X, me);
      for (int a = 0; a < dimension && leaf; a++) {
        int w = 1;
        for (int t = 0; t < a; t++)
          w *= 3;
        if (classify (p, TREE_DIR_CENTRE - w, &X) == 2)
          for (int c = 0; c < TREE_NCHILDREN; c++)
            if ((c >> a) & 1)
              add (&e, parents, TREE_MPI_FACES, p.level, index_child (X, c),
                   me);
      }
    }
    int lowest = parents ? 0 : tree_get_depth_min ();
    for (int q = 0; q < tc->nall; q++) {
      Index p = tc->all [q].p;
      BAMRNode * node = tc->all [q].node;
      if ((node->flags & BAMR_NODE_LEAF) || node->pid != me || p.level < lowest)
        continue;
      for (int c = 0; c < TREE_NCHILDREN; c++)
        add (&e, parents, TREE_MPI_RESTRICT, p.level, index_child (p, c), me);
    }
  }
  uniq (&e);
  return e;
}

/* the lists of a plan, send or receive, as entries */
static Entries entries (const TreeMpiPlan * plan, int send)
{
  Entries e = { NULL, 0, 0 };
  for (int parents = 0; parents < 2; parents++)
    for (int kind = 0; kind < TREE_MPI_KINDS; kind++)
      for (int l = 0; l < BAMR_MAX_DEPTH + 2; l++) {
        const TreeMpiExchange * ex = &plan->ex [parents][kind][l];
        int n = send ? ex->nsend : ex->nrecv;
        const TreeMpiList * list = send ? ex->send : ex->recv;
        for (int k = 0; k < n; k++)
          for (int q = 0; q < list [k].n; q++) {
            if (e.n == e.cap)
              e.cap = e.cap ? 2 * e.cap : 256,
              e.a = realloc (e.a, (size_t) e.cap * sizeof (Entry));
            e.a [e.n++] = (Entry) { parents, kind, l, list [k].peer,
                                    list [k].item [q].p };
          }
      }
  return e;
}

/* are a's lists, each in its order, b's: the same entries, in the same order */
static int same_lists (const TreeMpiList * a, const TreeMpiList * b)
{
  int ok = a->n == b->n;
  for (int q = 0; q < a->n && ok; q++) {
    TreeItem x = { a->item [q].p, NULL }, y = { b->item [q].p, NULL };
    ok = !item_cmp (&x, &y) && a->item [q].node == b->item [q].node;
  }
  return ok;
}

/* is node X (at x) below or at p? */
static int below (Index x, Index p)
{
  int s = x.level - p.level, in = s >= 0;
  for (int a = 0; a < dimension && in; a++)
    in = (index_coord (x, a) >> s) == index_coord (p, a);
  return in;
}

/* the plans of every process, against each other and the geometry */
static void check_plans (int npe)
{
  CHECK (tree_partition (npe, NULL) == 0);
  sort_nodes ();
  TreeMpiPlan ** plan = malloc ((size_t) npe * sizeof (TreeMpiPlan *));
  for (int r = 0; r < npe; r++)
    plan [r] = tree_mpi_plan_new (r);

  /* what each receives: what the geometry says, from the owners */
  int match = 1, peers = 1;
  long nrecv = 0;
  for (int r = 0; r < npe; r++) {
    Entries want = expected (r), got = entries (plan [r], 0);
    nrecv += got.n;
    uniq (&got);                    /* sorted the same way: each once */
    match &= want.n == got.n;
    for (int q = 0; q < want.n && match; q++)
      match = !entry_cmp (&want.a [q], &got.a [q]);
    for (int q = 0; q < got.n; q++)
      peers &= got.a [q].peer == find (got.a [q].p)->pid && got.a [q].peer != r;
    free (want.a), free (got.a);
  }
  CHECK (match && peers && geometry_ok);

  /* what r sends to s is what s receives from r, in the same order */
  int sym = 1, sends = 0, recvs = 0;
  for (int r = 0; r < npe; r++)
    for (int parents = 0; parents < 2; parents++)
      for (int kind = 0; kind < TREE_MPI_KINDS; kind++)
        for (int l = 0; l < BAMR_MAX_DEPTH + 2; l++) {
          const TreeMpiExchange * ex = &plan [r]->ex [parents][kind][l];
          for (int k = 0; k < ex->nsend; k++, sends++) {
            int s = ex->send [k].peer, found = 0;
            const TreeMpiExchange * ey = &plan [s]->ex [parents][kind][l];
            for (int j = 0; j < ey->nrecv; j++)
              if (ey->recv [j].peer == r)
                found += same_lists (&ex->send [k], &ey->recv [j]);
            sym &= found == 1 && s != r;
            for (int q = 0; q < ex->send [k].n; q++)
              sym &= ex->send [k].item [q].node->pid == r;
          }
          for (int k = 1; k < ex->nsend; k++)
            sym &= ex->send [k].peer > ex->send [k - 1].peer;
          for (int k = 1; k < ex->nrecv; k++)
            sym &= ex->recv [k].peer > ex->recv [k - 1].peer;
          recvs += ex->nrecv;
        }
  CHECK (sym && sends == recvs);

  printf ("  %d processes: %ld nodes received\n", npe, nrecv);
  for (int r = 0; r < npe; r++)
    tree_mpi_plan_free (plan [r]);
  free (plan);
}

/* the tree of the tests: level 1, then refined about the sphere */
static void build (int rounds)
{
  tree (5, 4);
  if (rounds)
    refine_sphere (rounds);
}

/*
.. pruning for each process r, on a tree made again for each, against the
.. geometry: needed, r's nodes and those its plan receives; a node kept if
.. it is the root or its parent has something needed strictly below it
.. (found by a search of the boxes); a kept node with children and nothing
.. needed below pruned (a parent received keeps its block, not its
.. children). Every node kept, at its place, with its flags and pid; the
.. rest freed; the counts, the lists, the scalars' blocks
*/
static void check_pruning (int npe, int rounds)
{
  long kept_all = 0, all = 0;
  for (int r = 0; r < npe; r++) {
    tree_destroy ();
    build (rounds);
    CHECK (tree_partition (npe, NULL) == 0);
    sort_nodes ();
    TreeMpiPlan * plan = tree_mpi_plan_new (r);
    Entries got = entries (plan, 0);
    int n = nsorted;
    TreeItem * node = malloc ((size_t) n * sizeof (TreeItem));
    int * flags = malloc ((size_t) n * sizeof (int));
    int * pid = malloc ((size_t) n * sizeof (int));
    char * needed = calloc ((size_t) n, 1), * some = calloc ((size_t) n, 1);
    for (int q = 0; q < n; q++) {
      node [q] = sorted [q];
      flags [q] = sorted [q].node->flags, pid [q] = sorted [q].node->pid;
      needed [q] = pid [q] == r;
    }
    for (int k = 0; k < got.n; k++) {
      TreeItem key = { got.a [k].p, NULL }, * f = bsearch (&key, sorted,
        (size_t) n, sizeof (TreeItem), item_cmp);
      geometry_ok &= f != NULL;
      if (f)
        needed [f - sorted] = 1;
    }
    for (int q = 0; q < n; q++)
      for (int k = 0; k < n && !some [q]; k++)
        some [q] = needed [k] && node [k].p.level > node [q].p.level
                   && below (node [k].p, node [q].p);
    int nkept = 0, nleaves = 0, nstubs = 0, lmin = BAMR_MAX_DEPTH, lmax = 0;
    char * kept = calloc ((size_t) n, 1), * stub = calloc ((size_t) n, 1);
    for (int q = 0; q < n; q++) {
      Index p = node [q].p;
      int up = -1;
      if (p.level) {
        TreeItem key = { index_parent (p), NULL }, * f = bsearch (&key,
          sorted, (size_t) n, sizeof (TreeItem), item_cmp);
        geometry_ok &= f != NULL;
        up = f ? (int) (f - sorted) : -1;
      }
      kept [q] = !p.level || (up >= 0 && some [up]);
      stub [q] = kept [q] && (flags [q] & BAMR_NODE_HAS_CHILDREN) && !some [q];
      nkept += kept [q], nstubs += stub [q];
      if (kept [q] && (flags [q] & BAMR_NODE_LEAF)) {
        nleaves++;
        lmin = p.level < lmin ? p.level : lmin;
        lmax = p.level > lmax ? p.level : lmax;
      }
    }

    int nroots = tree_mpi_prune (plan);
    int place = 1, freed = 1, marks = 1;
    for (int q = 0; q < n; q++) {
      Index p = node [q].p;
      int there = index_is_allocated (p);
      if (!kept [q]) {
        freed &= !there;
        continue;
      }
      BAMRNode * X = there ? index_node (p) : NULL;
      place &= X == node [q].node;
      if (X != node [q].node)
        continue;
      marks &= X->pid == pid [q];
      if (stub [q])
        marks &= X->flags == ((flags [q] & ~(BAMR_NODE_HAS_CHILDREN
                                             | BAMR_NODE_LEAF))
                              | BAMR_NODE_IS_PRUNED);
      else
        marks &= X->flags == flags [q];
    }
    CHECK (place && freed && marks && geometry_ok);
    CHECK (nroots == nstubs && tree_pruned () == nstubs);
    CHECK (tree_get_nblocks () == nkept && tree_get_nleaves () == nleaves);
    /* the levels with leaves: none left, as for a process owning nothing,
     * gives 0 */
    CHECK (nleaves ? tree_get_depth_min () == lmin && tree_get_depth () == lmax
                   : tree_get_depth_min () == 0 && tree_get_depth () == 0);
    const TreeCaches * tc = tree_caches ();
    CHECK (tc->nall == nkept && tc->nleaves == nleaves);
    /* the pruned nodes: in the lists of every node, in none other */
    int lists = 1, listed = 0;
    for (int q = 0; q < tc->nall; q++)
      listed += (tc->all [q].node->flags & BAMR_NODE_IS_PRUNED) != 0;
    lists &= listed == nstubs;
    for (int q = 0; q < tc->nleaves; q++)
      lists &= !(tc->leaves [q].node->flags & BAMR_NODE_IS_PRUNED);
    for (int q = 0; q < tc->parent_first [tc->depth]; q++)
      lists &= !(tc->parents_by_level [q].node->flags & BAMR_NODE_IS_PRUNED);
    CHECK (lists);
    /* the partition stale and refused on a pruned tree; a pruned node,
     * having no children, is not pruned again */
    CHECK (nstubs ? !tree_partition_valid () && tree_partition (npe, NULL) == -1
                  : tree_partition_valid ());
    int again = 0;
    for (int q = 0; q < n; q++)
      if (stub [q])
        again += tree_prune (node [q].p);
    CHECK (!again && tree_pruned () == nstubs);
    /* a scalar made and freed on the pruned tree; the blocks of every node */
    CHECK (grid_scalar_check ());
    {
      scalar t;
      (void) t;
    }
    CHECK (grid_scalar_check ());
    kept_all += nkept, all += n;
    free (got.a), free (node), free (flags), free (pid), free (needed);
    free (some), free (kept), free (stub);
    tree_mpi_plan_free (plan);
  }
  printf ("  %d processes, %d rounds: %.2f of the nodes kept by a pruned "
          "process\n", npe, rounds, (double) kept_all / all);
  tree_destroy ();
  build (rounds);
}

int main (void)
{
  int sizes [] = { 1, 2, 3, 5, 7, 16, 64 };
  int nsizes = sizeof (sizes) / sizeof (sizes [0]);

  /* one process, of rank 0: every node its own */
  CHECK (tree_npe () == 1 && tree_rank () == 0);
  tree (5, 4);
  {
    const TreeCaches * tc = tree_caches ();
    int own = 1;
    for (int q = 0; q < tc->nall; q++)
      own &= tc->all [q].node->pid == tree_rank ();
    CHECK (own && !tree_pruned ());
  }

  /* bad arguments: refused, the partition left as it was */
  CHECK (!tree_get_partition () && !tree_partition_valid ());
  CHECK (tree_partition (0, NULL) == -1 && !tree_get_partition ());
  CHECK (tree_partition (3, NULL) == 0);
  {
    int n = tree_caches ()->nleaves;
    double * w = malloc ((size_t) n * sizeof (double));
    for (int q = 0; q < n; q++)
      w [q] = 1;
    w [n - 1] = 0;
    CHECK (tree_partition (2, w) == -1);
    w [n - 1] = NAN;
    CHECK (tree_partition (2, w) == -1);
    CHECK (tree_get_partition ()->npe == 3 && tree_partition_valid ());
    free (w);
  }

  /* the uniform tree at level 1 (2^dimension leaves), then refined about
   * the sphere: unit weights, then weights of 1 to 7 */
  for (int rounds = 0; rounds <= 3; rounds += 3) {
    if (rounds)
      refine_sphere (rounds);
    const TreeCaches * tc = tree_caches ();
    int n = tc->nleaves;
    double * w = malloc ((size_t) n * sizeof (double));
    for (int q = 0; q < n; q++)
      w [q] = 1 + (q * 2654435761u >> 7) % 7;
    for (int s = 0; s < nsizes; s++) {
      CHECK (tree_partition (sizes [s], NULL) == 0);
      check (sizes [s], NULL);
      CHECK (tree_partition (sizes [s], w) == 0);
      check (sizes [s], w);
    }
    /* more processes than leaves */
    CHECK (tree_partition (n + 3, NULL) == 0);
    check (n + 3, NULL);
    free (w);
  }

  /* the bit maps, the plans, pruning */
  check_neighbours ();
  for (int npe = 1; npe <= 7; npe += 3) {
    check_plans (npe);
    check_pruning (npe, 3);
  }
  /* the uniform tree on more processes than leaves: some own nothing and
   * keep the root only */
  check_pruning (TREE_NCHILDREN + 3, 0);
  tree_destroy ();
  build (3);

  /* the surface of the pieces, measured (not checked) */
  int pairs = 0;
  for (int npe = 4; npe <= 16; npe *= 4) {
    tree_partition (npe, NULL);
    double f = cut_fraction (&pairs);
    printf ("  %d leaves, depth %d, %d processes: %.3f of %d face pairs cut, "
            "imbalance %.3f\n", tree_caches ()->nleaves, tree_get_depth (), npe,
            f, pairs, tree_get_partition ()->imbalance);
  }

  /* a change of the tree makes it stale; made again, valid. The children
   * of a leaf get its pid */
  {
    TreeItem last = tree_caches ()->leaves [tree_caches ()->nleaves - 1];
    CHECK (last.node->pid > 0 && tree_refine (last.p) > 0);
    int inherit = 1;
    for (int c = 0; c < TREE_NCHILDREN; c++)
      inherit &= index_node (index_child (last.p, c))->pid == last.node->pid;
    CHECK (inherit);
  }
  CHECK (!tree_partition_valid () && tree_get_partition ());
  CHECK (tree_partition (5, NULL) == 0);
  check (5, NULL);

  free (sorted);
  tree_partition_free ();
  CHECK (!tree_get_partition () && !tree_partition_valid ());
  tree_destroy ();
  printf ("dimension %d: %d checks, %d failed\n", dimension, checks, fails);
  return fails != 0;
}
