/* tree () on several processes (tree_mpi_start ()): each process makes its
 * part of a uniform level alone, from the cut's formula. Against the old
 * way, every process holding the whole level (tree_replicated (), refined
 * uniformly, then cut and pruned by tree_mpi_partition ()): every node
 * held, its index, flags, pid and blocks, in order; the partition (first
 * leaves, first keys, loads); the plan, list for list and entry for entry,
 * and the peers. Every block of both 0; every leaf owned once over the
 * processes. Levels 1 to a few (the old way needs the whole level on every
 * process) and tree ()'s default; with a scalar made before tree () and
 * one after. Build with -D_MPI, run on 1, 2, 3 ... processes. */
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>
#include "tree-mpi.h"
#include "block-ghosts.h"
#include "grid-scalar.h"

scalar a;                        /* made before tree () */

static int checks, fails, rank, npe;
#define CHECK(c) do { checks++; if (!(c)) { fails++; \
  printf ("rank %d: FAIL line %d: %s\n", rank, __LINE__, #c); } } while (0)

#define KEPT (BAMR_NODE_LEAF | BAMR_NODE_IS_PRUNED | BAMR_NODE_HAS_CHILDREN \
              | BAMR_NODE_IS_INACTIVE)

typedef struct { Index p; int flags, pid, blocks; } Node;
typedef struct { int parents, kind, level, send, peer; Index p; } Entry;
typedef struct
{
  int n, nleaves, nent, npeers, depth, depth_min, nowned, tdepth, tmin;
  Node * node;
  int * first;
  uint64_t * key;
  double * load;
  Entry * ent;
  int * peers;
  long nonzero;                  /* cells not 0, over every block */
} State;

static int same_index (Index a, Index b)
{
  return !memcmp (&a, &b, sizeof (Index));
}

/* the block of scalar s in table t, NULL if none */
static Real * blk (void *** t, int s)
{
  void ** g = t [s / SCALAR_ARRAY_WIDTH];
  return g ? g [s % SCALAR_ARRAY_WIDTH] : NULL;
}

/* what this process holds of the tree, the partition and the plan */
static State record (scalar b)
{
  State st;
  memset (&st, 0, sizeof (st));
  const TreeCaches * tc = tree_caches ();
  long cells = tree_block_shape ().cells;
  st.n = tc->nall, st.nowned = tc->nowned;
  st.tdepth = tree_get_depth (), st.tmin = tree_get_depth_min ();
  st.node = malloc ((size_t) st.n * sizeof (Node));
  for (int q = 0; q < st.n; q++)
  {
    BAMRNode * node = tc->all [q].node;
    Real * x = blk (node->s, a), * y = blk (node->s, b);
    st.node [q] = (Node) { tc->all [q].p, node->flags & KEPT, node->pid,
                           (x != NULL) | (y != NULL) << 1 };
    for (int k = 0; k < 2; k++)
    {
      Real * z = k ? y : x;
      for (long c = 0; z && c < cells; c++)
        st.nonzero += z [c] != 0;
    }
  }
  const TreePartition * P = tree_get_partition ();
  if (P)
  {
    st.nleaves = P->nleaves;
    st.first = malloc ((size_t) (npe + 1) * sizeof (int));
    st.key = malloc ((size_t) (npe + 1) * sizeof (uint64_t));
    st.load = malloc ((size_t) npe * sizeof (double));
    memcpy (st.first, P->leaf_first, (size_t) (npe + 1) * sizeof (int));
    memcpy (st.key, P->key_first, (size_t) (npe + 1) * sizeof (uint64_t));
    memcpy (st.load, P->load, (size_t) npe * sizeof (double));
  }
  const TreeMpiPlan * plan = tree_mpi_plan ();
  if (plan)
  {
    st.depth = plan->depth, st.depth_min = plan->depth_min;
    st.npeers = plan->npeers;
    st.peers = malloc ((size_t) (st.npeers + 1) * sizeof (int));
    memcpy (st.peers, plan->peers, (size_t) st.npeers * sizeof (int));
    int cap = 1024;
    st.ent = malloc ((size_t) cap * sizeof (Entry));
    for (int par = 0; par < 2; par++)
      for (int kind = 0; kind < TREE_MPI_KINDS; kind++)
        for (int l = 0; l < BAMR_MAX_DEPTH + 2; l++)
          for (int send = 0; send < 2; send++)
          {
            const TreeMpiExchange * ex = &plan->ex [par][kind][l];
            int nl = send ? ex->nsend : ex->nrecv;
            const TreeMpiList * li = send ? ex->send : ex->recv;
            for (int k = 0; k < nl; k++)
              for (int j = 0; j < li [k].n; j++)
              {
                if (st.nent == cap)
                  st.ent = realloc (st.ent, (size_t) (cap *= 2)
                                    * sizeof (Entry));
                st.ent [st.nent++] = (Entry) { par, kind, l, send,
                                               li [k].peer, li [k].item [j].p };
              }
          }
  }
  return st;
}

static void state_free (State * st)
{
  free (st->node), free (st->first), free (st->key), free (st->load);
  free (st->ent), free (st->peers);
}

/* how many nodes, partition entries and plan entries differ */
static void compare (const State * A, const State * B)
{
  long nodes = 0, part = 0, plan = 0;
  CHECK (A->n == B->n);
  for (int q = 0; q < A->n && q < B->n; q++)
    nodes += !same_index (A->node [q].p, B->node [q].p)
          || A->node [q].flags != B->node [q].flags
          || A->node [q].pid != B->node [q].pid
          || A->node [q].blocks != B->node [q].blocks;
  CHECK (nodes == 0);
  CHECK (A->nowned == B->nowned);
  CHECK (A->tdepth == B->tdepth && A->tmin == B->tmin);
  CHECK ((A->first != NULL) == (B->first != NULL));
  if (A->first && B->first)
  {
    part += A->nleaves != B->nleaves;
    for (int r = 0; r <= npe; r++)
      part += A->first [r] != B->first [r] || A->key [r] != B->key [r]
           || (r < npe && A->load [r] != B->load [r]);
  }
  CHECK (part == 0);
  CHECK (A->nent == B->nent && A->npeers == B->npeers
         && A->depth == B->depth && A->depth_min == B->depth_min);
  for (int k = 0; k < A->nent && k < B->nent; k++)
  {
    const Entry * x = &A->ent [k], * y = &B->ent [k];
    plan += x->parents != y->parents || x->kind != y->kind
         || x->level != y->level || x->send != y->send || x->peer != y->peer
         || !same_index (x->p, y->p);
  }
  for (int k = 0; k < A->npeers && k < B->npeers; k++)
    plan += A->peers [k] != B->peers [k];
  CHECK (plan == 0);
  if (nodes || part || plan)
    printf ("rank %d: %ld nodes, %ld partition entries, %ld plan entries "
            "differ (%d against %d nodes)\n", rank, nodes, part, plan, A->n,
            B->n);
}

/* tree ()'s default level: the first with 4 npe families at least */
static int default_level (int depthmax)
{
  int L = 1;
  if (npe > 1)
    while (L < depthmax && (1L << (dimension * (L - 1))) < 4L * npe)
      L++;
  return L;
}

int main (void)
{
  int top = dimension == 1 ? 8 : dimension == 2 ? 5 : 3, depthmax = top;
  long total = 0;
  for (int L = 0; L <= top + 1; L++)       /* top + 1: clamped */
  {
    /* the new way: tree () alone */
    tree_initial_depth (L);
    tree (depthmax, 4);
    rank = tree_rank (), npe = tree_npe ();
    int level = L ? L : default_level (depthmax);
    level = level < depthmax ? level : depthmax;     /* clamped */
    scalar b = scalar_new ();
    State A = record (b);
    if (A.nowned)
      CHECK (tree_get_depth () == level && tree_get_depth_min () == level);
    CHECK (grid_scalar_check ());
    long own = A.nowned, all;
    MPI_Allreduce (&own, &all, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
    CHECK (all == 1L << (dimension * level));
    scalar_free (b);
    tree_destroy ();

    /* the old way: the whole level on every process, then cut */
    tree_replicated (depthmax, 4);
    for (int l = 1; l < level; l++)
    {
      const TreeCaches * tc = tree_caches ();
      int n = tc->nleaves;
      Index * leaves = malloc ((size_t) n * sizeof (Index));
      for (int q = 0; q < n; q++)
        leaves [q] = tc->leaves [q].p;
      for (int q = 0; q < n; q++)
        tree_refine (leaves [q]);
      free (leaves);
    }
    if (npe > 1)
      CHECK (tree_mpi_partition (NULL) == 0);
    scalar c = scalar_new ();
    State B = record (c);
    scalar_free (c);
    tree_destroy ();

    compare (&A, &B);
    CHECK (A.nonzero == 0 && B.nonzero == 0);
    total += A.n;
    state_free (&A), state_free (&B);
  }
  tree_initial_depth (0);

  int ftot = 0;
  MPI_Allreduce (&fails, &ftot, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  printf ("dimension %d, %d processes: levels 0 (default) to %d, %ld nodes "
          "held here; %d checks, %d failed\n", dimension, npe, top, total,
          checks, fails);
  return ftot != 0;
}
