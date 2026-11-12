/* The tree: counts, the scalar hooks, Morton order, the loops of tree.h and
 * morton-curve.h, and refinement. Includes tree.c to refine. */
#include "tree.c"
#include "morton-curve.h"
#include <sys/resource.h>
#include <sys/wait.h>
#include <unistd.h>

static int checks, fails;
#define CHECK(c) do { checks++; if (!(c)) { fails++; \
  printf ("FAIL line %d: %s\n", __LINE__, #c); } } while (0)

scalar a, b;

static int ipow (int x, int e) { int r = 1; while (e--) r *= x; return r; }

/* the tree's lists (tree_caches ()) against a walk, node for node, in its
.. order; and by level, each level's nodes in that order */
static int caches_agree (void)
{
  const TreeCaches * tc = tree_caches ();
  int q = 0, ql = 0, bad = 0;
  morton (TREE_ROOT_INDEX) {
    bad += q >= tc->nall || tc->all [q].node != node || tc->all [q].p.level != level;
    q++;
    if (node->flags & BAMR_NODE_LEAF) {
      bad += ql >= tc->nleaves || tc->leaves [ql].node != node;
      ql++;
    }
  }
  bad += q != tc->nall || ql != tc->nleaves || q != tree_get_nblocks ()
         || ql != tree_get_nleaves ();
  for (int l = 0; l <= BAMR_MAX_DEPTH; l++) {
    for (int b = tc->level_first [l]; b < tc->level_first [l + 1]; b++)
      bad += tc->by_level [b].p.level != l;
    for (int b = tc->leaf_first [l]; b < tc->leaf_first [l + 1]; b++)
      bad += tc->leaves_by_level [b].p.level != l
             || !(tc->leaves_by_level [b].node->flags & BAMR_NODE_LEAF);
  }
  bad += tc->level_first [BAMR_MAX_DEPTH + 1] != q
         || tc->leaf_first [BAMR_MAX_DEPTH + 1] != ql;
  return !bad;
}

/* does refining in foreach_leaf () stop the run? in a child process */
static int change_in_loop_stops (void)
{
  fflush (stdout);
  pid_t pid = fork ();
  if (pid == 0) {
    struct rlimit none = { 0, 0 };
    setrlimit (RLIMIT_CORE, &none);
    if (!freopen ("/dev/null", "w", stderr))
      _exit (2);
    int r = 0;
    foreach_leaf ()
      if (r++ == 0)
        tree_refine ((Index) {IJKL});
    _exit (0);                        /* not reached: stopped */
  }
  int status = 0;
  waitpid (pid, &status, 0);
  return WIFSIGNALED (status) && WTERMSIG (status) == SIGABRT;
}

static int nblocks (int s)
{
  int n = 0;
  morton (TREE_ROOT_INDEX) { n += sblocks [s/8] && sblocks [s/8][s%8]; }
  return n;
}

/* the levels holding leaves, by a traversal: what level_bits must be */
static uint32_t scan_bits (void)
{
  uint32_t bits = 0;
  morton (TREE_ROOT_INDEX)
    if (node->flags & BAMR_NODE_LEAF)
      bits |= (uint32_t) 1 << level;
  return bits;
}

int main (void)
{
  int nc = TREE_NCHILDREN;
  int depth = 0, cur = 0, stack = 0, m = 0, c = 0;   /* names the loops use too */
  for (int round = 0; round < 2; round++) {
    tree (4, 4);
    CHECK (tree_get_depth () == 1 && tree_get_nblocks () == 1 + nc && tree_get_nleaves () == nc);
    CHECK (nblocks (a) == 1 + nc && nblocks (b) == 1 + nc && grid_scalar_check ());

    /* the scalar hooks: blocks for a local scalar, and back when it ends */
    {
      scalar t;
      CHECK (nblocks (t) == 1 + nc && grid_scalar_check ());
      (void) t;
    }
    CHECK (grid_scalar_check () && grid_scalar_blocks_attached () == 2 * (1 + nc));

    /* Morton order: the children c = i | j << 1 | k << 2, then the root */
    int seen = 0, ok = 1;
    morton (TREE_ROOT_INDEX) {
      if (level == 1) {
        int cc = i;
#if dimension > 1
        cc |= j << 1;
#endif
#if dimension > 2
        cc |= k << 2;
#endif
        ok &= cc == seen++;
      }
      else
        ok &= level == 0 && seen == nc;
      depth++;
    }
    CHECK (ok && depth == (round + 1) * (1 + nc));

    /* the indices in a cache */
    {
      cache (Index) ix;
      cache_push (ix, (Index) {TREE_ROOT});
      tree_block_indices (&ix);
      CHECK (ix.n == 1 + nc && ix.a [nc].level == 0);
    }

    /* refine the corner chain to level 4 */
    Index p = index_child (TREE_ROOT_INDEX, 0);
    for (int l = 1; l < 4; l++) {
      INDEX_UNPACK (p);
      CHECK (divide_leaf (IJKL));
      p = index_child (p, 0);
    }
    { INDEX_UNPACK (p); CHECK (!divide_leaf (IJKL)); }    /* depthmax */
    CHECK (tree_get_depth () == 4 && tree_get_nblocks () == 1 + 4 * nc);
    CHECK (tree_get_nleaves () == 1 + 4 * (nc - 1) && tree_get_depth_min () == 1);
    CHECK (nblocks (a) == tree_get_nblocks () && grid_scalar_check ());
    CHECK (G (level_bits) == scan_bits () && G (level_bits) == 0x1e);

    /* children_of, neighborhood_of, neighbors_of */
    int n = 0;
    children_of (TREE_ROOT_INDEX) n += node != NULL;
    CHECK (n == nc);
    Index centre = index_child (index_child (TREE_ROOT_INDEX, 0), nc - 1);
    int all = 0, in = 0, ring = 0;
    neighborhood_of (centre, 1) { all++; in += node != NULL; }
    neighbors_of (centre, 1) ring++;
    CHECK (all == ipow (3, dimension) && in == nc && ring == all - 1);
    n = 0;
    neighborhood_of (centre, 1) { if (n == 2) break; n++; }
    CHECK (n == 2);

    /* the tree's lists, after refining: as a walk */
    CHECK (caches_agree ());

    /* coarsening the chain back: depth follows the highest bit */
    Index q = p;
    for (int l = 3; l >= 1; l--) {
      q.level--;
      q.i >>= 1;
#if dimension > 1
      q.j >>= 1;
#endif
#if dimension > 2
      q.k >>= 1;
#endif
      CHECK (tree_coarsen (q) == 1);
      CHECK (G (level_bits) == scan_bits () && tree_get_depth () == l
             && tree_get_depth_min () == 1);
    }
    /* and after coarsening; a change of the tree in a loop over it stops */
    CHECK (caches_agree ());
    CHECK (change_in_loop_stops ());
    tree_destroy ();
    CHECK (G (level_bits) == 0 && tree_get_depth () == 0 && tree_get_depth_min () == 0);
  }
  CHECK (!cur && !stack && !m && !c);
  printf ("dimension %d: %d checks, %d failed\n", dimension, checks, fails);
  return fails != 0;
}
