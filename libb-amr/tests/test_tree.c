/* The tree: counts, the scalar hooks, Morton order, the loops of tree.h and
 * morton-curve.h, and refinement. Includes tree.c to refine. */
#include "tree.c"
#include "morton-curve.h"

static int checks, fails;
#define CHECK(c) do { checks++; if (!(c)) { fails++; \
  printf ("FAIL line %d: %s\n", __LINE__, #c); } } while (0)

scalar a, b;

static int ipow (int x, int e) { int r = 1; while (e--) r *= x; return r; }

static int nblocks (int s)
{
  int n = 0;
  morton (TREE_ROOT_INDEX) { n += sblocks [s/8] && sblocks [s/8][s%8]; }
  return n;
}

int main (void)
{
  int nc = TREE_NCHILDREN;
  int depth = 0, cur = 0, stack = 0, m = 0, c = 0;   /* names the loops use too */
  for (int round = 0; round < 2; round++) {
    tree (4, 4);
    CHECK (tree_depth () == 1 && tree_nblocks () == 1 + nc && tree_nleaves () == nc);
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
    CHECK (tree_depth () == 4 && tree_nblocks () == 1 + 4 * nc);
    CHECK (tree_nleaves () == 1 + 4 * (nc - 1) && G (depthmin) == 1);
    CHECK (nblocks (a) == tree_nblocks () && grid_scalar_check ());

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
    tree_destroy ();
  }
  CHECK (!cur && !stack && !m && !c);
  printf ("dimension %d: %d checks, %d failed\n", dimension, checks, fails);
  return fails != 0;
}
