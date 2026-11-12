/* Refining and coarsening by families, 2:1 balance, prolongation and
 * restriction, and adapt_wavelet () */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <cwrap/cache.h>
#include "adapt.h"
#include "block-ghosts.h"

static int checks, fails;
#define CHECK(c) do { checks++; if (!(c)) { fails++; \
  printf ("FAIL line %d: %s\n", __LINE__, #c); } } while (0)

scalar u;

static long ipow (long x, int e) { long r = 1; while (e--) r *= x; return r; }

/* the centre of interior cell (ci, cj, ck) of block p (0 to m - 1), in [0, 1] */
static void centre (Index p, int ci, int cj, int ck, double x [3])
{
  double n = (double) (tree_get_block_m () << p.level);
  int b [3] = { p.i, 0, 0 }, c [3] = { ci, cj, ck };
#if dimension > 1
  b [1] = p.j;
#endif
#if dimension > 2
  b [2] = p.k;
#endif
  for (int a = 0; a < 3; a++)
    x [a] = a < dimension ? (b [a] * tree_get_block_m () + c [a] + 0.5) / n : 0.5;
}

/* u = f on every leaf */
static void set (double (* f) (const double [3]))
{
  BlockShape st = tree_block_shape ();
  foreach_leaf () {
    Real * b = block_data (sblocks, u);
    block_cells (st, 0) {
      double x [3];
      int y = 0, z = 0;
#if dimension > 1
      y = cj - st.g;
#endif
#if dimension > 2
      z = ck - st.g;
#endif
      centre ((Index) {IJKL}, ci - st.ox, y, z, x);
      b [c] = (Real) f (x);
    }
  }
}

static double linear (const double x [3]) { return 1 + 2 * x [0] + 3 * x [1] + 5 * x [2]; }
static double one (const double x [3]) { (void) x; return 1; }
static double front (const double x [3])
{
  double r = 0;
  for (int a = 0; a < dimension; a++)
    r += (x [a] - 0.5) * (x [a] - 0.5);
  return tanh ((sqrt (r) - 0.3) / 0.02);
}

/* the integral of u over the leaves */
static double integral (void)
{
  BlockShape st = tree_block_shape ();
  double sum = 0;
  foreach_leaf () {
    const Real * b = block_data (sblocks, u);
    double v = pow (1. / (tree_get_block_m () << level), dimension);
    block_cells (st, 0) sum += b [c] * v;
  }
  return sum;
}

/* is q inside the domain? */
static int inside (Index q)
{
  int n = 1 << q.level, in = q.i >= 0 && q.i < n;
#if dimension > 1
  in &= q.j >= 0 && q.j < n;
#endif
#if dimension > 2
  in &= q.k >= 0 && q.k < n;
#endif
  return in;
}

/* 2:1: around every leaf, each index of its level holds a node, or its
 * parent index does; and the counts agree with the tree */
static int balanced (void)
{
  int ok = 1, leaves = 0, nodes = 0;
  morton (TREE_ROOT_INDEX) {
    nodes++;
    if (!(node->flags & BAMR_NODE_LEAF))
      continue;
    leaves++;
    Index p = {IJKL};
    neighbors_of (p, 1) {
      Index q = {IJKL};
      if (!inside (q) || node)
        continue;
      Index up = q;
      up.level--;
      up.i >>= 1;
#if dimension > 1
      up.j >>= 1;
#endif
#if dimension > 2
      up.k >>= 1;
#endif
      ok &= index_is_allocated (up);
    }
  }
  return ok && leaves == tree_get_nleaves () && nodes == tree_get_nblocks ();
}

/* the leaf holding the point x */
static Index leaf_at (const double x [3])
{
  for (int l = 0; ; l++) {
    int n = 1 << l;
    Index q = { (int) (x [0] * n),
#if dimension > 1
                (int) (x [1] * n),
#endif
#if dimension > 2
                (int) (x [2] * n),
#endif
                l };
    if (index_is_allocated (q) && (index_node (q)->flags & BAMR_NODE_LEAF))
      return q;
  }
}

/* refine every leaf once */
static void refine_all (void)
{
  cache (Index) all;
  foreach_leaf () cache_push (all, (Index) {IJKL});
  foreach (p in all) tree_refine (p);
}

int main (void)
{
  int nc = TREE_NCHILDREN;
  long cells;
  tree (6, 4);
  cells = tree_block_shape ().cells;

  /* 1. a family is allocated at once: nodes, and the blocks of a scalar */
  {
    BAMRNode * first = index_node (index_child (TREE_ROOT_INDEX, 0));
    Real * b0 = block_data (first->s, u);
    int ok = 1;
    for (int c = 0; c < nc; c++) {
      BAMRNode * n = index_node (index_child (TREE_ROOT_INDEX, c));
      ok &= n == first + c && block_data (n->s, u) == b0 + c * cells;
    }
    CHECK (ok);
  }

  /* 2. prolongation: exact for a linear field, away from the walls, and
   * keeping the mean everywhere */
  refine_all ();                                   /* uniform level 2 */
  CHECK (balanced () && tree_get_nleaves () == ipow (nc, 2));
  set (linear);
  ghosts_update (u);
  Index p = { 1,
#if dimension > 1
              1,
#endif
#if dimension > 2
              1,
#endif
              2 };
  Real * before = malloc (cells * sizeof (Real));
  memcpy (before, block_data (index_node (p)->s, u), cells * sizeof (Real));
  CHECK (tree_refine (p) == 1);
  {
    double err = 0;
    BlockShape st = tree_block_shape ();
    children_of (p) {
      Real * b = block_data (sblocks, u);
      Index q = {IJKL};
      block_cells (st, 0) {
        double x [3];
        int y = 0, z = 0;
#if dimension > 1
        y = cj - st.g;
#endif
#if dimension > 2
        z = ck - st.g;
#endif
        centre (q, ci - st.ox, y, z, x);
        err = fmax (err, fabs (b [c] - linear (x)));
      }
    }
    CHECK (err < 1e-5);
    /* the mean of the children is the parent as it was */
    block_restrict (p, u);
    double d = 0;
    const Real * now = block_data (index_node (p)->s, u);
    for (long q = 0; q < cells; q++)
      d = fmax (d, fabs (now [q] - before [q]));
    CHECK (d < 1e-5);
  }
  free (before);

  /* 3. coarsening gives the parent the mean, and frees the family */
  CHECK (tree_coarsen (p) == 1);
  CHECK ((index_node (p)->flags & BAMR_NODE_LEAF) && balanced () && grid_scalar_check ());
  CHECK (tree_coarsen (p) == 0);                  /* a leaf now */

  /* 4. balance: refining one corner deep refines around it; coarsening is
   * refused where it would break the balance */
  {
    double x [3] = { 0.01, 0.01, 0.01 };
    for (int l = 0; l < 4; l++)
      tree_refine (leaf_at (x));
    CHECK (leaf_at (x).level == 6 && balanced () && grid_scalar_check ());
    /* the parent of the deepest family: its neighbours have grandchildren */
    Index deep = leaf_at (x);
    deep.level -= 3;
    deep.i >>= 3;
#if dimension > 1
    deep.j >>= 3;
#endif
#if dimension > 2
    deep.k >>= 3;
#endif
    CHECK (tree_coarsen (deep) == 0);
    /* random refining and coarsening keep it balanced */
    srand (7);
    int ok = 1;
    for (int it = 0; it < 40; it++) {
      double y [3] = { rand () / (double) RAND_MAX, rand () / (double) RAND_MAX,
                       rand () / (double) RAND_MAX };
      Index q = leaf_at (y);
      if (it % 3 == 2 && q.level > 1) {
        q.level--;
        q.i >>= 1;
#if dimension > 1
        q.j >>= 1;
#endif
#if dimension > 2
        q.k >>= 1;
#endif
        tree_coarsen (q);
      }
      else
        tree_refine (q);
      ok &= balanced ();
    }
    CHECK (ok && grid_scalar_check ());
  }

  /* 5. the wavelet of a linear field is zero away from the walls */
  set (linear);
  ghosts_update (u);
  {
    double x [3] = { 0.5, 0.5, 0.5 };
    Index q = leaf_at (x);
    CHECK (block_wavelet (q, u) < 1e-5);
  }

  /* 6. adapt_wavelet on a front: fine on it, coarse away from it, balanced,
   * and the integral kept by each adaptation */
  {
    double max = 1e-3;
    int its = 0;
    AdaptStats st;
    do {
      set (front);
      double before = integral ();
      st = adapt_wavelet (1, (int []) { u }, &max, 1, 4);
      CHECK (fabs (integral () - before) < 1e-6);
      CHECK (balanced () && grid_scalar_check ());
    } while ((st.refined || st.coarsened) && ++its < 20);
    CHECK (its < 20);                              /* it settles */
    double on [3] = { 0.8, 0.5, 0.5 }, off [3] = { 0.02, 0.02, 0.02 };
    CHECK (leaf_at (on).level == 4 && tree_get_depth () == 4);   /* maxlevel */
    CHECK (leaf_at (off).level < 4);
    printf ("  front: %d leaves after %d adaptations, levels %d to %d\n",
            tree_get_nleaves (), its, tree_get_depth_min (), tree_get_depth ());
  }

  /* 7. a constant field coarsens down to minlevel */
  {
    double max = 0.01;
    for (int it = 0; it < 10; it++) {
      set (one);
      adapt_wavelet (1, (int []) { u }, &max, 1, 4);
    }
    CHECK (tree_get_nleaves () == nc && tree_get_depth () == 1 && balanced () && grid_scalar_check ());
  }

  /* 8. the root has no wavelet: even with minlevel 0, no family is
   * coarsened into it; a root made a leaf by hand is refined again */
  {
    double max = 0.01;
    set (one);
    adapt_wavelet (1, (int []) { u }, &max, 0, 4);
    CHECK (tree_get_nleaves () == nc && tree_get_depth () == 1);
    CHECK (tree_coarsen (TREE_ROOT_INDEX) && tree_get_nleaves () == 1);
    set (one);
    AdaptStats st = adapt_wavelet (1, (int []) { u }, &max, 0, 4);
    CHECK (st.refined == 1 && tree_get_nleaves () == nc && grid_scalar_check ());
  }

  tree_destroy ();
  printf ("dimension %d: %d checks, %d failed\n", dimension, checks, fails);
  return fails != 0;
}
