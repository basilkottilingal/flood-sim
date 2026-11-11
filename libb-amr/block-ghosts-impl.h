/*
.. The core of the ghosts (block-ghosts.h), dispatching to the operators
.. defined before it in the file that includes it, once in a program:
.. block-ghosts.c, or a program's own (see block-ghosts.h).
*/
#include <assert.h>
#include <stdlib.h>
#include <string.h>

#include "block-ghosts.h"

/* ---------------------------------------------------------------------- */
/* the operations of each scalar                                          */
/* ---------------------------------------------------------------------- */

static ScalarOps scalar_ops [NSCALAR_MAX];    /* zero: the defaults */

void scalar_set_ops (scalar s, ScalarOps ops)
{
  assert (s >= 0 && s < NSCALAR_MAX);
  scalar_ops [s] = ops;
  ghosts_mark_dirty (s);      /* made by the old ones: made again */
}

ScalarOps scalar_get_ops (scalar s)
{
  assert (s >= 0 && s < NSCALAR_MAX);
  ScalarOps o = scalar_ops [s];
  if (!o.prolongation)
    o.prolongation = prolong_linear;
  if (!o.restriction)
    o.restriction = restriction_mean;
  if (!o.boundary)
    o.boundary = bc_neumann;
  return o;
}

unsigned long long ghosts_dirty [GHOSTS_DIRTY_WORDS];

/* a new or freed scalar: the default operations, and no stale ghosts */
void scalar_ops_reset (scalar s)
{
  assert (s >= 0 && s < NSCALAR_MAX);
  scalar_ops [s] = (ScalarOps) {0};
  ghosts_dirty [s / 64] &= ~(1ull << (s % 64));
}

/* ---------------------------------------------------------------------- */
/* refining and coarsening one node                                       */
/* ---------------------------------------------------------------------- */

void block_prolong (Index p, int s)
{
  BlockShape sh = tree_block_shape ();
  int m = sh.m;
  assert (index_is_allocated (p) && tree_has_children (index_node (p)));
  ScalarOps o = scalar_get_ops (s);
  Prolongation pr = {
    .coarse = block_data (index_node (p)->s, s),
    .shape = sh, .coarse_index = p, .s = s,
  };
  for (int c = 0; c < TREE_NCHILDREN; c++)
  {
    pr.fine_index = index_child (p, c);
    pr.fine = block_data (index_node (pr.fine_index)->s, s);
    for (int a = 0; a < 3; a++)
    {
      int in = a < dimension;
      pr.off [a] = in ? ((c >> a) & 1) * m : 0;
      pr.lo [a] = 0;
      pr.hi [a] = in ? m : 1;
    }
    @prolongation (o.prolongation) (&pr);
  }
}

void block_restrict (Index p, int s)
{
  @restriction (scalar_get_ops (s).restriction) (p, s);
}

/* ---------------------------------------------------------------------- */
/* the ghosts of one node                                                 */
/* ---------------------------------------------------------------------- */

/*
.. The regions around the block: one for each of the 3^dimension - 1
.. directions, its faces, edges and corners. Those that reach outside the
.. domain are left to the boundary conditions, which come last, one axis
.. after the other: a corner across a wall then mirrors the ghosts along the
.. wall, filled already. Each region inside comes from one node: the one of
.. the same level in that direction, rows copied; or, if there is none,
.. the leaf one level coarser covering it, prolonged.
*/
void block_ghosts (Index p, int s)
{
  BlockShape sh = tree_block_shape ();
  int m = sh.m, g = sh.g, l = p.level, nb = 1 << l;
  assert (index_is_allocated (p));
  Real * own = block_data (index_node (p)->s, s);
  ScalarOps o = scalar_get_ops (s);

  int P [3] = { p.i, 0, 0 }, ndir = 1;
  #if dimension > 1
  P [1] = p.j;
  #endif
  #if dimension > 2
  P [2] = p.k;
  #endif
  for (int a = 0; a < dimension; a++)
    ndir *= 3;

  for (int dir = 0; dir < ndir; dir++)
  {
    int d [3] = { 0, 0, 0 }, q [3] = { 0, 0, 0 };
    int lo [3] = { 0, 0, 0 }, hi [3] = { 1, 1, 1 }, t = dir;
    int centre = 1, outside = 0;
    for (int a = 0; a < dimension; a++)
    {
      d [a] = t % 3 - 1;
      t /= 3;
      centre &= d [a] == 0;
      q [a] = P [a] + d [a];
      outside |= q [a] < 0 || q [a] >= nb;
      lo [a] = d [a] < 0 ? -g : d [a] == 0 ? 0 : m;
      hi [a] = d [a] < 0 ? 0 : d [a] == 0 ? m : m + g;
    }
    if (centre || outside)
      continue;

    Index b = { q [0],
  #if dimension > 1
                q [1],
  #endif
  #if dimension > 2
                q [2],
  #endif
                l };
    if (index_is_allocated (b))
    {
      /* the same level: rows copied, cell x from cell x - d m of b */
      const Real * src = block_data (index_node (b)->s, s);
      size_t len = (size_t) (hi [0] - lo [0]) * sizeof (Real);
      int x [3] = { lo [0], 0, 0 }, y [3];
      for (x [2] = lo [2]; x [2] < hi [2]; x [2]++)
      for (x [1] = lo [1]; x [1] < hi [1]; x [1]++)
      {
        for (int a = 0; a < 3; a++)
          y [a] = x [a] - d [a] * m;
        memcpy (own + block_at (sh, x), src + block_at (sh, y), len);
      }
      continue;
    }

    /* one level coarser, 2:1 balance: the leaf covering b */
    assert (l > 0);
    Index C = b;
    C.level = l - 1;
    C.i >>= 1;
  #if dimension > 1
    C.j >>= 1;
  #endif
  #if dimension > 2
    C.k >>= 1;
  #endif
    assert (index_is_allocated (C));
    int Cq [3] = { C.i, 0, 0 };
  #if dimension > 1
    Cq [1] = C.j;
  #endif
  #if dimension > 2
    Cq [2] = C.k;
  #endif
    Prolongation pr = {
      .coarse = block_data (index_node (C)->s, s), .fine = own, .shape = sh,
      .coarse_index = C, .fine_index = p, .s = s,
    };
    for (int a = 0; a < 3; a++)
    {
      pr.off [a] = a < dimension ? (P [a] - 2 * Cq [a]) * m : 0;
      pr.lo [a] = lo [a];
      pr.hi [a] = hi [a];
    }
    @prolongation (o.prolongation) (&pr);
  }

  /* the walls, one axis after the other */
  for (int a = 0; a < dimension; a++)
    for (int side = 0; side < 2; side++)
      if (side ? P [a] == nb - 1 : P [a] == 0)
      {
        Boundary bc = { own, sh, a, side, p, s };
        @boundary (o.boundary) (&bc);
      }
}

/* ---------------------------------------------------------------------- */
/* every node, level by level                                             */
/* ---------------------------------------------------------------------- */

/*
.. The nodes, by level: at [first [l] .. first [l + 1]) those of level l,
.. leaves only, or parents too. Each in Morton order.
*/
typedef struct
{
  Index * at;
  int first [BAMR_MAX_DEPTH + 2];
  int depth;
} Levels;

static Levels levels_new (int parents)
{
  Levels L = { .depth = tree_get_depth () };
  int count [BAMR_MAX_DEPTH + 1] = {0};
  morton (TREE_ROOT_INDEX)
    if (parents || (node->flags & BAMR_NODE_LEAF))
      count [level]++;
  L.first [0] = 0;
  for (int l = 0; l <= L.depth; l++)
    L.first [l + 1] = L.first [l] + count [l];
  L.at = malloc ((size_t) (L.first [L.depth + 1] + 1) * sizeof (Index));
  if (!L.at)
    abort ();
  int next [BAMR_MAX_DEPTH + 1];
  for (int l = 0; l <= L.depth; l++)
    next [l] = L.first [l];
  morton (TREE_ROOT_INDEX)
    if (parents || (node->flags & BAMR_NODE_LEAF))
      L.at [next [level]++] = (Index) {IJKL};
  return L;
}

/* each parent, finest first: in Morton order a parent comes after its
.. children, so one pass from the leaves up */
void tree_restrict (int s)
{
  int restriction = scalar_get_ops (s).restriction;
  morton (TREE_ROOT_INDEX)
    if (tree_has_children (node))
      @restriction (restriction) ((Index) {IJKL}, s);
}

void ghosts_update_list (int n, const int * s, int parents)
{
  if (n <= 0)
    return;
  Levels L = levels_new (1);
  ScalarOps * o = malloc ((size_t) n * sizeof (ScalarOps));
  if (!o)
    abort ();
  for (int q = 0; q < n; q++)
    o [q] = scalar_get_ops (s [q]);

  /* 1. restriction, finest first; within a level, at once */
  for (int l = L.depth - 1; l >= 0; l--)
  {
    BLOCK_PRAGMA (omp parallel for schedule (dynamic))
    for (int b = L.first [l]; b < L.first [l + 1]; b++)
      if (tree_has_children (index_node (L.at [b])))
        for (int q = 0; q < n; q++)
          @restriction (o [q].restriction) (L.at [b], s [q]);
  }

  /* 2. the ghosts, coarsest first; within a level, at once */
  for (int l = 0; l <= L.depth; l++)
  {
    BLOCK_PRAGMA (omp parallel for schedule (dynamic))
    for (int b = L.first [l]; b < L.first [l + 1]; b++)
      if (parents || (index_node (L.at [b])->flags & BAMR_NODE_LEAF))
        for (int q = 0; q < n; q++)
          block_ghosts (L.at [b], s [q]);
  }
  for (int q = 0; q < n; q++)
    ghosts_dirty [s [q] / 64] &= ~(1ull << (s [q] % 64));
  free (o);
  free (L.at);
}

void ghosts_update_dirty (int n, const int * s)
{
  int list [NSCALAR_MAX], m = 0;
  if (n < 0) {
    foreach_scalar (t)
      if (ghosts_is_dirty (t))
        list [m++] = t;
  }
  else
    for (int q = 0; q < n; q++)
    {
      int t = s [q], seen = 0;
      if (t < 0 || t >= NSCALAR_MAX || !scalar_in_use (t) || !ghosts_is_dirty (t))
        continue;
      for (int r = 0; r < m && !seen; r++)
        seen = list [r] == t;
      if (!seen)
        list [m++] = t;
    }
  ghosts_update_list (m, list, 0);
}

void ghosts_update_scalars (int parents)
{
  int list [NSCALAR_MAX], n = 0;
  foreach_scalar (s)
    list [n++] = s;
  ghosts_update_list (n, list, parents);
}
