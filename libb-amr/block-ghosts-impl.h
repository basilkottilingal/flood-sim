/*
.. The core of the ghosts (block-ghosts.h), dispatching to the operators
.. defined before it in the file that includes it, once in a program:
.. block-ghosts.c, or a program's own (see block-ghosts.h).
*/
#include <assert.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "block-ghosts.h"

/* ---------------------------------------------------------------------- */
/* the operations of each scalar                                          */
/* ---------------------------------------------------------------------- */

static ScalarOps scalar_ops [NSCALAR_MAX];    /* zero: the defaults */
unsigned long long ghost_dependents [NSCALAR_MAX][GHOSTS_DIRTY_WORDS];
int ghost_ops_abort = 1;

/* how many scalars the operator op of a kind depends on; -1: no such op */
static int ghost_op_uses (int kind, int op)
{
  if (!op)
    return 0;                 /* the default */
  switch (kind) {
  case GHOST_PROLONGATION: return @prolongation uses (op);
  case GHOST_RESTRICTION:  return @restriction uses (op);
  case GHOST_BOUNDARY:     return @boundary uses (op);
  }
  return -1;
}

/* the scalars that s's operators depend on, into d; how many */
/* the boundary condition of each face, and what it depends on: apart
.. from ScalarOps, which holds the other two */
static int boundary_ops [NSCALAR_MAX][GHOST_FACES];     /* 0: default */
static scalar boundary_uses [NSCALAR_MAX][GHOST_FACES][GHOST_USES_MAX];

static int ghost_uses_of (int s, int * d)
{
  const ScalarOps * o = &scalar_ops [s];
  int n = 0, c;

  if ((c = ghost_op_uses (GHOST_PROLONGATION, o->prolongation)) > 0)
    for (int k = 0; k < c; k++)
      d [n++] = o->prolongation_uses [k];
  if ((c = ghost_op_uses (GHOST_RESTRICTION, o->restriction)) > 0)
    for (int k = 0; k < c; k++)
      d [n++] = o->restriction_uses [k];
  for (int f = 0; f < GHOST_FACES; f++)
    if ((c = ghost_op_uses (GHOST_BOUNDARY, boundary_ops [s][f])) > 0)
      for (int k = 0; k < c; k++)
        d [n++] = boundary_uses [s][f][k];
  return n;
}

/* ghost_dependents again, from every scalar's operators, and through them:
.. over the scalars in use only, as an id not in use has none (reset when
.. made and when freed), and what one depends on is global, so in use */
static void ghost_dependents_rebuild (void)
{
  int d [GHOST_DEPS_MAX];

  memset (ghost_dependents, 0, sizeof ghost_dependents);
  foreach_scalar (s)
    for (int k = 0, n = ghost_uses_of (s, d); k < n; k++)
      ghost_dependents [d [k]][s / 64] |= 1ull << (s % 64);
  for (int changed = 1; changed; ) {
    changed = 0;
    foreach_scalar (t)
      foreach_scalar (s)
        if ((ghost_dependents [t][s / 64] >> (s % 64)) & 1)
          for (int w = 0; w < GHOSTS_DIRTY_WORDS; w++)
            if (ghost_dependents [s][w] & ~ghost_dependents [t][w]) {
              ghost_dependents [t][w] |= ghost_dependents [s][w];
              changed = 1;
            }
  }
}

/* a scalar that depends on itself, or -1 */
static int ghost_dependents_cycle (void)
{
  foreach_scalar (s)
    if ((ghost_dependents [s][s / 64] >> (s % 64)) & 1)
      return s;
  return -1;
}

static int ghost_ops_fail (const char * format, ...)
{
  va_list ap;

  va_start (ap, format);
  vfprintf (stderr, format, ap);
  va_end (ap);
  if (ghost_ops_abort)
    abort ();
  return 1;
}

static const char * const ghost_kinds [] = { "prolongation", "restriction",
                                             "boundary" };

/* what a binding is, for messages: prolongation, boundary.left, ... */
static const char * ghost_what (int kind, int face)
{
  static const char * const bounds [GHOST_FACES] = {
    "boundary.left", "boundary.right",
  #if dimension > 1
    "boundary.bottom", "boundary.top",
  #endif
  #if dimension > 2
    "boundary.back", "boundary.front",
  #endif
  };
  if (kind != GHOST_BOUNDARY)
    return ghost_kinds [kind];
  return face >= 0 && face < GHOST_FACES ? bounds [face] : "boundary.all";
}

/* a scalar in a message: its name (scalar.h), or "#unnamed", which no
.. scalar can be named (# is in no C name), with its id */
#define SNAME(s) (*scalar_name (s) ? scalar_name (s) : "#unnamed")

int ghost_ops_set (scalar s, int kind, int face, int op, const scalar * uses,
                   int n)
{
  assert (s >= 0 && s < NSCALAR_MAX && kind >= 0 && kind <= GHOST_BOUNDARY);
  if (kind == GHOST_BOUNDARY && face != GHOST_ALL_FACES
      && (face < 0 || face >= GHOST_FACES))
    return ghost_ops_fail ("ghost_ops: scalar %s (%d): no face %d in %d "
                           "dimension%s\n", SNAME (s), (int) s, face, dimension,
                           dimension == 1 ? "" : "s");
  int need = ghost_op_uses (kind, op);
  if (need < 0)
    return ghost_ops_fail ("ghost_ops: scalar %s (%d): no %s with the id %d\n",
                           SNAME (s), (int) s, ghost_kinds [kind], op);
  if (n != need)
    return ghost_ops_fail ("ghost_ops: scalar %s (%d): its %s %d depends on %d "
                           "scalar%s, not %d\n", SNAME (s), (int) s,
                           ghost_what (kind, face), op, need, need == 1 ? "" : "s",
                           n);
  if (n > GHOST_USES_MAX)
    return ghost_ops_fail ("ghost_ops: scalar %s (%d): %d scalars depended on, "
                           "%d at most (GHOST_USES_MAX)\n", SNAME (s), (int) s,
                           n, GHOST_USES_MAX);
  for (int k = 0; k < n; k++)
    if (uses [k] < 0 || uses [k] >= n_global_scalars)
      return ghost_ops_fail ("ghost_ops: scalar %s (%d): scalar %s (%d), which "
                             "its %s depends on, is not global\n", SNAME (s),
                             (int) s, SNAME (uses [k]), (int) uses [k],
                             ghost_what (kind, face));
  /* kept, to undo it all on a cycle */
  ScalarOps old = scalar_ops [s];
  int old_ops [GHOST_FACES];
  scalar old_uses [GHOST_FACES][GHOST_USES_MAX];
  memcpy (old_ops, boundary_ops [s], sizeof old_ops);
  memcpy (old_uses, boundary_uses [s], sizeof old_uses);
  if (kind == GHOST_BOUNDARY) {
    for (int f = 0; f < GHOST_FACES; f++)
      if (face == GHOST_ALL_FACES || face == f) {
        boundary_ops [s][f] = op;
        for (int k = 0; k < GHOST_USES_MAX; k++)
          boundary_uses [s][f][k] = k < n ? uses [k] : 0;
      }
  }
  else {
    int * id = kind == GHOST_PROLONGATION ? &scalar_ops [s].prolongation
                                          : &scalar_ops [s].restriction;
    scalar * to = kind == GHOST_PROLONGATION ? scalar_ops [s].prolongation_uses
                                             : scalar_ops [s].restriction_uses;
    *id = op;
    for (int k = 0; k < GHOST_USES_MAX; k++)
      to [k] = k < n ? uses [k] : 0;
  }
  ghost_dependents_rebuild ();
  int c = ghost_dependents_cycle ();
  if (c >= 0) {
    scalar_ops [s] = old;
    memcpy (boundary_ops [s], old_ops, sizeof old_ops);
    memcpy (boundary_uses [s], old_uses, sizeof old_uses);
    ghost_dependents_rebuild ();
    return ghost_ops_fail ("ghost_ops: scalar %s (%d): its operators would "
                           "depend on itself, through scalar %s (%d)\n",
                           SNAME (s), (int) s, SNAME (c), c);
  }
  ghosts_mark_dirty (s);      /* made by the old ones: made again */
  return 0;
}

void scalar_set_ops (scalar s, ScalarOps ops)
{
  assert (s >= 0 && s < NSCALAR_MAX);
  if (ghost_op_uses (GHOST_PROLONGATION, ops.prolongation) > 0
      || ghost_op_uses (GHOST_RESTRICTION, ops.restriction) > 0) {
    ghost_ops_fail ("scalar_set_ops: scalar %s (%d): an operator that depends "
                    "on scalars is bound with ghost_ops (), which names them\n",
                    SNAME (s), (int) s);
    return;
  }
  for (int k = 0; k < GHOST_USES_MAX; k++)
    ops.prolongation_uses [k] = ops.restriction_uses [k] = 0;
  scalar_ops [s] = ops;
  ghost_dependents_rebuild ();
  ghosts_mark_dirty (s);      /* made by the old ones: made again */
}

/* each scalar after those it depends on, the stale ones of those added */
static void ghost_order_visit (int s, char * state, const char * want,
                               int * out, int * n)
{
  int d [GHOST_DEPS_MAX];

  if (state [s])
    return;
  state [s] = 1;
  for (int k = 0, nd = ghost_uses_of (s, d); k < nd; k++)
    if (want [d [k]] || ghosts_is_dirty (d [k]))
      ghost_order_visit (d [k], state, want, out, n);
  out [(*n)++] = s;
}

int ghost_order (int n, const int * in, int * out)
{
  char state [NSCALAR_MAX] = {0}, want [NSCALAR_MAX] = {0};
  int m = 0;

  for (int q = 0; q < n; q++)
    if (in [q] >= 0 && in [q] < NSCALAR_MAX)
      want [in [q]] = 1;
  for (int q = 0; q < n; q++)
    if (in [q] >= 0 && in [q] < NSCALAR_MAX)
      ghost_order_visit (in [q], state, want, out, &m);
  return m;
}

int ghost_order_all (int * out)
{
  int in [NSCALAR_MAX], n = 0;

  foreach_scalar (t)
    in [n++] = t;
  return ghost_order (n, in, out);
}

ScalarOps scalar_get_ops (scalar s)
{
  assert (s >= 0 && s < NSCALAR_MAX);
  ScalarOps o = scalar_ops [s];
  if (!o.prolongation)
    o.prolongation = prolong_linear;
  if (!o.restriction)
    o.restriction = restriction_mean;
  return o;
}

int scalar_get_boundary (scalar s, int face)
{
  assert (s >= 0 && s < NSCALAR_MAX && face >= 0 && face < GHOST_FACES);
  return boundary_ops [s][face] ? boundary_ops [s][face] : bc_neumann;
}

unsigned long long ghosts_dirty [GHOSTS_DIRTY_WORDS];

/* a new or freed scalar: the default operations, and no stale ghosts */
void scalar_ops_reset (scalar s)
{
  assert (s >= 0 && s < NSCALAR_MAX);
  scalar_ops [s] = (ScalarOps) {0};
  memset (boundary_ops [s], 0, sizeof boundary_ops [s]);
  memset (boundary_uses [s], 0, sizeof boundary_uses [s]);
  ghost_dependents_rebuild ();
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
    pr.off [0] = (c & 1) * m, pr.lo [0] = 0, pr.hi [0] = m;
    pr.off [1] = pr.off [2] = pr.lo [1] = pr.lo [2] = 0;
    pr.hi [1] = pr.hi [2] = 1;
  #if dimension > 1
    pr.off [1] = ((c >> 1) & 1) * m, pr.hi [1] = m;
  #endif
  #if dimension > 2
    pr.off [2] = ((c >> 2) & 1) * m, pr.hi [2] = m;
  #endif
    @prolongation (o.prolongation, o.prolongation_uses) (&pr);
  }
}

/*
.. p's interior as its parent's prolongation makes it, into target, p's
.. block left as it is. For the call, target is p's block of s: an operator
.. writes there, through pr->fine or s[], while coarse (s) reads the parent
.. and other scalars p's own blocks. A prolongation only writes its fine
.. cells, never reads them (it fills uninitialised blocks in divide_leaf),
.. so the cells of p's block it would have read are none.
*/
void block_prolong_into (Index p, int s, Real * target)
{
  BlockShape sh = tree_block_shape ();
  int m = sh.m;
  assert (p.level > 0 && index_is_allocated (p) && target);
  Index C = p;
  C.level--;
  C.i >>= 1;
  #if dimension > 1
  C.j >>= 1;
  #endif
  #if dimension > 2
  C.k >>= 1;
  #endif
  ScalarOps o = scalar_get_ops (s);
  void ** slot = &index_node (p)->s [s / SCALAR_ARRAY_WIDTH][s % SCALAR_ARRAY_WIDTH];
  void * own = *slot;
  *slot = target;
  Prolongation pr = {
    .coarse = block_data (index_node (C)->s, s), .fine = target, .shape = sh,
    .coarse_index = C, .fine_index = p, .s = s,
  };
  pr.off [0] = (p.i & 1) * m, pr.lo [0] = 0, pr.hi [0] = m;
  pr.off [1] = pr.off [2] = pr.lo [1] = pr.lo [2] = 0;
  pr.hi [1] = pr.hi [2] = 1;
  #if dimension > 1
  pr.off [1] = (p.j & 1) * m, pr.hi [1] = m;
  #endif
  #if dimension > 2
  pr.off [2] = (p.k & 1) * m, pr.hi [2] = m;
  #endif
  @prolongation (o.prolongation, o.prolongation_uses) (&pr);
  *slot = own;
}

void block_restrict (Index p, int s)
{
  ScalarOps o = scalar_get_ops (s);
  @restriction (o.restriction, o.restriction_uses) (p, s);
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

  int P [3] = { p.i, 0, 0 };
  #if dimension > 1
  P [1] = p.j;
  #endif
  #if dimension > 2
  P [2] = p.k;
  #endif
  enum { ndir = dimension == 1 ? 3 : dimension == 2 ? 9 : 27 };   /* 3^dim */

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
      int x [3] = { lo [0], 0, 0 }, y [3] = { lo [0] - d [0] * m, 0, 0 };
  #if dimension > 2
      for (x [2] = lo [2]; x [2] < hi [2]; x [2]++)
  #endif
  #if dimension > 1
      for (x [1] = lo [1]; x [1] < hi [1]; x [1]++)
  #endif
      {
  #if dimension > 1
        y [1] = x [1] - d [1] * m;
  #endif
  #if dimension > 2
        y [2] = x [2] - d [2] * m;
  #endif
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
      pr.off [a] = (P [a] - 2 * Cq [a]) * m;   /* 0 beyond the dimension */
      pr.lo [a] = lo [a];
      pr.hi [a] = hi [a];
    }
    @prolongation (o.prolongation, o.prolongation_uses) (&pr);
  }

  /* the walls, one axis after the other */
  for (int a = 0; a < dimension; a++)
    for (int side = 0; side < 2; side++)
      if (side ? P [a] == nb - 1 : P [a] == 0)
      {
        Boundary bc = { own, sh, a, side, p, s };
        int f = 2 * a + side;             /* the face: GHOST_LEFT, ... */
        @boundary (scalar_get_boundary (s, f), boundary_uses [s][f]) (&bc);
      }
}

/* ---------------------------------------------------------------------- */
/* every node, level by level                                             */
/* ---------------------------------------------------------------------- */

/* each parent, finest first: in Morton order a parent comes after its
.. children, so one pass over the tree's list from the leaves up */
void tree_restrict (int s)
{
  ScalarOps o = scalar_get_ops (s);
  const TreeCaches * tc = tree_caches ();
  for (int q = 0; q < tc->nall; q++)
    if (tree_has_children (tc->all [q].node))
      @restriction (o.restriction, o.restriction_uses) (tc->all [q].p, s);
}

/*
.. The ghosts of the grid at depth l, foreach_level_or_leaf (l), for s
.. alone: what the levels of a multigrid cycle need, without the
.. restriction of the leaves that ghosts_update () does first, as on level
.. l s is a field of its own (a correction).
..   1. The parents of levels depth_min .. l - 1, finest first, restricted
..      from the grid: a coarse leaf of the grid reads its neighbours of its
..      own level, which may be parents. Below depth_min no leaf is coarser
..      than l, and nothing is restricted. Their values are overwritten.
..   2. The ghosts of the grid's leaves coarser than l, coarsest first, as
..      a finer one at a level jump prolongs from their first ghost layer;
..      then of its nodes of level l.
.. Then s is marked stale: its leaves' ghosts are the grid's now, and a
.. foreach () reading them makes them again. What s's operators depend on
.. is not updated: s's alone.
*/
void ghosts_update_level (int s, int l)
{
  const TreeCaches * tc = tree_caches ();
  int L = block_level_clamp (l);
  if (L < 0)
    return;
  for (int k = L - 1; k >= tree_get_depth_min (); k--)
    for (int q = tc->level_first [k]; q < tc->level_first [k + 1]; q++)
      if (tree_has_children (tc->by_level [q].node))
        block_restrict (tc->by_level [q].p, s);
  for (int q = tc->leaf_first [0]; q < tc->leaf_first [L]; q++)
    block_ghosts (tc->leaves_by_level [q].p, s);
  for (int q = tc->level_first [L]; q < tc->level_first [L + 1]; q++)
    block_ghosts (tc->by_level [q].p, s);
  ghosts_mark_dirty (s);
}

void ghosts_update_list (int n, const int * given, int parents)
{
  if (n <= 0)
    return;
  /* each after those it depends on, the stale ones of those with them */
  int s [NSCALAR_MAX];
  n = ghost_order (n, given, s);
  /* the nodes by level, each in Morton order: the tree's list */
  const TreeCaches * tc = tree_caches ();
  const TreeItem * at = tc->by_level;
  const int * first = tc->level_first;
  ScalarOps * o = malloc ((size_t) n * sizeof (ScalarOps));
  if (!o)
    abort ();
  for (int q = 0; q < n; q++)
    o [q] = scalar_get_ops (s [q]);

  /* 1. restriction, finest first; within a level, at once */
  for (int l = tc->depth - 1; l >= 0; l--)
  {
    BLOCK_PRAGMA (omp parallel for schedule (dynamic))
    for (int b = first [l]; b < first [l + 1]; b++)
      if (tree_has_children (at [b].node))
        for (int q = 0; q < n; q++)
          @restriction (o [q].restriction, o [q].restriction_uses) (at [b].p, s [q]);
  }

  /* 2. the ghosts, coarsest first; within a level, at once */
  for (int l = 0; l <= tc->depth; l++)
  {
    BLOCK_PRAGMA (omp parallel for schedule (dynamic))
    for (int b = first [l]; b < first [l + 1]; b++)
      if (parents || (at [b].node->flags & BAMR_NODE_LEAF))
        for (int q = 0; q < n; q++)
          block_ghosts (at [b].p, s [q]);
  }
  for (int q = 0; q < n; q++)
    ghosts_dirty [s [q] / 64] &= ~(1ull << (s [q] % 64));
  free (o);
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
