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
.. made and when freed), and what one depends on is in use: a scalar others
.. depend on is not freed (ghost_ops_freeing ()) */
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

/* the types an operator depends on, a letter each (cwrap's @prolongation
.. signature (id)), and its name; "" and "the default" for 0 */
static const char * ghost_op_signature (int kind, int op)
{
  if (!op)
    return "";
  switch (kind) {
  case GHOST_PROLONGATION: return @prolongation signature (op);
  case GHOST_RESTRICTION:  return @restriction signature (op);
  case GHOST_BOUNDARY:     return @boundary signature (op);
  }
  return NULL;
}

static const char * ghost_op_name (int kind, int op)
{
  const char * n = NULL;
  if (!op)
    return "the default";
  switch (kind) {
  case GHOST_PROLONGATION: n = @prolongation name (op); break;
  case GHOST_RESTRICTION:  n = @restriction name (op); break;
  case GHOST_BOUNDARY:     n = @boundary name (op); break;
  }
  return n ? n : "?";
}

/* a signature in words, (vector, scalar), into buf */
static const char * ghost_sig_text (const char * sig, char * buf, size_t size)
{
  size_t k = 0;
  k += (size_t) snprintf (buf + k, size - k, "(");
  for (const char * c = sig; *c && k < size; c++)
    k += (size_t) snprintf (buf + k, size - k, "%s%s", c == sig ? "" : ", ",
                            *c == 's' ? "scalar" : *c == 'x' ? "vertex_scalar"
                            : *c == 'v' ? "vector" : *c == 'f' ? "face_vector"
                            : *c == 't' ? "tensor" : "?");
  if (k < size)
    snprintf (buf + k, size - k, ")");
  return buf;
}

int ghost_ops_bind (scalar s, int kind, int face, int op, const scalar * uses,
                    int n, const char * sig)
{
  assert (s >= 0 && s < NSCALAR_MAX && kind >= 0 && kind <= GHOST_BOUNDARY);
  const char * want = ghost_op_signature (kind, op);
  if (sig && want && strcmp (want, sig) != 0) {
    char a [256], b [256];
    return ghost_ops_fail ("ghost_ops: scalar %s (%d): its %s %s depends on "
                           "%s, given %s\n", SNAME (s), (int) s,
                           ghost_what (kind, face), ghost_op_name (kind, op),
                           ghost_sig_text (want, a, sizeof a),
                           ghost_sig_text (sig, b, sizeof b));
  }
  return ghost_ops_set (s, kind, face, op, uses, n);
}

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
    if (uses [k] < 0 || uses [k] >= NSCALAR_MAX || !scalar_in_use (uses [k]))
      return ghost_ops_fail ("ghost_ops: scalar %s (%d): scalar %d, which its "
                             "%s depends on, is not in use\n", SNAME (s),
                             (int) s, (int) uses [k], ghost_what (kind, face));
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

int scalar_restore_ops (scalar s, ScalarOps o)
{
  int np = ghost_op_uses (GHOST_PROLONGATION, o.prolongation);
  int nr = ghost_op_uses (GHOST_RESTRICTION, o.restriction);
  if (np < 0 || nr < 0)
    return ghost_ops_fail ("scalar_restore_ops: scalar %s (%d): no such "
                           "operator\n", SNAME (s), (int) s);
  return ghost_ops_set (s, GHOST_PROLONGATION, GHOST_ALL_FACES, o.prolongation,
                        o.prolongation_uses, np)
      || ghost_ops_set (s, GHOST_RESTRICTION, GHOST_ALL_FACES, o.restriction,
                        o.restriction_uses, nr);
}

/* s about to be freed: none may depend on it. One that does stops the run;
.. or, with ghost_ops_abort 0, gets the default operators back. */
void ghost_ops_freeing (scalar s)
{
  assert (s >= 0 && s < NSCALAR_MAX);
  int d [GHOST_DEPS_MAX];
  foreach_scalar (t) {
    if (t == s)
      continue;
    for (int k = 0, n = ghost_uses_of (t, d); k < n; k++)
      if (d [k] == s) {
        ghost_ops_fail ("scalar_free: scalar %s (%d): the operators of "
                        "scalar %s (%d) depend on it; bind theirs again "
                        "first (scalar_restore_ops ())\n", SNAME (s), (int) s,
                        SNAME (t), (int) t);
        scalar_ops_reset (t);
        break;
      }
  }
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
  int loc = scalar_location (s);
  int face = loc >= CWRAP_FACE_X && loc <= CWRAP_FACE_Z, vertex = loc == CWRAP_VERTEX;
  if (!o.prolongation)
    o.prolongation = face ? prolong_face : vertex ? prolong_vertex : prolong_linear;
  if (!o.restriction)
    o.restriction = face ? restriction_face : vertex ? restriction_vertex
                                                     : restriction_mean;
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
  if (o.prolongation == prolong_none)
    return;                   /* filled by another's: no child visited */
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
  /*
  .. A face scalar along a: the children whose upper faces along a lie on
  .. the parent's (the upper half) take them too, at their face m, in the
  .. ghost layer: where no node of their level is above, those are theirs
  .. (block_ghosts () keeps them), the parent's face there, each fine face
  .. the coarse one it lies on.
  */
  int loc = scalar_location (s);
  if (loc == CWRAP_VERTEX && block_upper_walls (p))
  {
    /*
    .. A vertex scalar: the children on the parent's upper walls take the
    .. vertices they own there, at m along those axes, in their ghost
    .. layer, from the parent's, interpolated as prolong_vertex does.
    */
    int walls = block_upper_walls (p);
    const Real * coarse = block_data (index_node (p)->s, s);
    for (int c = 0; c < TREE_NCHILDREN; c++)
    {
      int cw = walls & c;                 /* the child's upper walls */
      if (!cw)
        continue;
      Real * fine = block_data (index_node (index_child (p, c))->s, s);
      int hi [3] = { 1, 1, 1 }, x [3] = { 0, 0, 0 };
      for (int d = 0; d < dimension; d++)
        hi [d] = m + 1;
      for (x [2] = 0; x [2] < hi [2]; x [2]++)
      for (x [1] = 0; x [1] < hi [1]; x [1]++)
      for (x [0] = 0; x [0] < hi [0]; x [0]++)
      {
        int at_m = 0;
        for (int d = 0; d < dimension; d++)
          at_m |= (x [d] == m) << d;
        if (!at_m || (at_m & ~cw))
          continue;
        /* the coarse vertices around: T, or T and T + 1 where between */
        int F [3] = { 0, 0, 0 };
        for (int d = 0; d < dimension; d++)
          F [d] = x [d] + ((c >> d) & 1) * m;
        RealSum v = 0;
        int n = 0;
        for (int q = 0; q < (1 << dimension); q++)
        {
          int X [3] = { 0, 0, 0 }, ok = 1;
          for (int d = 0; d < dimension; d++)
          {
            int up = (q >> d) & 1;
            if (up && !(F [d] & 1))
              ok = 0;                     /* at a coarse vertex: one only */
            X [d] = F [d] / 2 + up;
          }
          if (ok)
            v += coarse [block_at (sh, X)], n++;
        }
        fine [block_at (sh, x)] = (Real) (v / n);
      }
    }
  }
  if (loc >= CWRAP_FACE_X && loc <= CWRAP_FACE_Z)
  {
    int a = loc - CWRAP_FACE_X;
    const Real * coarse = block_data (index_node (p)->s, s);
    for (int c = 0; c < TREE_NCHILDREN; c++)
    {
      if (!((c >> a) & 1))
        continue;
      Real * fine = block_data (index_node (index_child (p, c))->s, s);
      int hi [3] = { 1, 1, 1 }, x [3] = { 0, 0, 0 };
      for (int d = 0; d < dimension; d++)
        hi [d] = d == a ? m + 1 : m;
      for (x [2] = a == 2 ? m : 0; x [2] < hi [2]; x [2]++)
      for (x [1] = a == 1 ? m : 0; x [1] < hi [1]; x [1]++)
      for (x [0] = a == 0 ? m : 0; x [0] < hi [0]; x [0]++)
      {
        int X [3] = { 0, 0, 0 };
        for (int d = 0; d < dimension; d++)
          X [d] = d == a ? m : (x [d] + ((c >> d) & 1) * m) / 2;
        fine [block_at (sh, x)] = coarse [block_at (sh, X)];
      }
    }
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
  if (o.prolongation == prolong_none) {
    fprintf (stderr, "block_prolong_into: scalar %s (%d): its prolongation "
             "is prolong_none, so it has no wavelet of its own: adapt on the "
             "scalar whose operators fill it\n", SNAME (s), (int) s);
    abort ();
  }
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
  if (o.restriction != restriction_none)
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
static void block_ghosts_cells (Index p, int s);
static void block_walls (Index p, int s);

/* the faces along a at X of block b, out to buffer (put 0) or back (put 1):
.. every transverse interior place */
static void face_row (BlockShape sh, Real * b, int a, int X, Real * buf, int put)
{
  int hi [3] = { 1, 1, 1 }, x [3] = { 0, 0, 0 }, k = 0;
  for (int d = 0; d < dimension; d++)
    hi [d] = d == a ? X + 1 : sh.m;
  for (x [2] = a == 2 ? X : 0; x [2] < hi [2]; x [2]++)
  for (x [1] = a == 1 ? X : 0; x [1] < hi [1]; x [1]++)
  for (x [0] = a == 0 ? X : 0; x [0] < hi [0]; x [0]++)
  {
    long c = block_at (sh, x);
    if (put)
      b [c] = buf [k++];
    else
      buf [k++] = b [c];
  }
}

/* the vertices of a block at m along some axes, all of them upper walls
.. (walls, a mask): those it owns there; out to buf (put 0) or back (1) */
static int vertex_walls_row (BlockShape sh, Real * b, int walls, Real * buf,
                             int put)
{
  int hi [3] = { 1, 1, 1 }, x [3] = { 0, 0, 0 }, k = 0;
  for (int d = 0; d < dimension; d++)
    hi [d] = sh.m + 1;
  for (x [2] = 0; x [2] < hi [2]; x [2]++)
  for (x [1] = 0; x [1] < hi [1]; x [1]++)
  for (x [0] = 0; x [0] < hi [0]; x [0]++)
  {
    int at_m = 0;
    for (int d = 0; d < dimension; d++)
      at_m |= (x [d] == sh.m) << d;
    if (!at_m || (at_m & ~walls))
      continue;
    long c = block_at (sh, x);
    if (put)
      b [c] = buf [k];
    else
      buf [k] = b [c];
    k++;
  }
  return k;
}

static Index index_step (Index p, int a, int dir)
{
  if (a == 0)
    p.i += dir;
#if dimension > 1
  else if (a == 1)
    p.j += dir;
#endif
#if dimension > 2
  else
    p.k += dir;
#endif
  return p;
}

/*
.. The ghosts of p for s, in two phases. Cell scalars have all of theirs in
.. the first. Faces and vertices have those of a cell scalar, but where
.. blocks meet:
..
.. A face scalar (a face_vector's component along a), first phase: its face
.. m, the upper one, kept as the block made it where no node of its level
.. is above (a wall, a coarser leaf: the fine side's faces); else that
.. node's face 0, copied as any ghost: at a finer region a parent's, the
.. mean of its children's (restriction_face). Second phase: a leaf's face 0
.. where the node below is a parent (a finer region): the mean of the fine
.. faces on it, its children's faces m.
..
.. A vertex scalar, first phase: the vertices it owns on the upper walls
.. kept. Second phase: those on an upper wall at m along axes that are not
.. walls too, from their owner, the node of this level across those axes
.. (the walls' conditions made them here as a cell's).
..
.. The second phase reads only what the first has finished: ghosts_update
.. runs it once every node of the level has had the first, at once.
*/
static int block_ghosts_located (int s)
{
  int loc = scalar_location (s);
  return loc == CWRAP_VERTEX || (loc >= CWRAP_FACE_X && loc <= CWRAP_FACE_Z);
}

static void block_ghosts_first (Index p, int s)
{
  int loc = scalar_location (s);
  BlockShape sh = tree_block_shape ();
  Real * b = block_data (index_node (p)->s, s);
  if (loc == CWRAP_VERTEX)
  {
    int walls = block_upper_walls (p);
    if (!walls)
    {
      block_ghosts_cells (p, s);
      return;
    }
    long n = 1;
    for (int d = 0; d < dimension; d++)
      n *= sh.m + 1;
    Real * buf = malloc ((size_t) n * sizeof (Real));
    if (!buf)
      abort ();
    vertex_walls_row (sh, b, walls, buf, 0);
    block_ghosts_cells (p, s);
    vertex_walls_row (sh, b, walls, buf, 1);
    free (buf);
    return;
  }
  if (loc < CWRAP_FACE_X || loc > CWRAP_FACE_Z)
  {
    block_ghosts_cells (p, s);
    return;
  }
  int a = loc - CWRAP_FACE_X, row = 1;
  assert (a < dimension);
  for (int d = 0; d < dimension; d++)
    row *= d == a ? 1 : sh.m;
  Real saved [row];
  int keep = !block_neighbour_at_level (p, a, 1);
  if (keep)
    face_row (sh, b, a, sh.m, saved, 0);
  block_ghosts_cells (p, s);
  if (keep)
    face_row (sh, b, a, sh.m, saved, 1);
}

static void block_ghosts_second (Index p, int s)
{
  int loc = scalar_location (s);
  BlockShape sh = tree_block_shape ();
  int m = sh.m;
  Real * b = block_data (index_node (p)->s, s);
  if (loc == CWRAP_VERTEX)
  {
    int walls = block_upper_walls (p), hi [3] = { 1, 1, 1 }, x [3] = { 0, 0, 0 };
    if (!walls)
      return;
    for (int d = 0; d < dimension; d++)
      hi [d] = m + 1;
    for (x [2] = 0; x [2] < hi [2]; x [2]++)
    for (x [1] = 0; x [1] < hi [1]; x [1]++)
    for (x [0] = 0; x [0] < hi [0]; x [0]++)
    {
      int at_m = 0;
      for (int d = 0; d < dimension; d++)
        at_m |= (x [d] == m) << d;
      int across = at_m & ~walls;
      if (!(at_m & walls) || !across)
        continue;
      Index q = p;
      int y [3] = { x [0], x [1], x [2] };
      for (int d = 0; d < dimension; d++)
        if ((across >> d) & 1)
          q = index_step (q, d, 1), y [d] = 0;
      if (index_is_allocated (q))
        b [block_at (sh, x)] = block_data (index_node (q)->s, s) [block_at (sh, y)];
    }
    return;
  }
  if (loc < CWRAP_FACE_X || loc > CWRAP_FACE_Z)
    return;
  int a = loc - CWRAP_FACE_X;
  if (!(index_node (p)->flags & BAMR_NODE_LEAF) || !block_neighbour_at_level (p, a, -1))
    return;
  Index lo = index_step (p, a, -1);
  if (!tree_has_children (index_node (lo)))
    return;
  /* each face of row 0: the mean of the fine faces on it, at m in lo's
   * upper children along a (restriction_face's, there) */
  int hi [3] = { 1, 1, 1 }, T [3] = { 0, 0, 0 };
  for (int d = 0; d < dimension; d++)
    hi [d] = d == a ? 1 : m;
  for (T [2] = 0; T [2] < hi [2]; T [2]++)
  for (T [1] = 0; T [1] < hi [1]; T [1]++)
  for (T [0] = 0; T [0] < hi [0]; T [0]++)
  {
    RealSum sum = 0;
    int n = 0;
    for (int q = 0; q < (1 << (dimension - 1)); q++)
    {
      int local [3] = { 0, 0, 0 }, child = 1 << a, bit = 0;
      for (int d = 0; d < dimension; d++)
      {
        if (d == a)
        {
          local [d] = m;
          continue;
        }
        int F = 2 * T [d] + ((q >> bit++) & 1), up = F >= m;
        child |= up << d;
        local [d] = F - up * m;
      }
      sum += block_data (index_node (index_child (lo, child))->s, s)
             [block_at (sh, local)];
      n++;
    }
    int X [3] = { T [0], T [1], T [2] };
    X [a] = 0;
    b [block_at (sh, X)] = (Real) (sum / n);
  }
}

void block_ghosts (Index p, int s)
{
  block_ghosts_first (p, s);
  if (block_ghosts_located (s))
    block_ghosts_second (p, s);
}

static void block_ghosts_cells (Index p, int s)
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
  int prolonged = 0;                  /* a region from a coarser leaf */

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

    /* one level coarser, 2:1 balance: the leaf covering b; nothing to do
     * for prolong_none (another's prolongation fills the region) */
    if (o.prolongation == prolong_none)
      continue;
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
    prolonged = 1;
  }

  block_walls (p, s);
  /* a scalar s's prolongation depends on, whose own does nothing: s's may
   * have written its ghosts at level jumps, after its walls were made from
   * them; made again (prolong_none, block-ghosts.h) */
  for (int k = 0, n = prolonged ? ghost_op_uses (GHOST_PROLONGATION, o.prolongation)
                                 : 0; k < n; k++)
    if (scalar_get_ops (o.prolongation_uses [k]).prolongation == prolong_none)
      block_walls (p, o.prolongation_uses [k]);
}

/* the walls of block p, one axis after the other */
static void block_walls (Index p, int s)
{
  BlockShape sh = tree_block_shape ();
  int nb = 1 << p.level;
  Real * own = block_data (index_node (p)->s, s);
  int P [3] = { p.i, 0, 0 };
  #if dimension > 1
  P [1] = p.j;
  #endif
  #if dimension > 2
  P [2] = p.k;
  #endif
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
  if (o.restriction == restriction_none)
    return;                   /* no parent visited */
  /* the parents, finest level first: each after its children */
  const TreeCaches * tc = tree_caches ();
  for (int l = tc->depth - 1; l >= 0; l--)
    for (int q = tc->parent_first [l]; q < tc->parent_first [l + 1]; q++)
      @restriction (o.restriction, o.restriction_uses)
        (tc->parents_by_level [q].p, s);
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
  if (scalar_get_ops (s).restriction != restriction_none)
    for (int k = L - 1; k >= tree_get_depth_min (); k--)
      for (int q = tc->parent_first [k]; q < tc->parent_first [k + 1]; q++)
        block_restrict (tc->parents_by_level [q].p, s);
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

  /* 1. restriction, finest first; within a level, at once: the scalars
   * whose restriction does something (restriction_none: no parent visited),
   * over the parents only (the tree's list of them). Without the parents'
   * own ghosts, down to the coarsest leaf's level only: a parent is read
   * by a leaf of its level beside it, and below that level there is none
   * (on a uniform tree, no level at all) */
  int r [NSCALAR_MAX], nr = 0;
  for (int q = 0; q < n; q++)
    if (o [q].restriction != restriction_none)
      r [nr++] = q;
  const TreeItem * par = tc->parents_by_level;
  const int * pfirst = tc->parent_first;
  int lowest = parents ? 0 : tree_get_depth_min ();
  for (int l = tc->depth - 1; l >= lowest && nr; l--)
  {
    BLOCK_PRAGMA (omp parallel for schedule (dynamic))
    for (int b = pfirst [l]; b < pfirst [l + 1]; b++)
      for (int k = 0; k < nr; k++)
        @restriction (o [r [k]].restriction, o [r [k]].restriction_uses)
          (par [b].p, s [r [k]]);
  }

  /* 2. the ghosts, coarsest first; within a level, at once */
  for (int l = 0; l <= tc->depth; l++)
  {
    BLOCK_PRAGMA (omp parallel for schedule (dynamic))
    for (int b = first [l]; b < first [l + 1]; b++)
      if (parents || (at [b].node->flags & BAMR_NODE_LEAF))
        for (int q = 0; q < n; q++)
          block_ghosts_first (at [b].p, s [q]);
    /* faces and vertices: their second phase, once the level's first is
     * done (it reads what the first leaves) */
    int located = 0;
    for (int q = 0; q < n; q++)
      located |= block_ghosts_located (s [q]);
    if (located)
    {
      BLOCK_PRAGMA (omp parallel for schedule (dynamic))
      for (int b = first [l]; b < first [l + 1]; b++)
        if (parents || (at [b].node->flags & BAMR_NODE_LEAF))
          for (int q = 0; q < n; q++)
            if (block_ghosts_located (s [q]))
              block_ghosts_second (at [b].p, s [q]);
    }
  }
  for (int q = 0; q < n; q++)
    ghosts_dirty [s [q] / 64] &= ~(1ull << (s [q] % 64));
  free (o);
}

void ghosts_update_dirty_located (int n, const int * s)
{
  int list [NSCALAR_MAX], k = 0;
  for (int q = 0; q < n && k < NSCALAR_MAX; q++)
    if (s [q] >= 0 && s [q] < NSCALAR_MAX && scalar_in_use (s [q])
        && ghosts_is_dirty (s [q]) && scalar_location (s [q]) != CWRAP_CENTRE)
      list [k++] = s [q];
  if (k)
    ghosts_update_dirty (k, list);
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
