/*
.. The ghost operators libb-amr gives, as cwrap definitions (block-ghosts.h):
.. included once, by the one file of a program that dispatches to its
.. operators, before any of its own, so that these keep the ids 1, 2, ...
.. in every program. The defaults come first: prolong_linear,
.. restriction_mean, bc_neumann.
*/
#ifndef _BAMR_BLOCK_OPS_H_
#define _BAMR_BLOCK_OPS_H_

#include <assert.h>
#include <stdlib.h>
#include <string.h>

#include "block-ghosts.h"

/* ---------------------------------------------------------------------- */
/* prolongations, restriction, boundary conditions                        */
/* ---------------------------------------------------------------------- */

@prolongation prolong_linear (const Prolongation * pr)
{
  scalar s = pr->s;
  foreach_prolongation (pr) {
    RealSum v = coarse (s);
    v += child.x * (coarse (s, 1) - coarse (s, -1)) / 8.;
#if dimension > 1
    v += child.y * (coarse (s, 0, 1) - coarse (s, 0, -1)) / 8.;
#endif
#if dimension > 2
    v += child.z * (coarse (s, 0, 0, 1) - coarse (s, 0, 0, -1)) / 8.;
#endif
    s[] = (Real) v;
  }
}

static inline RealSum minmod (RealSum a, RealSum b)
{
  return a * b <= 0 ? 0 : a > 0 ? (a < b ? a : b) : (a > b ? a : b);
}

@prolongation prolong_minmod (const Prolongation * pr)
{
  scalar s = pr->s;
  foreach_prolongation (pr) {
    RealSum v = coarse (s);
    v += child.x * minmod ((RealSum) coarse (s, 1) - coarse (s),
                           (RealSum) coarse (s) - coarse (s, -1)) / 4.;
#if dimension > 1
    v += child.y * minmod ((RealSum) coarse (s, 0, 1) - coarse (s),
                           (RealSum) coarse (s) - coarse (s, 0, -1)) / 4.;
#endif
#if dimension > 2
    v += child.z * minmod ((RealSum) coarse (s, 0, 0, 1) - coarse (s),
                           (RealSum) coarse (s) - coarse (s, 0, 0, -1)) / 4.;
#endif
    s[] = (Real) v;
  }
}

@prolongation prolong_injection (const Prolongation * pr)
{
  scalar s = pr->s;
  foreach_prolongation (pr)
    s[] = coarse (s);
}

@restriction restriction_mean (Index p, scalar s)
{
  assert (index_is_allocated (p) && tree_has_children (index_node (p)));
  foreach_restriction (p) {
    RealSum sum = 0;
    foreach_child ()
      sum += s[];
    s[] = (Real) (sum / TREE_NCHILDREN);
  }
}

@boundary bc_neumann (const Boundary * bc)
{
  boundary_cells (bc)
    bc->b [ghost] = bc->b [inner];
}

@boundary bc_odd (const Boundary * bc)
{
  boundary_cells (bc)
    bc->b [ghost] = - bc->b [inner];
}

#endif
