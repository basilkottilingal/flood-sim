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
  const Real * u = pr->coarse;
  long st [3];
  for (int a = 0; a < 3; a++)
    st [a] = block_stride (pr->shape, a);
  prolong_cells (pr)
  {
    RealSum v = u [cc];
    for (int a = 0; a < dimension; a++)
      v += sgn [a] * (u [cc + st [a]] - u [cc - st [a]]) / 8.;
    pr->fine [f] = (Real) v;
  }
}

static inline RealSum minmod (RealSum a, RealSum b)
{
  return a * b <= 0 ? 0 : a > 0 ? (a < b ? a : b) : (a > b ? a : b);
}

@prolongation prolong_minmod (const Prolongation * pr)
{
  const Real * u = pr->coarse;
  long st [3];
  for (int a = 0; a < 3; a++)
    st [a] = block_stride (pr->shape, a);
  prolong_cells (pr)
  {
    RealSum v = u [cc];
    for (int a = 0; a < dimension; a++)
      v += sgn [a] * minmod ((RealSum) u [cc + st [a]] - u [cc],
                             (RealSum) u [cc] - u [cc - st [a]]) / 4.;
    pr->fine [f] = (Real) v;
  }
}

@prolongation prolong_injection (const Prolongation * pr)
{
  prolong_cells (pr)
    pr->fine [f] = pr->coarse [cc];
}

@restriction restriction_mean (Index p, scalar s)
{
  BlockShape sh = tree_block_shape ();
  int m = sh.m, h = m / 2;
  assert (index_is_allocated (p) && tree_has_children (index_node (p)));
  Real * parent = block_data (index_node (p)->s, s);
  const Real * child [TREE_NCHILDREN];
  for (int c = 0; c < TREE_NCHILDREN; c++)
    child [c] = block_data (index_node (index_child (p, c))->s, s);
  /* the 2^dimension fine cells of a coarse one, from the first */
  long sub [TREE_NCHILDREN];
  for (int k = 0; k < TREE_NCHILDREN; k++)
  {
    sub [k] = 0;
    for (int a = 0; a < dimension; a++)
      sub [k] += ((k >> a) & 1) * block_stride (sh, a);
  }
  /* child c fills the parent's cells c's half along each axis: by rows */
  int hy = dimension > 1 ? h : 1, hz = dimension > 2 ? h : 1;
  for (int c = 0; c < TREE_NCHILDREN; c++)
    for (int Z = 0; Z < hz; Z++)
    for (int Y = 0; Y < hy; Y++)
    {
      int px [3] = { (c & 1) * h, ((c >> 1) & 1) * h + Y, ((c >> 2) & 1) * h + Z };
      int fx [3] = { 0, 2 * Y, 2 * Z };
      Real * prow = parent + block_at (sh, px);
      const Real * frow = child [c] + block_at (sh, fx);
      for (int X = 0; X < h; X++)
      {
        RealSum sum = 0;
        for (int k = 0; k < TREE_NCHILDREN; k++)
          sum += frow [2 * X + sub [k]];
        prow [X] = (Real) (sum / TREE_NCHILDREN);
      }
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
