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
    foreach_dimension ()
      v += child.x * (coarse (s, 1) - coarse (s, -1)) / 8.;
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
    foreach_dimension ()
      v += child.x * minmod ((RealSum) coarse (s, 1) - coarse (s),
                             (RealSum) coarse (s) - coarse (s, -1)) / 4.;
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

/*
.. Faces: a face_vector's component along axis a holds, at a cell, its lower
.. face along a (libc-pp's scalar.h). Their defaults, by location
.. (scalar_location ()):
..
.. restriction_face: a parent's faces 0 .. m along a (its upper face m too,
.. in its ghost layer), each the mean of the 2^(d-1) fine faces on it: the
.. flux through a coarse face is the fluxes through its fine faces.
..
.. prolong_face: a fine cell's lower face along a lies on the coarse cell's
.. lower face (the lower half along a), which it takes; or across the
.. coarse cell's middle, the mean of its two faces along a. The same along
.. the other axes: a flux per unit area, as the coarse face's.
*/
@restriction restriction_face (Index p, scalar s)
{
  int a = scalar_location (s) - CWRAP_FACE_X;
  BlockShape sh = tree_block_shape ();
  int m = sh.m, hi [3] = { 1, 1, 1 };
  assert (a >= 0 && a < dimension && tree_has_children (index_node (p)));
  Real * b = block_data (index_node (p)->s, s);
  for (int d = 0; d < dimension; d++)
    hi [d] = d == a ? m + 1 : m;
  int X [3];
  for (X [2] = 0; X [2] < hi [2]; X [2]++)
  for (X [1] = 0; X [1] < hi [1]; X [1]++)
  for (X [0] = 0; X [0] < hi [0]; X [0]++)
  {
    /* the fine faces, in the parent's fine indices: 2 X along a, 2 X and
     * 2 X + 1 along the others; then their child and place in it */
    RealSum sum = 0;
    int n = 0;
    for (int q = 0; q < (1 << (dimension - 1)); q++)
    {
      int F [3] = { 0, 0, 0 }, local [3] = { 0, 0, 0 }, child = 0, bit = 0;
      for (int d = 0; d < dimension; d++)
      {
        F [d] = 2 * X [d];
        if (d != a)
          F [d] += (q >> bit++) & 1;
        int up = F [d] >= m;
        child |= up << d;
        local [d] = F [d] - up * m;
      }
      sum += block_data (index_node (index_child (p, child))->s, s)
             [block_at (sh, local)];
      n++;
    }
    b [block_at (sh, X)] = (Real) (sum / n);
  }
}

enum { face_location_x = CWRAP_FACE_X, face_location_y = CWRAP_FACE_Y,
       face_location_z = CWRAP_FACE_Z };

@prolongation prolong_face (const Prolongation * pr)
{
  scalar s = pr->s;
  int location = scalar_location (s);
  foreach_prolongation (pr) {
    foreach_dimension ()
      if (location == face_location_x)
        s[] = child.x < 0 ? coarse (s)
              : (Real) (((RealSum) coarse (s) + coarse (s, 1)) / 2);
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
