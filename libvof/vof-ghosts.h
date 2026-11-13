#ifndef _VOF_VOF_GHOSTS_H_
#define _VOF_VOF_GHOSTS_H_

  /*
  .. The ghost operators of libvof, for a volume fraction c:
  ..
  ..   prolong_fraction   each coarse cell cut by its interface, a plane
  ..                      (Youngs' normal from the coarse cells around,
  ..                      its position by the cell's fraction), and each
  ..                      child given the fraction of it below the plane
  ..
  ..   ghost_ops (c, prolongation (prolong_fraction));
  ..
  .. The children's mean is the parent's fraction (to round-off), each in
  .. [0, 1], and the interface stays as sharp as the parent's: for the
  .. ghosts at level jumps (the fine faces' upwind cells, vof_advection ())
  .. and for refinement. The restriction is the mean, the default. It reads
  .. no other scalar, so it can be bound to any c, local ones too.
  ..
  .. A program that binds it builds vof-ghosts.c, the file that dispatches
  .. to libb-amr's operators and these, in place of libb-amr's
  .. block-ghosts.c (block-ghosts.h).
  */

  #include "block-ghosts.h"

  @prolongation prolong_fraction;

  /*
  .. c with its plane, n (Youngs' normal, minus the gradient) and alpha
  .. (n . x <= alpha in the unit cell), as vof_advection () binds them for
  .. its call: c's operators fill all three, n's and alpha's do nothing
  .. (prolong_none, restriction_none), so one ghost update of c makes them
  .. all:
  ..   prolong_plane      the fine cells cut by the coarse cell's plane, n
  ..                      the coarse cell's, alpha moved to the fine cell
  ..   restriction_plane  c the mean, n the mean of the children's, alpha
  ..                      the plane of n holding the parent's c
  .. Declared with what they depend on, so that cwrap checks a binding:
  ..   ghost_ops (c, prolongation (prolong_plane, n, alpha),
  ..                 restriction (restriction_plane, n, alpha));
  */
  @prolongation prolong_plane (const Prolongation * pr) (vector n, scalar alpha);
  @restriction  restriction_plane (Index p, scalar s) (vector n, scalar alpha);

#endif
