#ifndef _VOF_VOF_H_
#define _VOF_VOF_H_

  /*
  .. Volume of fluid on libb-amr's tree: c, the fraction of each cell that
  .. is fluid, carried by a face velocity uf.
  ..
  ..   vof_advection (c, uf, dt)
  ..
  .. One step of the conservative split scheme of Weymouth and Yue (2010),
  .. as libvof/vof/vof_advect.cpp: a sweep along each axis, the order
  .. rotated at every call (x y z, then y z x, then z x y; in 2D x y, y x).
  .. A sweep along x:
  ..
  ..   - the plane of every leaf, n . x <= alpha: n Youngs' normal, alpha
  ..     its position holding c; at a level jump the coarse cell's plane,
  ..     with c (prolong_plane);
  ..   - the flux through each face, made once per face (foreach_face_x):
  ..     the volume that crosses it in dt, cut from the upwind cell by its
  ..     interface, a plane (PLIC): the upwind cell's;
  ..   - each cell:
  ..       c += (flux.x[] - flux.x[1] + cc (uf.x[1] - uf.x[]) dt) / Delta
  ..     with cc = 1 where c > 1/2 at the start of the step, 0 elsewhere,
  ..     the same in every sweep: the sweeps' cc terms add up to
  ..     cc div (uf) dt, zero for a solenoidal uf.
  ..
  .. The fluxes are a face vector, so at a level jump a coarse face's flux
  .. is the mean of the fine ones on it (refluxing, libb-amr): the integral
  .. of c over the leaves changes only by what crosses the walls. A fine
  .. face's upwind cell may be a coarse leaf's ghost, c's prolongation.
  ..
  .. Conditions:
  ..   - |uf| dt / Delta <= 1/2 on every face, for 0 <= c <= 1 without the
  ..     clamp (below); the scheme takes a slab of the upwind cell at most
  ..     a cell wide.
  ..   - uf divergence free, face by face (on each leaf the sum of
  ..     uf.x[1] - uf.x[] over the axes zero), the coarse faces at level
  ..     jumps read as the mean of the fine (as libb-amr makes them), for c
  ..     to be conserved.
  ..   - Outside a call, c's prolongation, for refinement, should make no
  ..     new extrema: prolong_fraction (vof-ghosts.h) or prolong_minmod.
  ..     During a call c's operators are prolong_plane and
  ..     restriction_plane (vof-ghosts.h), which fill c, n and alpha
  ..     together; c's own are given back on return. vof_advection () warns
  ..     once of the linear default.
  ..   - The tree must not change during a call.
  ..
  .. c is clamped to [0, 1] after each sweep, which removes round-off only
  .. under the conditions above. A cell with c <= VOF_EPS gives nothing, one
  .. with c >= 1 - VOF_EPS gives its slab whole (no reconstruction).
  */

  #include "tree.h"

  #ifndef VOF_EPS
  # define VOF_EPS (sizeof (Real) < 8 ? 1e-6 : 1e-10)
  #endif

  void vof_advection (scalar c, face_vector uf, double dt);

  /*
  .. The geometry of a plane in the unit cell [0, 1]^dimension, the fluid
  .. on the side n . x <= alpha, n of any length and signs (not zero); the
  .. members beyond the dimension are not read:
  ..
  ..   vof_plane_alpha (n, f)        alpha such that the fluid fills the
  ..                                 fraction f of the cell
  ..   vof_plane_volume (n, alpha, lo, hi)   the fraction of the box
  ..                                 [lo, hi] (within the cell) that is
  ..                                 fluid
  ..
  .. A normal with its components permuted gives the same values, bit for
  .. bit (the components are summed in sorted order).
  */
  double vof_plane_alpha (coord n, double f);
  double vof_plane_volume (coord n, double alpha, coord lo, coord hi);

#endif
