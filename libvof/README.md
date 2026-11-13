# libvof

Volume of fluid on libb-amr's block tree, written with cwrap's stencil
syntax. For now, the advection of a volume fraction:

```c
#include "vof.h"
#include "vof-ghosts.h"

scalar c;                 /* the fraction of each cell that is fluid */
face_vector uf;           /* uf.x at a cell's lower face along x     */
...
ghost_ops (c, prolongation (prolong_fraction));   /* see below */
...
vof_advection (c, uf, dt);
```

A program that binds `prolong_fraction` builds `vof-ghosts.c` in place of
libb-amr's `block-ghosts.c`: it is the file that dispatches to libb-amr's
operators and libvof's (libb-amr's `block-ghosts.h`). Only one such file
goes into a program, so a program with operators of its own, or another
library's, writes one that includes them all.

`vof/` and `AMR-2D-vofbranch/` hold the earlier C++ code that this one
follows (`vof/vof_advect.cpp`, `vof/vof_utils.cpp`).

## The scheme

The split, conservative scheme of Weymouth and Yue (2010), as
`VOFSolver::weymouth_advection ()`: one sweep along each axis, the order
turned at every call, x y z, then y z x, then z x y (in 2D, x y then y x),
as `toggle` turned it. A sweep along x:

1. **The flux through each face** (`foreach_face_x ()`, each face once):
   the volume that crosses it in dt, `u dt` times the fraction of fluid in
   the slab of the upwind cell that crosses, `|u| dt / Delta` of it wide.
   The slab's fluid is cut by the cell's interface, a plane (PLIC),
   `n . x <= alpha`: its normal by Youngs' stencil on the 3^d cells around
   (`calc_slope ()`), its position so that the cell holds its fraction
   (`calc_alpha ()`), both made on every leaf before the sweep
   (`reconstruct ()`, all the branching of the reconstruction, once per
   cell). An empty cell gives nothing; a full one its slab whole.
2. **Each cell:**

       c += (flux.x[] - flux.x[1] + cc (uf.x[1] - uf.x[]) dt) / Delta

   with `cc` 1 where c > 1/2 at the start of the step, 0 elsewhere, the
   same in every sweep. Its terms add up over the sweeps to `cc div (uf)
   dt`, zero for a solenoidal `uf`, and they keep c within [0, 1] for a
   Courant number up to 1/2 along each axis. c is clamped to [0, 1] after
   each sweep, which then removes only round-off.

The plane geometry (`vof_plane_alpha ()`, `vof_plane_volume ()`, on
`coord`s) is that of
Scardovelli and Zaleski (2000), as `vof_utils.cpp` (Ahmed Basil Kottilingal
- old VOF code. link missing), written to divide by the
smallest component of the normal only where the quotient is bounded: so a
normal close to an axis needs no threshold. Where the volume is a cubic in
alpha (3D), alpha is found by Newton's method, bisection where a step
would leave the bracket, in place of the trigonometric root.

**On the tree:** the fluxes are a face vector, so at a level jump the
coarse face's flux is the mean of the fine faces' on it (refluxing,
libb-amr): the integral of c over the leaves changes only by what crosses
the walls. `uf` is refluxed the same way when it is read.

## Conditions

- `|uf| dt / Delta <= 1/2` on every face.
- `uf` divergence free on each leaf, the coarse faces at level jumps as
  the mean of the fine.
- **Outside a call, bind `prolong_fraction` to c** (or `prolong_minmod`),
  for refinement and adaptation: the linear default makes children beyond
  [0, 1] at an interface. During a call c's ghosts at level jumps are
  `prolong_plane`'s, whatever c's own operators (below). With linear ghosts
  there, the coarse cell gave more than it held, the clamp made up the
  difference, and c was lost: 12% through the 2D vortex below.
  `vof_advection ()` warns once if c has the default.
- The tree must not change during a call.

## prolong_fraction

`vof-ghosts.h`: each coarse cell is cut by its interface, Youngs' normal from
the coarse cells around (`coarse (c, i, j)`, one cell away) and its position
by its fraction, and each child gets the fraction of its part of the cell
below the plane, as Basilisk's `fraction_refine`. The children's mean is the
parent's fraction, each is within [0, 1], and the interface stays as sharp as
the parent's, at level jumps and on refinement. The restriction is the mean.

It reads no other scalar, so it can stay bound to c: for refinement and
adaptation outside `vof_advection ()`.

## c, n and alpha together

During a call, `vof_advection ()` holds the plane of every leaf, `n` (a
vector) and `alpha`, and binds c's operators so that they fill all three:

    @prolongation prolong_plane (const Prolongation * pr) (vector n, scalar alpha)
    @restriction  restriction_plane (Index p, scalar s) (vector n, scalar alpha)

    ghost_ops (c, prolongation (prolong_plane, n, alpha),
                  restriction (restriction_plane, n, alpha));
    n.x, n.y (, n.z), alpha: prolong_none, restriction_none

One definition and one binding in every dimension: what they depend on is a
vector and a scalar, declared so in `vof-ghosts.h`, and cwrap checks the
binding against that declaration (libc-pp's README, "Coupled operators").

- **`prolong_plane`:** a fine cell takes the coarse cell's plane: `n` the
  same, `alpha` moved to its own unit cell (`2 alpha - n . o`, `o` its
  corner in the coarse cell), and c the fraction of it below the plane.
- **`restriction_plane`:** a parent's c is the mean, `n` the mean of its
  children's, `alpha` the plane of that `n` holding the parent's c.
- **One update makes the three.** Writing `n` and `alpha` (`reconstruct
  ()`) marks c stale, as c depends on them. The sweep reads all three at an
  offset, so they are updated at once: every restriction first (c's fills
  the parents' `n` and `alpha`), then the ghosts level by level, `n` and
  `alpha` copied where a node of their level is beside, then c's
  prolongation filling the three at level jumps. libb-amr makes `n`'s and
  `alpha`'s walls again after it (`prolong_none`, libb-amr's README).
- **The upwind ghost at a level jump** is thus a cut of the coarse cell's
  own plane: its c, `n` and `alpha` agree, and the slab is cut by the
  plane the coarse cell has.
- **The binding is for the call:** `n` and `alpha` are local scalars, which
  libb-amr lets an operator depend on while they live; c's own operators
  (`scalar_get_ops ()` on entry) are given back by `scalar_restore_ops ()`
  before they go.

## Tests

    sh tests/run.sh           # 1, 2 and 3D, double and float, sanitizers
    DIMS=2 sh tests/run.sh    # one dimension

`tests/test_vof.c`:
- **the geometry:** alpha and the volume inverse to each other (to 1e-15,
  random normals, some with components of 1e-9 or 0, fractions of 1e-12
  and 1 - 1e-12); a cell cut in two adds up; against sampled points;
- **a slab carried along x** across a level jump: exact (1e-14), c kept;
- **`prolong_fraction`:** a plane's fractions, refined below the middle: c
  kept to round-off, within [0, 1], exact for a plane along an axis (minmod:
  0.13 off), at most 0.018 off for an oblique one in 2D, 0.23 in 3D
  (minmod: 0.42, 0.54);
- **a sphere (circle) through a vortex, reversed** (the flux of each face
  the circulation of a vector potential, so `uf` is divergence free face
  by face and zero at the walls): c conserved to 1e-15 at every step, on
  uniform trees and across level jumps, within [0, 1], and the shape back.
  2D, 64^2: L1 error 9.4e-4; 128^2: 3.8e-4; 64^2 with half refined to
  128^2: 5.2e-4 (6.1e-4 with minmod ghosts). 3D, 32^3: 5.4e-4; with half refined
  to 64^3: 4.7e-4 (5.7e-4 with minmod ghosts);
- **the order turned:** a problem turned as the axes are (x to y to z, in
  2D x and y swapped), advanced by the call after the unturned one, is the
  unturned one's result turned, bit for bit; by a call with the unturned
  one's order, not.

The whole takes about 3.5 minutes, 3D most of it.

## How it is written

The sweep is written once, for x, and cwrap writes it for y and z
(`foreach_dimension ()` before the function, libc-pp's README):

```c
foreach_dimension ()
static void sweep_x (scalar c, scalar cc, vector n, scalar alpha,
                     face_vector uf, face_vector flux, double dt)
{
  foreach_face_x () {               /* sweep_y: foreach_face_y (), c[0,-1] */
    double u = uf.x[], f = u > 0 ? c[-1] : c[];
    ...
  }
  foreach ()
    c[] += (flux.x[] - flux.x[1] + cc[] * (uf.x[1] - uf.x[]) * dt) / Delta;
}
...
foreach_dimension ()
  if (axis == axis_x)
    sweep_x (c, cc, n, alpha, uf, flux, dt);
```

Youngs' normal is computed on every leaf before each sweep, by one stencil
written for x and turned (`foreach_dimension ()` in the loop), with the
plane's `alpha`; the face loop cuts the upwind cell's plane. The plane
geometry takes `coord`s, and sums the components in sorted order, so that
the sweep along y of a problem turned is the sweep along x, bit for bit (the
order test).

## SIMD and GPU

What the code does for vector lanes and GPU threads, and what it measures.

- **The geometry has no branch on the data** (`vof_plane_alpha ()`,
  `vof_plane_volume ()`): the components are sorted by a network of `fmin`
  and `fmax`, every case of the volume is computed and the right one
  selected, the cubic of 3D solved by six Newton steps from the upper end of
  its bracket (the volume is convex there, so the steps come down to the
  root: six reach round-off over two million random planes), and nothing is
  an array indexed as it runs (on a GPU, such arrays leave the registers).
  The divisions are by reciprocals made once per plane.
- **The branches left are coherent:** "is the upwind cell, or the cell, an
  interface one" is the same for whole regions (the bulk on either side),
  so SIMD lanes and GPU warps mostly agree on it; the work under it is the
  geometry, for the few interface cells. Computing it for every cell, to have
  no branch at all, would cost the geometry on all of them.
- **`reconstruct ()` is two loops:** Youngs' normal, a stencil with no
  branch, on every cell (vectorised along rows), then alpha on the interface
  cells.

Measured (`sh tests/bench.sh [flags]`: a sphere in the vortex, uniform
trees, 128^2 and 64^3 cells, 0.9% and 0.7% of them interface cells; one
core, gcc -O3, a loaded laptop, so rough):

| | 2D ns/cell/step | 3D ns/cell/step |
|---|--:|--:|
| before (branching geometry, one reconstruction loop) | 72 | 223-227 |
| now | 62 | 203-226 |
| now, `-fopenmp-simd -DBLOCK_OPENMP_SIMD -march=native` | 43 | 190 |

Per call, as scalar code, the geometry without branches is about 10% slower
than with them (alpha 194 ns against 174, a slab's volume 56 against 50):
it computes every case. It is called on interface cells only, so that
does not show; it is what lets a vector or warp run it in step.

Where the time goes now (gprof, 3D): the loops of `vof_advection ()` 47%
(Youngs' stencil, the face loops, the updates); the restrictions of the
ghost updates 32% (`restriction_face` for the fluxes and `uf`,
`restriction_plane` for c, n, alpha), on a uniform tree that reads no
parent; the ghost copies and walls 11%; the geometry 1%.

What is left, outside libvof:
- **Restrict only the parents that are read** (libb-amr): those a leaf of
  their level is beside, and the parents below them. A uniform tree reads
  none; near level jumps a few. A third of the time here.
- **Face loops by rows** (libb-amr's `foreach_face_cells`): each face's
  cell is found by `block_at ()` from an index array, which gcc does not
  vectorise ("complicated access pattern"); rows as `block_cells_simd` has
  them would. libvof's face loop would stay scalar (its interface faces
  call the geometry), but the flux loops of multigrid would gain.
- **A device back end** for `foreach ()` and the ghost updates: until then
  none of this runs on a GPU.

## Next

- Adaptation that follows the interface (`adapt_wavelet ()` on c, with
  `prolong_fraction` refining).
- Height functions and curvature (`vof/vof_curvature.cpp`).
- Tracers carried with c, as Basilisk's `tracers`.
