# liblinear-algebra

Solvers on libb-amr's block tree, written with cwrap's stencil syntax. For
now one: geometric multigrid for

    div (alpha grad a) + lambda a = b

```c
#include "multigrid.h"

scalar a, b;
...
ghost_ops (a, boundary.all (bc_odd));      /* a = 0 at the walls */
MGStats st = mg_poisson (a, b, (MGParams) { .tolerance = 1e-6 });
/* or: mg_helmholtz (a, b, alpha, lambda, p) */
```

`MGParams` zeroed gives the defaults; `MGStats` has the cycles made and the
largest |residual| before and after. `mg_helmholtz_operator (a, alpha,
lambda, out)` applies the operator as the solver discretises it.

## The discretisation

In flux form. `alpha` is a `face_vector`: `alpha.x` at a cell's lower face
along x. The fluxes `g.x = alpha.x (a[] - a[-1]) / Delta` are made once per
face (`foreach_face ()`), and on each leaf the operator is what they bring
in, over its size, plus `lambda a`:

    foreach_face ()
      g.x[] = alpha.x[] * (a[] - a[-1]) / Delta;
    foreach () {
      double div = 0;
      foreach_dimension ()
        div += g.x[1] - g.x[];
      out[] = div / Delta + lambda[] * a[];
    }

At a level jump a coarse face's flux is the mean of the fine faces' on it
(refluxing, made by libb-amr when `g` is read). So the operator's integral
over the leaves is what crosses the walls: none with zero-gradient walls,
to round-off, level jumps included (`tests/test_poisson.c`, check 0). At a
jump the fine side's flux reads the coarse leaf's ghost, as its
prolongation makes it.

## The cycle

As Basilisk's `mg_cycle ()`:
1. The residual on the leaves, restricted to every level (`tree_restrict`).
2. From the coarsest level up, on the grid at each depth l
   (`foreach_level_or_leaf (l)`: level l and the leaves coarser):
   - the correction `da` set to zero on the coarsest level, or prolonged
     into level l from level l - 1 by its own prolongation
     (`block_prolong`);
   - smoothed, `nrelax` sweeps, its ghosts made after each pass
     (`ghosts_update_level (da, l)`).
3. `da` added to `a` on the leaves.

**The smoother** (`MGParams.smoother`):
- `MG_RED_BLACK`, the default: red-black Gauss-Seidel. The red cells are
  updated in place from the black ones, then the ghosts, the black cells
  from the red ones, and the ghosts again
  (`foreach_level_or_leaf_colour (l, BAMR_RED / BAMR_BLACK)`). cwrap checks
  that a pass reads `da` only at odd offset sums, the other colour, which
  the pass does not change (see libc-pp's README). 2 sweeps by default.
- `MG_JACOBI`: weighted Jacobi. The new values go to a field of their own,
  then are mixed in by `omega` (2d / (2d + 1): 2/3, 4/5, 6/7). 4 sweeps.
- `MG_GAUSS_SEIDEL`: in place, cell after cell within each block
  (`BAMR_BOTH`), the ghosts fixed during the pass. Across blocks that is
  Jacobi, so it falls between the two others.

`nrelax` must be 2 at least: with one sweep a level, the cycle diverges in
2D and 3D for every smoother (measured: a residual ratio of 1.08 to 1.5 per
cycle).

**The coarsest level is a whole block,** m^dimension cells at least, which
the tree cannot coarsen. It gets `nrelax_coarse` sweeps: m^2 / 2 by default
(m^2 for Jacobi, whose spectral radius is the square root of
Gauss-Seidel's).

The design, and why red-black, are in `doc/red-black-gauss-seidel.pdf`
(made by `make` in `doc/red-black/`, not tracked).

## Conditions

- `alpha > 0` and `lambda <= 0`, so that the smoothers converge.
- `a`'s boundary conditions, face by face, are the correction's too: they
  must be homogeneous, as `bc_neumann` (zero gradient) and `bc_odd` (zero
  at the wall) are, and depend on no scalar.
- With zero-gradient walls everywhere and `lambda = 0`, `a` is known up to
  a constant: `b` must have zero mean, or the residual stops at it.
- The tree must not change during a solve. The parents of `alpha` and
  `lambda` are restricted from their leaves, and `alpha`'s ghosts made at
  every node.
- **The tolerance is an absolute residual.** Its floor is round-off in an
  operator scaled by 1 / Delta^2: about epsilon 4d max(alpha) max|a| /
  Delta_min^2. In single precision that is close to what a cycle reaches
  (the tests use twice the floor there, 100 times in double).

## Measured (tests/test_poisson.c)

| | 1D | 2D | 3D |
|---|---|---|---|
| worst residual ratio per cycle (adaptive, alpha.x = 1 + x, lambda = -1; flux form) | 0.25 | 0.39 | 0.41 |
| its own discrete solution recovered, double | 2e-13 | 7e-13 | 5e-13 |
| the operator's integral, zero-gradient walls, random field | 2e-17 | 3e-18 | 3e-18 |
| error ratio against sin pi x ..., uniform level 1 to 2 | 4.00 | 3.99 | 3.99 |

**The smoothers against each other** (`tests/bench_smoothers.c`, -O2, double,
an adaptive tree, from zero to the same tolerance; time in ms):

| flux form | Jacobi, 4 sweeps | red-black, 2 | in-block Gauss-Seidel, 3 |
|---|---|---|---|
| 2D, 268 leaves | 20 cycles, 165 | 20 cycles, **131** | 19 cycles, 202 |
| 3D, 92 leaves (377k cells) | 29 cycles, 2589 | 28 cycles, **1913** | 22 cycles, 1832 |

**Flux form costs about three times the cell-centred form** it replaced
(2D: 131 ms against 46.9; 3D: 1913 against 602), for two reasons:

- **more cycles,** 20 against 14 in 2D, 28 against 20 in 3D. The residual
  is in flux form but the smoother relaxes the cell stencil, so at level
  jumps the two operators differ;
- **a costlier residual,** 1.48 ms against 0.39 (2D, 268 leaves): fluxes
  0.22, divergence 0.31, and **the face ghost update of the fluxes 0.95**.
  It makes every ghost layer of every component, and restricts every
  parent, where the divergence needs only each block's upper face and the
  reflux.

Both can be helped (see Not done yet); what the form buys is the integral
kept at level jumps.

## Build and tests

    make                      # liblinear-algebra.a, 2D, as libb-amr's
    sh tests/run.sh           # 1-3D, double and float, sanitizers

The Makefile is the others' (`../lib.mk`, `DEPS := b-amr utils`): it
builds with the `cwrap` on the PATH, which must know `face_vector` and
`foreach_face ()` (libc-pp of 2026-11-12 or later); `make
CC=../libc-pp/build/cwrap` builds with this tree's own.

The smoothers' benchmark, by hand (any dimension):

    ../libc-pp/build/cwrap -dimension=2 -D_GNU_SOURCE -std=gnu11 -O2 -I. -I.. \
      -I../libb-amr tests/bench_smoothers.c multigrid.c ../libb-amr/tree.c \
      ../libb-amr/grid-scalar.c ../libb-amr/memindex.c ../libb-amr/voxel-memindex.c \
      ../libb-amr/nano-vdb.c ../libb-amr/pool.c ../libb-amr/block-ghosts.c \
      ../libb-amr/adapt.c ../libutils/bits.c -lm -o bench && ./bench

`tests/test_poisson.c`:
- exact for its own discretisation, with each smoother: `b` from the
  operator of a known field on an adaptive tree (level jumps), variable
  `alpha`, `lambda = -1`;
- the rate per cycle, while the residual is well above its floor;
- `itermax` honoured;
- second order against the continuous solution of a Poisson problem.

## Not done yet

- A generic solver for other operators: the Helmholtz form first, as
  agreed; `@relax` and `@residual` kinds dispatched by id, as the ghost
  operators are, when a second operator is needed.
- Lazy ghosts per level: the solver updates them explicitly now
  (`ghosts_update_level`); cwrap could, as for `foreach ()`
  (`libb-amr/multigrid.md`).
- **Flux form's cost:**
  - a face ghost update that makes only the shared faces and the reflux,
    when a loop reads a face vector only at the cell and one face up;
  - a smoother that relaxes the flux-form operator at level jumps, as the
    residual is, to bring the cycles back to the cell-centred form's.

  Both are planned in libb-amr's README, Future plans, "Faces, fluxes and
  vertices" (items 1, 2 and 6).
- Coarsest level: an exact (or in-block multigrid) solve instead of m^2
  sweeps.
