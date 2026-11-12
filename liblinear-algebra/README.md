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

Cell-centred, on the leaves: the sum over a cell's faces of alpha_f (a_n -
a) / Delta^2, plus lambda a, alpha_f the mean of alpha in the face's two
cells. At a level jump the neighbour is the leaf's ghost, as the scalar's
prolongation makes it.

Not in flux form: the integral of the solution is not kept exactly across
level jumps. That needs face loops and flux correction (refluxing), which
libb-amr does not have yet.

## The cycle

As Basilisk's `mg_cycle ()`:
1. The residual on the leaves, restricted to every level (`tree_restrict`).
2. From the coarsest level up, on the grid at each depth l
   (`foreach_level_or_leaf (l)`: level l and the leaves coarser):
   - the correction `da` set to zero on the coarsest level, or prolonged
     into level l from level l - 1 by its own prolongation
     (`block_prolong`);
   - relaxed by weighted Jacobi, `nrelax` sweeps (4), its ghosts made after
     each (`ghosts_update_level (da, l)`).
3. `da` added to `a` on the leaves.

Jacobi writes its new values to a field of their own, then mixes them in by
`omega` (2d / (2d + 1): 2/3, 4/5, 6/7), so no loop reads a neighbour it
writes, as cwrap asks.

**The coarsest level is a whole block,** m^dimension cells at least, which
the tree cannot coarsen. It gets `nrelax_coarse` sweeps, m^2 by default:
about a 0.03 reduction there for m = 16.

## Conditions

- `alpha > 0` and `lambda <= 0`, so that Jacobi converges.
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
| worst residual ratio per cycle (adaptive, alpha = 1 + x, lambda = -1) | 0.13 | 0.21 | 0.28 |
| its own discrete solution recovered, double | 7e-12 | 1e-12 | 1e-12 |
| error ratio against sin pi x ..., uniform level 1 to 2 | 4.00 | 3.99 | 3.99 |

## Build and tests

    make DIMENSION=2          # liblinear-algebra.a
    sh tests/run.sh           # 1-3D, double and float, sanitizers

`tests/test_poisson.c`:
- exact for its own discretisation: `b` from the operator of a known field
  on an adaptive tree (level jumps), variable `alpha`, `lambda = -1`;
- the rate per cycle, while the residual is well above its floor;
- `itermax` honoured;
- second order against the continuous solution of a Poisson problem.

## Not done yet

- A generic solver for other operators: the Helmholtz form first, as
  agreed; `@relax` and `@residual` kinds dispatched by id, as the ghost
  operators are, when a second operator is needed.
- Red-black Gauss-Seidel: faster, but cwrap warns on a loop reading a
  neighbour it writes; it needs an explicit opt-out.
- Lazy ghosts per level: the solver updates them explicitly now
  (`ghosts_update_level`); cwrap could, as for `foreach ()`
  (`libb-amr/multigrid.md`).
- Flux form and refluxing at level jumps.
- Coarsest level: an exact (or in-block multigrid) solve instead of m^2
  sweeps.
