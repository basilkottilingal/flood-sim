#ifndef _LA_MULTIGRID_H_
#define _LA_MULTIGRID_H_

  /*
  .. Geometric multigrid on libb-amr's tree, for
  ..
  ..   div (alpha grad a) + lambda a = b
  ..
  .. cell-centred: on each leaf, sum over its faces of alpha_f (a_n - a) /
  .. Delta^2, plus lambda a, alpha_f the mean of alpha in the two cells of
  .. the face. At a level jump the neighbour is the leaf's ghost, as the
  .. scalar's prolongation makes it: not in flux form, so the integral of
  .. the solution is not kept exactly across level jumps (flux correction is
  .. not done; see README.md).
  ..
  .. The cycle, as Basilisk's mg_cycle (): the residual restricted to every
  .. level; then from the coarsest level up, the correction da set to zero
  .. on it, or prolonged into each finer level by the scalar's own
  .. prolongation, and relaxed on the grid of that depth
  .. (foreach_level_or_leaf (l), block-loops.h) by weighted Jacobi, its
  .. ghosts made after each sweep (ghosts_update_level (), block-ghosts.h);
  .. last, da added to a on the leaves.
  ..
  .. The coarsest level is a block of m^dimension cells at least, which the
  .. tree cannot coarsen: it gets nrelax_coarse sweeps, m^2 by default.
  ..
  .. Conditions:
  ..   - alpha > 0 and lambda <= 0, so that Jacobi converges.
  ..   - a's boundary conditions, face by face, are those of the
  ..     correction too: they must be homogeneous, as bc_neumann (zero
  ..     gradient) and bc_odd (zero at the wall) are, and depend on no
  ..     scalar.
  ..   - With zero-gradient walls everywhere and lambda = 0, a is known up
  ..     to a constant only: b must have zero mean, or the residual stops
  ..     at it.
  ..   - The tree must not change during a solve.
  ..
  .. Fields are read on the leaves (a, b, alpha, lambda); a is the first
  .. guess on entry, the solution on return. The parents of alpha and lambda
  .. are restricted from their leaves, alpha's ghosts made at every node.
  */

  #include "tree.h"

  /*
  .. The smoother, a sweep on a level:
  ..   MG_RED_BLACK     red-black Gauss-Seidel: the red cells in place from
  ..                    the black ones, the ghosts, the black from the red,
  ..                    the ghosts (foreach_level_or_leaf_colour, block-loops.h)
  ..   MG_JACOBI        weighted Jacobi: new values in a field of their own,
  ..                    mixed in by omega, the ghosts
  ..   MG_GAUSS_SEIDEL  Gauss-Seidel in place in each block, cell after cell
  ..                    (BAMR_BOTH), the ghosts: across blocks, Jacobi
  ..   0                the default: MG_RED_BLACK
  */
  enum { MG_SMOOTHER_DEFAULT, MG_RED_BLACK, MG_JACOBI, MG_GAUSS_SEIDEL };

  typedef struct
  {
    int    smoother;        /* see above: 0 gives MG_RED_BLACK             */
    int    nrelax;          /* sweeps on each level of a cycle: 0 gives 2,
                               4 for Jacobi; 2 at least: with 1 the cycle
                               diverges in 2D and 3D                       */
    int    nrelax_coarse;   /* on the coarsest level: 0 gives m^2 / 2, m^2
                               for Jacobi                                  */
    double omega;           /* Jacobi's weight: 0 gives 2d / (2d + 1)      */
    double tolerance;       /* the largest |residual| to stop at: 0 gives
                               1e-3                                        */
    int    itermax;         /* cycles at most: 0 gives 100                 */
    int    minlevel;        /* the coarsest level of a cycle: 0, the root  */
  } MGParams;

  typedef struct
  {
    int    i;               /* cycles made                                  */
    double resb, resa;      /* the largest |residual| before and after      */
  } MGStats;

  /* div (alpha grad a) + lambda a = b; p zeroed for the defaults */
  MGStats mg_helmholtz (scalar a, scalar b, scalar alpha, scalar lambda,
                        MGParams p);

  /* the Laplacian: alpha 1, lambda 0 */
  MGStats mg_poisson (scalar a, scalar b, MGParams p);

  /* out = div (alpha grad a) + lambda a, on the leaves, as the solver
  .. discretises it */
  void mg_helmholtz_operator (scalar a, scalar alpha, scalar lambda,
                              scalar out);

#endif
