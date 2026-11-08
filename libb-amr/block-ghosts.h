#ifndef _BAMR_BLOCK_GHOSTS_H_
#define _BAMR_BLOCK_GHOSTS_H_

  /*
  .. Ghost cells, stored with every block.
  ..
  .. A block stores its m^dimension interior cells and GHOSTS layers of ghost
  .. cells around them (tree_block_shape (), block-loops.h). A stencil reads
  .. the ghosts in place, once they are filled:
  ..
  ..   ghosts_update (u, v);          the ghosts of every leaf, for u and v
  ..   foreach_leaf () { ... cell_at (a, -1, 0) ... }
  ..
  .. ghosts_update () does what block-AMR codes do (PARAMESH, Athena++,
  .. Parthenon), in this order:
  ..
  ..   1. restriction, finest level first: each parent gets its scalar's
  ..      restriction of its children (the mean by default);
  ..   2. the ghosts, coarsest level first, each block's from
  ..        - the node of its own level beside it, leaf or parent: a copy;
  ..        - else the leaf one level coarser that covers it (2:1 balance:
  ..          never two): its scalar's prolongation, which reads the coarse
  ..          leaf's interior and its first ghost layer, filled already, as
  ..          that level came first;
  ..        - outside the domain: its scalar's boundary condition, one axis
  ..          after the other, so that corners come out right.
  ..
  .. Within one level the blocks are independent: they run on all threads
  .. with OpenMP.
  ..
  .. Per scalar, three operations: ScalarOps, set with scalar_set_ops (). The
  .. same prolongation fills the children of a refined block and the ghosts
  .. at a level jump, so the two agree.
  ..
  .. Derived scalars. Fill the ghosts of the scalars a derived one is
  .. computed from, then compute it on the ghosts too, halo GHOSTS:
  ..
  ..   ghosts_update (h, zb);
  ..   BlockShape shape = tree_block_shape ();
  ..   foreach_leaf () {
  ..     Real * e = block_data (sblocks, eta);
  ..     const Real * hh = block_data (sblocks, h), * z = block_data (sblocks, zb);
  ..     block_cells (shape, GHOSTS) e [c] = hh [c] + z [c];
  ..   }
  ..
  .. which costs no prolongation of eta, and gives the ghosts of eta that
  .. h and zb have (what Parthenon's FillDerived does).
  ..
  .. Refining and coarsening (tree_refine (), tree_coarsen ()) keep the ghosts
  .. of every leaf valid, if they were: see tree.h.
  */

  #include "block-loops.h"

  /*
  .. A prolongation: fills the cells lo to hi - 1 (interior coordinates of
  .. the fine block, along each axis; ghosts are < 0 or >= m) of the fine
  .. block from the coarse one. Fine cell x lies in the fine cell
  .. X = off + x of the coarse block, counted at the fine level from the
  .. coarse block's first interior cell: its coarse cell is X >> 1 (0 to
  .. m - 1, the coarse block's interior), and X & 1 says which half. The
  .. coarse block's first ghost layer is valid, so X >> 1 +- 1 can be read.
  ..
  .. prolong_cells (pr) below visits every fine cell to fill: write a
  .. prolongation with it. It may read other scalars of the same nodes
  .. (index_node (pr->coarse_index)->s), for prolongations that need them,
  .. as the depth for the water height. It writes only pr->fine, at those
  .. cells: blocks of one level run at once.
  */
  typedef struct
  {
    const Real * coarse;     /* the coarse block, stored               */
    Real * fine;             /* the fine block, stored                 */
    BlockShape shape;        /* of both: tree_block_shape ()           */
    int off [3];             /* see above; 0 beyond the dimension      */
    int lo [3], hi [3];      /* the fine cells to fill; 0, 1 beyond    */
    Index coarse_index, fine_index;
    int s;                   /* the scalar                             */
  } Prolongation;

  /*
  .. A boundary condition: fills the ghost layers of block b outside the
  .. domain on one side (0: low, 1: high) of one axis, at every cell of the
  .. other axes, ghosts included. Those are already filled, ghosts from the
  .. neighbours and the walls of the axes before, so a corner gets the wall
  .. condition of each of its axes in turn. boundary_cells (bc) below visits
  .. every ghost cell to fill: write a condition with it.
  */
  typedef struct
  {
    Real * b;                /* the block, stored                      */
    BlockShape shape;        /* tree_block_shape ()                    */
    int axis, side;
    Index p;                 /* the block's node                       */
    int s;                   /* the scalar                             */
  } Boundary;

  /*
  .. The operations of a scalar. NULL takes the default: linear
  .. prolongation, the mean, a zero gradient at the walls. A scalar gets the
  .. defaults back when it is freed while the tree exists.
  */
  typedef struct
  {
    void (* prolongation) (const Prolongation * pr);
    void (* restriction)  (Index parent, int s);
    void (* boundary)     (const Boundary * bc);
  } ScalarOps;

  void      scalar_set_ops (int s, ScalarOps ops);
  ScalarOps scalar_get_ops (int s);     /* the defaults where NULL */
  void      scalar_ops_reset (int s);   /* back to the defaults */

  /*
  .. The prolongations given:
  ..   prolong_linear     u (C) + sum over axes of sgn / 8 (u (C+1) - u (C-1)):
  ..                      a quarter of a coarse cell times the central
  ..                      gradient. Exact for a linear field, and keeps the
  ..                      mean, as the 2^dimension fine cells of C have both
  ..                      signs on every axis. The default.
  ..   prolong_minmod     the same, the gradient limited (minmod): no new
  ..                      extrema, so a positive field stays positive (a
  ..                      water depth at a wet/dry front). Keeps the mean.
  ..   prolong_injection  u (C): first order, keeps the mean.
  .. The restriction given: restriction_mean, the mean of the 2^dimension
  .. children's cells. The boundary conditions given:
  ..   bc_neumann         the cell mirrored inside: zero gradient. The default.
  ..   bc_odd             minus that: zero at the wall, for the normal
  ..                      velocity at a solid wall.
  */
  void prolong_linear    (const Prolongation * pr);
  void prolong_minmod    (const Prolongation * pr);
  void prolong_injection (const Prolongation * pr);
  void restriction_mean  (Index parent, int s);
  void bc_neumann        (const Boundary * bc);
  void bc_odd            (const Boundary * bc);

  /*
  .. prolong_cells (pr) { ... }: every fine cell pr asks for, with f its
  .. index in pr->fine, cc the index of its coarse cell in pr->coarse, and
  .. sgn [3] the half it is in along each axis, -1 or +1 (0 beyond the
  .. dimension). A neighbour of the coarse cell along axis a is
  .. cc +- block_stride (pr->shape, a). No break.
  ..
  ..   prolong_cells (pr)
  ..     pr->fine [f] = pr->coarse [cc];          // injection
  */
  @loop prolong_cells (const Prolongation * pr) : f, cc, sgn
  {
    int lo [3] = { 0, 0, 0 }, hi [3] = { 1, 1, 1 };
    for (int a = 0; a < dimension; a++)
      lo [a] = pr->lo [a], hi [a] = pr->hi [a];
    int x [3];
    for (x [2] = lo [2]; x [2] < hi [2]; x [2]++)
    for (x [1] = lo [1]; x [1] < hi [1]; x [1]++)
    for (x [0] = lo [0]; x [0] < hi [0]; x [0]++)
    {
      int C [3] = { 0, 0, 0 }, sgn [3] = { 0, 0, 0 };
      for (int a = 0; a < dimension; a++)
      {
        int X = pr->off [a] + x [a];
        assert (X >= 0 && X < 2 * pr->shape.m);
        C [a] = X >> 1;
        sgn [a] = X & 1 ? 1 : -1;
      }
      long f = block_at (pr->shape, x), cc = block_at (pr->shape, C);
      (void) f, (void) cc, (void) sgn;
      {...}
    }
  }

  /*
  .. boundary_cells (bc) { ... }: every ghost cell bc asks for, with ghost
  .. its index in bc->b, inner the index of the cell mirrored inside the
  .. domain, layer 0 for the ghosts next to the wall and up to GHOSTS - 1,
  .. and inward the stride from a cell to the next one towards the inside.
  .. The interior cell next to the wall is ghost + (layer + 1) inward.
  ..
  ..   boundary_cells (bc)                       // u = 0 at the wall face
  ..     bc->b [ghost] = - bc->b [inner];
  */
  @loop boundary_cells (const Boundary * bc) : ghost, inner, layer, inward
  {
    BlockShape sh = bc->shape;
    int m = sh.m, g = sh.g, ax = bc->axis;
    long st = block_stride (sh, ax), inward = bc->side ? -st : st;
    int lo [3] = { 0, 0, 0 }, hi [3] = { 1, 1, 1 };
    for (int a = 0; a < dimension; a++)
      lo [a] = -g, hi [a] = m + g;
    for (int layer = 0; layer < g; layer++)
    {
      lo [ax] = bc->side ? m + layer : -1 - layer;
      hi [ax] = lo [ax] + 1;
      int x [3];
      for (x [2] = lo [2]; x [2] < hi [2]; x [2]++)
      for (x [1] = lo [1]; x [1] < hi [1]; x [1]++)
      for (x [0] = lo [0]; x [0] < hi [0]; x [0]++)
      {
        long ghost = block_at (sh, x), inner = ghost + (2 * layer + 1) * inward;
        (void) ghost, (void) inner, (void) inward;
        {...}
      }
    }
  }

  /*
  .. API
  ..
  ..   ghosts_update (s, ...)         the ghosts of every leaf, for the
  ..                                  scalars given; parents restricted first
  ..   ghosts_update_levels (s, ...)  the same for every node, parents too: for
  ..                                  a stencil on a coarse level (multigrid)
  ..   ghosts_update_scalars (parents) every scalar in use; parents too if
  ..                                  parents is not 0
  ..   ghosts_update_list (n, s, parents)   the function under them all
  ..
  ..   tree_restrict (s)    each parent: its restriction, finest first
  ..   block_ghosts (p, s)  the ghosts of the node p only. The coarser leaves
  ..                        around it must have valid ghosts, and the parents
  ..                        around it their restriction
  ..   block_prolong (p, s) the interiors of p's children, from p
  ..   block_restrict (p, s) p from its children, by the scalar's restriction
  */
  #define ghosts_update(...)                                            \
    ghosts_update_list (sizeof ((int []) {__VA_ARGS__}) / sizeof (int), \
                        (int []) {__VA_ARGS__}, 0)
  #define ghosts_update_levels(...)                                     \
    ghosts_update_list (sizeof ((int []) {__VA_ARGS__}) / sizeof (int), \
                        (int []) {__VA_ARGS__}, 1)

  void ghosts_update_list    (int n, const int * s, int parents);
  void ghosts_update_scalars (int parents);
  void tree_restrict  (int s);
  void block_ghosts   (Index p, int s);
  void block_prolong  (Index p, int s);
  void block_restrict (Index p, int s);

#endif
