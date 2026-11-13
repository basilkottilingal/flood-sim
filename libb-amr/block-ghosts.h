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
  ..      restriction of its children (the mean by default); the parents
  ..      of the coarsest leaf's level and finer only, as a parent is read
  ..      by a leaf of its level beside it (every parent for
  ..      ghosts_update_levels (), which makes their ghosts too);
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
  .. Per scalar, three operations: a prolongation and a restriction
  .. (ScalarOps, set with scalar_set_ops ()), and a boundary condition per
  .. face of the domain (ghost_ops (), boundary.left (op), ...), anywhere
  .. and at any time. The same prolongation fills the children of a
  .. refined block and the ghosts at a level jump, so the two agree.
  ..
  .. An operation is an integer, not a function pointer: cwrap's
  .. @prolongation, @restriction and @boundary (libc-pp/README.md) define
  .. them as inline functions, each with an id, and ghosts_update () calls
  .. the one a scalar has through a switch, each case inlined. The file that
  .. dispatches must see every definition: block-ghosts.c, which includes the
  .. ones given here (block-ops.h) and the core (block-ghosts-impl.h). A
  .. program with operators of its own builds its own such file instead:
  ..
  ..   #include "block-ops.h"             // the given ones, ids first
  ..   @boundary bc_inflow (const Boundary * bc) { ... }
  ..   #include "block-ghosts-impl.h"     // the core, dispatching to all
  ..
  .. and any other file attaches one by its name, an int:
  ..
  ..   @boundary bc_inflow;
  ..   ghost_ops (h, boundary.left (bc_inflow));
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
    scalar s;                /* the scalar                             */
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
    scalar s;                /* the scalar                             */
  } Boundary;

  /*
  .. The operations of a scalar, by their ids (see above). 0 takes the
  .. default: linear prolongation, the mean, a zero gradient at the walls. A
  .. scalar gets the defaults back when it is freed while the tree exists.
  .. Setting any marks s's ghosts stale: those made by the operations it had
  .. are made again by the new ones, as foreach () next needs them.
  ..
  .. A boundary condition is bound per face of the domain, with ghost_ops ()
  .. (cwrap): boundary.left (op), .right, .bottom, .top, .back, .front (as
  .. the dimension has them), or boundary.all (op) for every face:
  ..
  ..   ghost_ops (u.x, boundary.all (bc_neumann), boundary.left (bc_odd),
  ..                   boundary.right (bc_odd));
  ..
  .. ScalarOps and scalar_set_ops () hold the prolongation and restriction;
  .. scalar_get_boundary (s, face) reads a face's condition.
  .. Their signatures:
  ..   @prolongation name (const Prolongation * pr) { ... }
  ..   @restriction  name (Index parent, scalar s)   { ... }
  ..   @boundary     name (const Boundary * bc)      { ... }
  ..
  .. Coupled operators: one that also reads other scalars of the same nodes
  .. names them in a second parameter list, as the water depth h prolonged
  .. from the surface h + zb, the bed zb prolonged first:
  ..
  ..   @prolongation prolong_lake (const Prolongation * pr) (scalar zb) { ... }
  ..   ghost_ops (h, prolongation (prolong_lake, zb));
  ..
  .. What an operator depends on may be scalars, vertex scalars, vectors,
  .. face vectors or tensors, declared so: (vector n, scalar alpha). It gets
  .. them as declared; the binding flattens them to their scalars, n.x, n.y,
  .. alpha (GHOST_USES_MAX at most), with their types (the signature), and
  .. a binding of other types is refused: by cwrap, if the operator is
  .. declared so in the file, @prolongation f (const Prolongation * pr)
  .. (vector n, scalar alpha); and by ghost_ops_bind () as it runs.
  ..
  .. ghost_ops () (cwrap) binds an operator and the scalars it depends on,
  .. which must be in use. It checks, as it runs, that the
  .. operator takes that many, and that no scalar's operators depend on it
  .. itself, through others: either stops the run with a message (unless
  .. ghost_ops_abort is 0: then ghost_ops_set () returns 1 and binds
  .. nothing). From then on:
  ..   - a scalar's ghosts are made after those of the scalars it depends
  ..     on, wherever several are made (ghosts_update (), refining,
  ..     coarsening), and a stale one it depends on is updated with it;
  ..   - ghosts_mark_dirty (zb) marks those that depend on zb too.
  .. The scalars depended on may be local: freeing one that a scalar in use
  .. depends on stops the run (ghost_ops_freeing (), from the free hook of
  .. grid-scalar.c). Bind the dependent's operators again first: to bind a
  .. coupled operator for a while, keep scalar_get_ops (s) before, and give
  .. it back with scalar_restore_ops () after.
  .. scalar_set_ops () binds a prolongation and a restriction that depend
  .. on nothing.
  */
  #define GHOST_USES_MAX 16
  enum { GHOST_PROLONGATION, GHOST_RESTRICTION, GHOST_BOUNDARY };

  /* the faces of the domain, those of the dimension only: 2 axis + side,
  .. side 0 the low one; GHOST_FACES of them */
  enum { GHOST_LEFT, GHOST_RIGHT,
  #if dimension > 1
         GHOST_BOTTOM, GHOST_TOP,
  #endif
  #if dimension > 2
         GHOST_BACK, GHOST_FRONT,
  #endif
         GHOST_FACES, GHOST_ALL_FACES = -1 };

  /* what any scalar's operations depend on, at most: see ghost_ops () */
  #define GHOST_DEPS_MAX ((2 + GHOST_FACES) * GHOST_USES_MAX)

  typedef struct
  {
    int prolongation, restriction;
    scalar prolongation_uses [GHOST_USES_MAX];  /* what each depends on */
    scalar restriction_uses  [GHOST_USES_MAX];
  } ScalarOps;

  void      scalar_set_ops (scalar s, ScalarOps ops);
  ScalarOps scalar_get_ops (scalar s);  /* the defaults where 0 */
  int       scalar_get_boundary (scalar s, int face);  /* the default where
                                                          none is bound */
  void      scalar_ops_reset (scalar s);  /* back to the defaults */
  int       scalar_restore_ops (scalar s, ScalarOps ops);
                                        /* a prolongation and a restriction
                                           as scalar_get_ops () gave them,
                                           what they depend on with them:
                                           to bind them again, after a
                                           ghost_ops () for a while */
  void      ghost_ops_freeing (scalar s);  /* s is being freed: none may
                                              depend on it */
  int       ghost_ops_set (scalar s, int kind, int face, int op,
                           const scalar * uses, int n);
                                        /* face only for GHOST_BOUNDARY,
                                           else GHOST_ALL_FACES; the types
                                           of uses not checked */
  int       ghost_ops_bind (scalar s, int kind, int face, int op,
                            const scalar * uses, int n, const char * sig);
                                        /* what ghost_ops () calls: as
                                           ghost_ops_set (), uses the
                                           scalars of the dependencies
                                           flattened, sig their types (s
                                           scalar, x vertex_scalar, v
                                           vector, f face_vector, t tensor),
                                           checked against the operator's
                                           (NULL: not checked) */
  extern int ghost_ops_abort;           /* 1, the default: stop on errors */

  /*
  .. The scalars, ordered so that each comes after those it depends on:
  .. ghost_order (n, in, out) of in [0..n), with the stale ones they depend
  .. on added; ghost_order_all (out) of every scalar in use. Each returns
  .. how many. foreach_scalar_ordered (s) { ... } runs over the second.
  */
  int ghost_order (int n, const int * in, int * out);
  int ghost_order_all (int * out);
  #define foreach_scalar_ordered(s)                                       \
    for (int __cwrap_ol [NSCALAR_MAX], __cwrap_on = ghost_order_all (__cwrap_ol), \
             __cwrap_oq = 0, s = 0;                                         \
         __cwrap_oq < __cwrap_on && ((s = __cwrap_ol [__cwrap_oq]), 1);     \
         __cwrap_oq++)

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
  .. And two that do nothing:
  ..   prolong_none       the fine cells left as they are: for a scalar that
  ..                      another's prolongation fills, a coupled one that
  ..                      writes the scalars it depends on, at the same fine
  ..                      cells (a volume fraction's, filling its normal and
  ..                      plane with it). Its own ghosts at level jumps, and
  ..                      its children, are then made only with that scalar's:
  ..                      update that one to have them. It is never called:
  ..                      refinement, a ghost update at level jumps and
  ..                      block_prolong () skip it; adapt_wavelet () and
  ..                      block_prolong_into () refuse such a scalar, which
  ..                      has no wavelet of its own.
  ..   restriction_none   the parents left as they are, for a scalar whose
  ..                      parents nothing reads, or another's restriction
  ..                      fills. Never called either: no parent is visited
  ..                      for it (ghost updates, tree_restrict (),
  ..                      coarsening, block_restrict ()).
  */
  @prolongation prolong_linear;
  @prolongation prolong_minmod;
  @prolongation prolong_injection;
  @restriction  restriction_mean;
  @boundary     bc_neumann;
  @boundary     bc_odd;
  @prolongation prolong_none;
  @restriction  restriction_none;

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
  /* the coarse cell C and the half sgn of the fine cell x, along axis a:
  .. in range, as prolongation_in_range () checked for the whole region */
  #define BAMR_FINE_TO_COARSE(off, x, C, sgn, m, a) do {                    \
      int X_ = (off) [a] + (x) [a];                                       \
      (C) [a] = X_ >> 1;                                                  \
      (sgn) [a] = X_ & 1 ? 1 : -1;                                        \
    } while (0)

  /* every fine cell pr asks for has its coarse cell in the coarse block:
  .. off + x within 0 .. 2 m - 1 along each axis, for all of lo <= x < hi,
  .. so at its two ends; checked once per call, not per cell */
  static inline int prolongation_in_range (const Prolongation * pr)
  {
    for (int a = 0; a < dimension; a++)
      if (pr->lo [a] < pr->hi [a]
          && (pr->off [a] + pr->lo [a] < 0
              || pr->off [a] + pr->hi [a] > 2 * pr->shape.m))
        return 0;
    return 1;
  }

  @loop prolong_cells (const Prolongation * pr) : f, cc, sgn
  {
    assert (prolongation_in_range (pr));
    const int * lo = pr->lo, * hi = pr->hi;
    int x [3] = { 0, 0, 0 };
  #if dimension > 2
    for (x [2] = lo [2]; x [2] < hi [2]; x [2]++)
  #endif
  #if dimension > 1
    for (x [1] = lo [1]; x [1] < hi [1]; x [1]++)
  #endif
    for (x [0] = lo [0]; x [0] < hi [0]; x [0]++)
    {
      int C [3] = { 0, 0, 0 }, sgn [3] = { 0, 0, 0 };
      BAMR_FINE_TO_COARSE (pr->off, x, C, sgn, pr->shape.m, 0);
  #if dimension > 1
      BAMR_FINE_TO_COARSE (pr->off, x, C, sgn, pr->shape.m, 1);
  #endif
  #if dimension > 2
      BAMR_FINE_TO_COARSE (pr->off, x, C, sgn, pr->shape.m, 2);
  #endif
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
    int lo [3] = { -g, -g, -g }, hi [3] = { m + g, m + g, m + g };
    for (int layer = 0; layer < g; layer++)
    {
      lo [ax] = bc->side ? m + layer : -1 - layer;
      hi [ax] = lo [ax] + 1;
      int x [3] = { 0, 0, 0 };
  #if dimension > 2
      for (x [2] = lo [2]; x [2] < hi [2]; x [2]++)
  #endif
  #if dimension > 1
      for (x [1] = lo [1]; x [1] < hi [1]; x [1]++)
  #endif
      for (x [0] = lo [0]; x [0] < hi [0]; x [0]++)
      {
        long ghost = block_at (sh, x), inner = ghost + (2 * layer + 1) * inward;
        (void) ghost, (void) inner, (void) inward;
        {...}
      }
    }
  }

  /*
  .. Stencil loops over two levels, for the operators (cwrap: see
  .. libc-pp/README.md), written as foreach () is:
  ..
  ..   foreach_prolongation (pr) { ... }  every fine cell pr asks for: s[] is
  ..                     the fine cell, coarse (s), coarse (s, i, j) its coarse
  ..                     cell and that cell's neighbours, child.x, .y, .z the
  ..                     half it is in (-1 or +1, 0 beyond the dimension)
  ..   foreach_restriction (p) { ... }   every cell of the parent p, child by
  ..                     child: s[] is the parent cell, and foreach_child ()
  ..                     moves the point to each of its children, in order
  ..
  .. In both, no offset on the point's own level (s[1]): those cells are
  .. being filled. coarse () reaches BAMR_COARSE_REACH at most: only the
  .. first ghost layer of the coarse block is valid when a level jump is
  .. prolonged.
  */
  enum { BAMR_COARSE_REACH = 1 };

  /* the fine node's blocks, and the coarse one's */
  @loop prolong_block (const Prolongation * pr) : sblocks, coarse_sblocks
  {
    void *** sblocks = index_node (pr->fine_index)->s;
    void *** coarse_sblocks = index_node (pr->coarse_index)->s;
    (void) sblocks, (void) coarse_sblocks;
    {...}
  }

  /* prolong_cells (pr), the fine cell c, its coarse cell cc, sgn its half */
  @loop prolong_points (const Prolongation * pr) : c, shape, cc, sgn
  {
    BlockShape shape = pr->shape;
    assert (prolongation_in_range (pr));
    const int * lo = pr->lo, * hi = pr->hi;
    int x [3] = { 0, 0, 0 };
  #if dimension > 2
    for (x [2] = lo [2]; x [2] < hi [2]; x [2]++)
  #endif
  #if dimension > 1
    for (x [1] = lo [1]; x [1] < hi [1]; x [1]++)
  #endif
    for (x [0] = lo [0]; x [0] < hi [0]; x [0]++)
    {
      int C [3] = { 0, 0, 0 }, sgn [3] = { 0, 0, 0 };
      BAMR_FINE_TO_COARSE (pr->off, x, C, sgn, shape.m, 0);
  #if dimension > 1
      BAMR_FINE_TO_COARSE (pr->off, x, C, sgn, shape.m, 1);
  #endif
  #if dimension > 2
      BAMR_FINE_TO_COARSE (pr->off, x, C, sgn, shape.m, 2);
  #endif
      long c = block_at (shape, x), cc = block_at (shape, C);
      (void) c, (void) cc, (void) sgn;
      {...}
    }
  }

  /* the parent's blocks, and each child's, child by child */
  @loop family_children (Index parent) : sblocks, child_sblocks, child
  {
    void *** sblocks = index_node (parent)->s;
    for (int child = 0; child < TREE_NCHILDREN; child++)
    {
      void *** child_sblocks = index_node (index_child (parent, child))->s;
      (void) sblocks, (void) child_sblocks;
      {...}
    }
  }

  /* the parent cells c of child's quadrant, by rows, and the first of
  .. their children fc, in child's block */
  @loop quadrant_cells (BlockShape shape, int child) : c, shape, fc
  {
    int h = shape.m / 2;
  #if dimension > 2
    for (int Z = 0; Z < h; Z++)
  #else
    const int Z = 0;
  #endif
  #if dimension > 1
    for (int Y = 0; Y < h; Y++)
  #else
    const int Y = 0;
  #endif
    {
      int px [3] = { (child & 1) * h, ((child >> 1) & 1) * h + Y,
                     ((child >> 2) & 1) * h + Z };
      int fx [3] = { 0, 2 * Y, 2 * Z };
      long prow = block_at (shape, px), frow = block_at (shape, fx);
      for (int X = 0; X < h; X++)
      {
        long c = prow + X, fc = frow + 2 * X;
        (void) c, (void) fc;
        {...}
      }
    }
  }

  @foreach foreach_prolongation (const Prolongation * pr)
    = prolong_block (pr) prolong_points (pr);
  @foreach foreach_restriction (Index p)
    = family_children (p) quadrant_cells (tree_block_shape (), child);

  /*
  .. API (and, below, the stale ghosts of foreach ())
  ..
  ..   ghosts_update (s, ...)         the ghosts of every leaf, for the
  ..                                  scalars given; parents restricted first
  ..   ghosts_update_levels (s, ...)  the same for every node, parents too: for
  ..                                  a stencil on a coarse level (multigrid)
  ..   ghosts_update_scalars (parents) every scalar in use; parents too if
  ..                                  parents is not 0
  ..   ghosts_update_list (n, s, parents)   the function under them all
  ..
  ..   ghosts_update_level (s, l)    the ghosts of the grid at depth l
  ..                        (foreach_level_or_leaf (l), block-loops.h), for
  ..                        multigrid: the parents coarser than l restricted
  ..                        from that grid, not from the leaves; s alone
  ..   tree_restrict (s)    each parent: its restriction, finest first
  ..   block_ghosts (p, s)  the ghosts of the node p only. The coarser leaves
  ..                        around it must have valid ghosts, and the parents
  ..                        around it their restriction
  ..   block_prolong (p, s) the interiors of p's children, from p
  ..   block_prolong_into (p, s, target)   p's interior from its parent,
  ..                        into target (a block), p's own left as it is:
  ..                        what p would be if its parent were refined
  ..   block_restrict (p, s) p from its children, by the scalar's restriction
  */
  #define ghosts_update(...)                                            \
    ghosts_update_list (sizeof ((int []) {__VA_ARGS__}) / sizeof (int), \
                        (int []) {__VA_ARGS__}, 0)
  #define ghosts_update_levels(...)                                     \
    ghosts_update_list (sizeof ((int []) {__VA_ARGS__}) / sizeof (int), \
                        (int []) {__VA_ARGS__}, 1)

  void ghosts_update_list    (int n, const int * s, int parents);
  void ghosts_update_level   (int s, int l);
  void ghosts_update_scalars (int parents);

  /*
  .. Stale ghosts, for cwrap's foreach (): a bit per scalar, set while its
  .. ghosts may be out of date. foreach () marks the scalars it writes, and
  .. before a foreach () that reads a scalar at an offset (s[1], not s[]),
  .. updates those of them that are marked, in one ghosts_update_list ().
  .. Every update clears the bits of the scalars it updates, and a new or
  .. freed scalar has its bit clear. Code that writes cells by hand, through
  .. block_data (), marks nothing: call ghosts_update () as before, or
  .. ghosts_mark_dirty ().
  ..
  ..   ghosts_mark_dirty (s)        s's ghosts may be stale
  ..   ghosts_mark_all_dirty ()     every scalar's may be
  ..   ghosts_is_dirty (s)
  ..   ghosts_update_dirty (n, s)   updates those of s [0..n) that are marked
  ..                                (and in use), at once; n < 0: every
  ..                                scalar marked
  */
  #define GHOSTS_DIRTY_WORDS ((NSCALAR_MAX + 63) / 64)
  extern unsigned long long ghosts_dirty [GHOSTS_DIRTY_WORDS];
  /* the scalars that depend on s, through others too: see ghost_ops () */
  extern unsigned long long ghost_dependents [NSCALAR_MAX][GHOSTS_DIRTY_WORDS];

  static inline void ghosts_mark_dirty (int s)
  {
    ghosts_dirty [s / 64] |= 1ull << (s % 64);
    for (int w = 0; w < GHOSTS_DIRTY_WORDS; w++)
      ghosts_dirty [w] |= ghost_dependents [s][w];
  }

  static inline void ghosts_mark_all_dirty (void)
  {
    for (int w = 0; w < GHOSTS_DIRTY_WORDS; w++)
      ghosts_dirty [w] = ~0ull;
  }

  static inline int ghosts_is_dirty (int s)
  {
    return (ghosts_dirty [s / 64] >> (s % 64)) & 1;
  }

  void ghosts_update_dirty (int n, const int * s);
  /* those of s [0..n) stale and not at the centre (faces, vertices): for a
  .. read at the cell (cwrap) */
  void ghosts_update_dirty_located (int n, const int * s);
  void tree_restrict  (int s);
  void block_ghosts   (Index p, int s);
  void block_prolong  (Index p, int s);
  void block_prolong_into (Index p, int s, Real * target);
  void block_restrict (Index p, int s);

#endif
