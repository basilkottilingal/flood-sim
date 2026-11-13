/*
.. The ghosts, with libb-amr's operators and libvof's (vof-ghosts.h): the
.. file that dispatches to them, in place of libb-amr's block-ghosts.c.
*/
#include "block-ops.h"
#include "vof.h"
#include "vof-ghosts.h"

@prolongation prolong_fraction (const Prolongation * pr)
{
  scalar s = pr->s;
  const double eps = VOF_EPS;
  foreach_prolongation (pr) {
    double f = coarse (s);
    if (f <= eps || f >= 1 - eps)
      s[] = (Real) f;
    else {
      /* Youngs' normal at the coarse cell, minus the gradient (vof.c) */
      coord m = { 0, 0, 0 }, lo = { 0, 0, 0 }, hi = { 1, 1, 1 };
      foreach_dimension () {
#if dimension == 1
        m.x = (double) coarse (s, -1) - coarse (s, 1);
#elif dimension == 2
        m.x = ((double) coarse (s, -1, 1) + 2 * coarse (s, -1)
               + coarse (s, -1, -1))
            - ((double) coarse (s, 1, 1) + 2 * coarse (s, 1)
               + coarse (s, 1, -1));
#else
        m.x = (4 * (double) coarse (s, -1)
               + 2 * ((double) coarse (s, -1, 1) + coarse (s, -1, -1)
                      + coarse (s, -1, 0, 1) + coarse (s, -1, 0, -1))
               + ((double) coarse (s, -1, 1, 1) + coarse (s, -1, 1, -1)
                  + coarse (s, -1, -1, 1) + coarse (s, -1, -1, -1)))
            - (4 * (double) coarse (s, 1)
               + 2 * ((double) coarse (s, 1, 1) + coarse (s, 1, -1)
                      + coarse (s, 1, 0, 1) + coarse (s, 1, 0, -1))
               + ((double) coarse (s, 1, 1, 1) + coarse (s, 1, 1, -1)
                  + coarse (s, 1, -1, 1) + coarse (s, 1, -1, -1)));
#endif
        /* the child's half of the coarse cell along x */
        lo.x = child.x > 0 ? 0.5 : 0;
        hi.x = lo.x + 0.5;
      }
      if (fabs (m.x) + fabs (m.y) + fabs (m.z) == 0)
        s[] = (Real) f;
      else
        s[] = (Real) vof_plane_volume (m, vof_plane_alpha (m, f), lo, hi);
    }
  }
}

/*
.. c with its plane, n and alpha (vof_advection ()): the fine cells from the
.. coarse cell's plane, n the same, alpha in the fine cell's unit cell: a
.. point x of the coarse cell is (o + X) / 2 of the fine cell's X, o the
.. fine cell's corner (0 or 1 along each axis), so n . x <= alpha is
.. n . X <= 2 alpha - n . o.
*/
@prolongation prolong_plane (const Prolongation * pr) (vector n, scalar alpha)
{
  scalar s = pr->s;
  const double eps = VOF_EPS;
  foreach_prolongation (pr) {
    double f = coarse (s), a = 2 * (double) coarse (alpha);
    coord m = { 0, 0, 0 }, lo = { 0, 0, 0 }, hi = { 1, 1, 1 };
    foreach_dimension () {
      m.x = coarse (n.x);
      a -= child.x > 0 ? m.x : 0;
    }
    if (f <= eps || f >= 1 - eps || fabs (m.x) + fabs (m.y) + fabs (m.z) == 0)
      s[] = (Real) f;
    else
      s[] = (Real) vof_plane_volume (m, a, lo, hi);
    foreach_dimension ()
      n.x[] = (Real) m.x;
    alpha[] = (Real) a;
  }
}

/* the parents: c the mean, n the mean of the children's, alpha its plane
 * holding the parent's c */
@restriction restriction_plane (Index p, scalar s) (vector n, scalar alpha)
{
  foreach_restriction (p) {
    RealSum f = 0;
    coord m = { 0, 0, 0 };
    foreach_child () {
      f += s[];
      foreach_dimension ()
        m.x += n.x[];
    }
    f /= TREE_NCHILDREN;
    s[] = (Real) f;
    foreach_dimension () {
      m.x /= TREE_NCHILDREN;
      n.x[] = (Real) m.x;
    }
    alpha[] = (Real) (fabs (m.x) + fabs (m.y) + fabs (m.z) == 0 ? 0
                      : vof_plane_alpha (m, f));
  }
}

#include "block-ghosts-impl.h"
