#ifndef _HMESH_IMPLICIT_H_
#define _HMESH_IMPLICIT_H_

  /*
  .. The manifold from a function: the interface f = 0
  .. --------------------------------------------------
  ..
  .. f (x, y) in 2D, f (x, y, z) in 3D: f < 0 inside the phase, f > 0
  .. outside. hm_from_function () adds to M (m = D - 1: a curve in 2D, a
  .. surface in 3D) the interface f = 0 inside the box [lo, lo + size]^D,
  .. sampled on a lattice of n cells a side, h = size / n:
  ..
  ..   - each lattice cell split into simplices the same way everywhere:
  ..     2D a square in two triangles along its diagonal (0,0)-(1,1), 3D a
  ..     cube in six tetrahedra about (0,0,0)-(1,1,1) (Kuhn); neighbours
  ..     split their shared faces alike, so the interface is conforming, and
  ..     no case is ambiguous (as marching squares and cubes are)
  ..   - a point on each lattice edge whose ends differ in sign (f <= 0 is
  ..     in), found once and shared by every simplex on that edge: the root
  ..     of f along the edge, by bracketed regula falsi (Illinois), to
  ..     round-off; a root at a lattice node (f = 0 there) is that node's
  ..     point, shared by its edges, and a simplex it would make flat is
  ..     left out
  ..   - in each simplex of the lattice, the interface's piece: a segment
  ..     (2D), a triangle or a quadrangle in two (3D), oriented with its
  ..     normal towards f > 0 (out of the phase, as hm_vof () wants)
  ..
  .. An edge of the interface is in one lattice simplex: at most sqrt (2) h
  .. in 2D, sqrt (3) h in 3D. So h = the cell size / 2 keeps every edge
  .. shorter than the cells (-DEULERIAN_LAGRANGIAN): hm_tree_from_function ()
  .. (hmesh-tree.h) takes the box of the tree's domain and h half its finest
  .. cell, then binds the mesh to the tree.
  ..
  .. Where the phase reaches the box, the interface ends on its faces, as a
  .. mesh cut by hm_clip () does: the phase lies along the boundary there.
  .. Features smaller than h may be missed (a bubble between lattice nodes).
  .. f is evaluated (n + 1)^D times, a layer of nodes at a time.
  ..
  .. Returns the number of simplices of level m added, or -1 (m is not
  .. D - 1, n < 1, a point refused: outside M's domain, or an edge too long
  .. with -DEULERIAN_LAGRANGIAN).
  */

  #include "hmesh.h"

  #if HM_DIMENSION == 1
    typedef double (* HmFunction) (double x);
  #elif HM_DIMENSION == 2
    typedef double (* HmFunction) (double x, double y);
  #else
    typedef double (* HmFunction) (double x, double y, double z);
  #endif

  long hm_from_function (Hmesh * M, HmFunction f, const double * lo,
                         double size, int n);

#endif
