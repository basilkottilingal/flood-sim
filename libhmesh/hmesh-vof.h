#ifndef _HMESH_VOF_H_
#define _HMESH_VOF_H_

  /*
  .. Volume fractions from the manifold: Front2VOF (F2V), Pan et al. 2026
  .. ---------------------------------------------------------------------
  ..
  .. J. Pan, D.-A. Koffi Bi, A. B. Kottilingal, S. Costanzo, J. Lu, Y. Ling,
  .. R. Scardovelli, G. Tryggvason, S. Zaleski, Exact computation of the
  .. color function for triangular element interfaces, Int. J. Numer. Meth.
  .. Fluids 98, 1181-1194, 2026, doi:10.1002/fld.70086.
  ..
  .. Only for a manifold that bounds a region: m = D - 1, a curve in 2D or a
  .. surface in 3D (hm_vof () refuses anything else). Its simplices are
  .. oriented with their normals out of the phase (the reference phase,
  .. whose fraction is computed):
  ..
  ..   2D  an edge (v0, v1): normal (t_y, -t_x), t = v1 - v0: a curve going
  ..       counterclockwise about the phase
  ..   3D  a triangle (v0, v1, v2): normal (v1 - v0) x (v2 - v0): its points
  ..       counterclockwise seen from outside the phase
  ..
  .. (for the other phase, 1 - c). The manifold must close the phase, inside
  .. the domain or on its boundary: a closed curve or surface, or one cut to
  .. the domain by hm_clip () (the phase may then lie along the domain's
  .. boundary, as water does along the bottom).
  ..
  .. One cell: the divergence theorem
  .. --------------------------------
  ..
  .. For the region R of the phase in the box [lo, lo + h]^D, with F = (x -
  .. lo_x) e_x, div F = 1:
  ..
  ..   |R| = sum over the pieces P of the manifold in the box of
  ..         A_x(P) (c_x(P) - lo_x)  +  h |R on the face x = lo_x + h|
  ..
  .. A(P) the area vector of the piece (its area times its normal), c(P) its
  .. centroid; the face x = lo_x contributes nothing (x - lo_x = 0), nor the
  .. faces along y and z (n_x = 0). The area on the face is the same
  .. identity one dimension down: the traces of the pieces on the face,
  .. their normals in the face (the pieces' normals projected), and the
  .. length on one edge of the face; the length on an edge, from the points
  .. where the traces cross it, with their signs. (The paper: equations
  .. (7) and (10), Algorithm 3.)
  ..
  .. The pieces: each element clipped to the box (Sutherland-Hodgman, axis by
  .. axis, the paper's Algorithm 2), in the box's unit coordinates. A point
  .. within HM_VOF_EPS (1e-12) of a face of the box is moved onto it, and
  .. flagged with the faces it is on (Algorithm 1): the tests "on a face",
  .. "an edge on a face", "a point on an edge" are then on bits, exact. A
  .. piece lying on a face of the box (every point on it) is left out, as
  .. are traces on an edge of the face and points at a corner: each would be
  .. counted twice.
  ..
  .. What the pieces cannot tell. A face, or an edge, that no piece reaches
  .. is wholly in the phase or wholly out. On an edge the points tell it
  .. exactly: the corner at its end is in the phase if the last crossing
  .. before it leaves the phase behind (its sign). Above, when the face (or
  .. the edge) below gave nothing, its term is taken as 0, and a negative
  .. result means it was full: 1 is added (the paper's F_v). A box with no
  .. piece of positive measure strictly inside is "unknown": its fraction is
  .. 0 or 1, and the pieces cannot tell which.
  ..
  ..   hm_vof_box (M, lo, h, elem, n, &known)   the fraction of the box
  ..                                            [lo, lo + h]^D, from the
  ..                                            elements elem [0 ... n - 1]
  ..                                            (ids of level m), all that
  ..                                            may touch it; known 0 for a
  ..                                            box no piece is inside
  ..   hm_vof_face (M, lo, h, elem, n, a, s,    the fraction of its face
  ..                &known)                     along axis a, side s (0 the
  ..                                            lower, 1 the upper)
  ..   hm_vof_inside (M, x, elem, n)            the paper's Algorithm 4: 1
  ..                                            if x is in the phase, by the
  ..                                            side of the element closest
  ..                                            to x it is on; -1 for no
  ..                                            element. Right when the
  ..                                            closest point of the whole
  ..                                            manifold is among them.
  ..
  .. On the tree: hm_vof (M, c)
  .. --------------------------
  ..
  .. c, a scalar of libb-amr, gets the fraction of every leaf cell.
  ..
  ..   - each element is binned to the leaves its bounding box meets (the box
  ..     halved until each part's corners are in one leaf), then within a
  ..     leaf to the cells its bounding box meets
  ..   - each cell with elements: hm_vof_box ()
  ..   - the cells left unknown (no piece inside: wholly in or out) take
  ..     the state of a face they share with a cell that is cut: the exact
  ..     fraction of that face (hm_vof_face (), on the part of it the unknown
  ..     cell has, at a level jump), which no element crosses, so 0 or 1; or,
  ..     if that face is unknown too, hm_vof_inside () at its centre. Then
  ..     from cell to cell, through faces no element touches, across leaves
  ..     and levels (a breadth-first walk). This replaces the paper's
  ..     Algorithm 4 on the grid: it needs no search for the closest element.
  ..   - a region of cells no cut cell borders (no manifold in the domain)
  ..     is left at 0; the count of such cells is returned.
  ..
  .. Returns that count, or -1 if m is not D - 1. The elements outside the
  .. domain are left out.
  */

  #include "hmesh.h"

  #define HM_VOF_EPS 1e-12

  double hm_vof_box (const Hmesh * M, const double * lo, double h,
                     const uint32_t * elem, uint32_t n, int * known);
  double hm_vof_face (const Hmesh * M, const double * lo, double h,
                      const uint32_t * elem, uint32_t n, int a, int s,
                      int * known);
  int    hm_vof_inside (const Hmesh * M, const double * x,
                        const uint32_t * elem, uint32_t n);
  long   hm_vof (Hmesh * M, scalar c);

#endif
