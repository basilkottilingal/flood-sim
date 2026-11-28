#ifndef _HMESH_QUALITY_H_
#define _HMESH_QUALITY_H_

  /*
  .. The quality of the mesh: shape measures of its simplices, and their
  .. statistics (min, mean, max, and which simplex has the min and the max).
  .. For remeshing (which edges to split or collapse, which simplices are
  .. bad), for the calculus on the mesh (cotan weights go wrong on obtuse
  .. triangles, curvature on uneven edges), and to watch a mesh as it moves.
  ..
  .. The measures of a k-simplex, k >= 2 (hm_shape ()), each 1 for the
  .. regular simplex and larger the worse it is, as in Verdict and in
  .. Frey and George:
  ..
  ..   edge ratio      the longest edge over the shortest
  ..   aspect ratio    the longest edge over the inradius, scaled:
  ..                   l_max / (2 sqrt (3) r)  triangle, r = 2 A / P
  ..                   l_max / (2 sqrt (6) r)  tetrahedron, r = 3 V / S
  ..   radius ratio    the circumradius over the inradius, scaled:
  ..                   R / (2 r)  triangle, R = a b c / (4 A)
  ..                   R / (3 r)  tetrahedron, R from the products of
  ..                              opposite edges, (p+q+r)(p+q-r)(p-q+r)
  ..                              (-p+q+r) = (24 V R)^2
  ..   angles          the smallest and largest, in degrees: a triangle's
  ..                   interior angles (60 regular), a tetrahedron's
  ..                   dihedral angles (70.53 regular)
  ..
  .. A, P: area and perimeter; V, S: volume and the faces' total area. The
  .. angles by atan2 of the sine and the cosine, accurate near 0 and 180.
  .. Every point is in R^D, so a triangle of a surface in 3D is measured in
  .. its own plane. A simplex whose measure is below HM_QUALITY_EPS times
  .. l_max^k is degenerate: its ratios are INFINITY.
  ..
  .. hm_quality (M, Q): the statistics over the whole mesh, m = M->dim:
  ..
  ..   measure                  the top simplices' length, area or volume
  ..   edge                     every edge's length (level 1)
  ..   edge_ratio, aspect,      over the top simplices, m >= 2; the
  ..   radius_ratio,            degenerate ones counted in degenerate and
  ..   angle_min, angle_max     left out (angle_min is the statistics of
  ..                            each simplex's smallest angle)
  ..   edge_ratio, turning      m = 1, over the points with two edges: the
  ..                            longer edge over the shorter, and the
  ..                            turning angle (0 straight), in degrees
  ..   valence                  the edges at each point, m >= 1
  ..   edge_cell                -DEULERIAN_LAGRANGIAN with a cell size set:
  ..                            each edge over min (cell (a), cell (b)),
  ..                            below 1 if it fits
  ..
  .. An HmStat's argmin and argmax are ids at its level: top simplices,
  .. edges, or points (valence; edge_ratio and turning for m = 1); n is the
  .. number of values, mean their mean; with n = 0 the rest is 0.
  ..
  ..   hm_shape (M, k, s, &q)      the measures of the k-simplex s, k >= 1
  ..                               (k = 1: its length, ratios 1, angles
  ..                               NaN); -1 if s is not alive, or k > m
  ..   hm_quality (M, &Q)          0, or -1 (out of memory)
  ..   hm_quality_print (fp, &Q)   a table, a line a statistic with values
  ..
  .. Plain C on hmesh.h, built as the library is (cwrap -dimension=D, or
  .. -DHM_DIMENSION=D). Geometry in double whatever HmReal is.
  ..
  .. References: P. M. Knupp, Algebraic mesh quality metrics, SIAM J. Sci.
  .. Comput. 23 (2001); C. J. Stimpson, C. D. Ernst, P. Knupp, P. P.
  .. Pebay, D. Thompson, The Verdict geometric quality library, Sandia
  .. report SAND2007-1751 (2007); P. J. Frey, P.-L. George, Mesh
  .. generation, Hermes (2000), chapter 10.
  */

  #include <stdio.h>
  #include <stdint.h>

  #include "hmesh.h"

  #define HM_QUALITY_EPS 1e-12

  typedef struct {
    double   min, mean, max;
    uint32_t argmin, argmax;   /* ids at the statistic's level */
    uint32_t n;                /* values counted */
  } HmStat;

  typedef struct {
    double measure;            /* length, area, volume */
    double lmin, lmax;         /* the shortest and longest edges */
    double edge_ratio, aspect, radius_ratio;
    double angle_min, angle_max;   /* degrees */
  } HmShape;

  typedef struct {
    int      dim;              /* m */
    uint32_t points, simplices;    /* counts of levels 0 and m */
    uint32_t degenerate;       /* top simplices left out of the ratios */
    HmStat   measure, edge;
    HmStat   edge_ratio, aspect, radius_ratio, angle_min, angle_max;
    HmStat   turning, valence, edge_cell;
  } HmQuality;

  int  hm_shape (const Hmesh * M, int k, uint32_t s, HmShape * q);
  int  hm_quality (const Hmesh * M, HmQuality * Q);
  void hm_quality_print (FILE * fp, const HmQuality * Q);

#endif
