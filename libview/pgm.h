#ifndef _VIEW_PGM_H
#define _VIEW_PGM_H

  /*
  .. Grayscale (binary PGM, P5) output of 2D rasters.
  .. All functions return 0 on success, -1 on failure.
  .. Values are stretched linearly from [min, max] (finite cells) to [0, 255];
  .. NaN cells are black.
  */

  #include "liballoc/grid.h"

  /* any uint8_t / int / double grid from grid (). Cropped to 1024 x 1024 */
  #define pgm(g,out) pgm_general ((void *) g, out, grid_data (g))
  int pgm_general    (void * g, const char * out, GridData gd);

  int pgm_grayscale  (double **, const char * out, int width, int height);

  /* ncontours iso-lines, evenly spaced inside (min, max), on white */
  int pgm_contour_grayscale (double **, const char * out, int width, int height, int ncontours);

  /* for easy reference use azimuth = 315 deg, alt = 45 deg, cellsize = 30 m, z_factor = 1. */
  int pgm_hillshade  (double **, const char * out, int width, int height,
    double azimuth, double altitude, double cellsize, double z_factor);

#endif
