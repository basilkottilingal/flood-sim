#ifndef _VIEW_IMAGE_H
#define _VIEW_IMAGE_H

  /*
  .. Shared helpers of pgm.c and ppm.c . Not meant as a user API.
  ..
  .. Rasters are row pointers (raster [y][x]), as returned by grid ().
  .. Flat buffers returned by the helpers below are row-major (w*h),
  .. allocated with malloc () and to be released with free ().
  */

  #include <stddef.h>
  #include <stdint.h>

  /* pgm_general () / ppm_general () crop larger grids to this size */
  #define IMAGE_MAX_SIZE 1024

  /* 1 if arguments are usable : non-NULL raster & out, w,h >= minsize */
  int image_check (const void * raster, const char * out, int w, int h, int minsize);

  /* value of cell (x,y) of a uint8_t / int / double grid (by element size) */
  double image_value (const void * g, size_t size, int x, int y);

  /* min/max over finite cells. Returns the number of finite cells */
  long image_limits (const void * g, size_t size, int w, int h, double * min, double * max);

  /* (v - min)/(max - min) in [0,1] for each cell. NaN for non-finite cells */
  double * image_normalize (const void * g, size_t size, int w, int h);

  /*
  .. Horn (1981) hillshade. Illumination in [0,1], NaN where dem is NaN.
  .. Out of range arguments fall back to azimuth 315, altitude 45, cellsize 30.
  .. Border cells use clamped neighbours (ghost cells are not read).
  */
  double * image_hillshade (double ** dem, int w, int h, double azimuth,
    double altitude, double cellsize, double z_factor);

  /*
  .. Marching squares. Returns, for each pixel, the index of the contour
  .. level drawn on it, or -1. Levels are evenly spaced strictly inside the
  .. finite range of the raster. Returns NULL on failure.
  */
  int * image_contour (double ** raster, int w, int h, int nlevels);

  /* writes binary netpbm : channels 1 (P5, pgm) or 3 (P6, ppm). 0 on success */
  int image_write (const char * out, int channels, int w, int h, const uint8_t * pixels);

#endif
