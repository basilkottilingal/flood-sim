#ifndef _VIEW_PPM_H
#define _VIEW_PPM_H

  /*
  .. Colour (binary PPM, P6) output of 2D rasters.
  .. All functions return 0 on success, -1 on failure.
  .. Values are stretched linearly from [min, max] (finite cells) onto the
  .. colormap; NaN cells are white.
  */

  #include <stdint.h>
  #include "liballoc/grid.h"

  typedef enum
  {
    PPM_VIRIDIS = 0,  /* perceptually uniform, default            */
    PPM_TURBO,        /* rainbow-like, high contrast              */
    PPM_TERRAIN,      /* sea blue - green - yellow - brown - snow */
    PPM_BLUES,        /* white to dark blue (water depth)         */
    PPM_GRAY
  } PPMColormap;

  /* colour of t in [0,1] (clamped; NaN => white) */
  void ppm_color (PPMColormap cmap, double t, uint8_t rgb [3]);

  /* any uint8_t / int / double grid from grid (). Cropped to 1024 x 1024 */
  #define ppm(g,out)           ppm_general ((void *) g, out, grid_data (g), PPM_VIRIDIS)
  #define ppm_cmap(g,out,cmap) ppm_general ((void *) g, out, grid_data (g), cmap)
  int ppm_general  (void * g, const char * out, GridData gd, PPMColormap cmap);

  /*
  .. categorical colours for integer tags (basins, pits, masks ..):
  .. 0 is black, every other value gets a distinct hue.
  */
  #define ppm_labels(g,out)    ppm_labels_general ((void *) g, out, grid_data (g))
  int ppm_labels_general (void * g, const char * out, GridData gd);

  int ppm_colormap (double **, const char * out, int width, int height, PPMColormap cmap);

  /* colormap filled raster with ncontours black iso-lines */
  int ppm_contour  (double **, const char * out, int width, int height,
    int ncontours, PPMColormap cmap);

  /*
  .. shaded relief : elevation coloured with cmap (PPM_TERRAIN suits DEMs)
  .. and modulated by the hillshade. Same parameters as pgm_hillshade ().
  */
  int ppm_hillshade (double **, const char * out, int width, int height,
    double azimuth, double altitude, double cellsize, double z_factor, PPMColormap cmap);

#endif
