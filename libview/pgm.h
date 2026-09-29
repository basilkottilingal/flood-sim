#ifndef _VIEW_PGM_H
#define _VIEW_PGM_H

  int pgm_grayscale (double **, const char * out, int width, int height);
  int pgm_contour_grayscale (double **, const char * out, int width, int height, int ncontours);
  /* for easy reference use azimuth = 315 deg, alt = 45 deg, cellsize = 30 m, z_factor = 1. */
  int pgm_hillshade (double **, const char * out, int width, int height, double azimuth, double altitude, double cellsize, double z_factor);

#endif
