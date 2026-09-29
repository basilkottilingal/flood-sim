#ifndef _VIEW_PGM_H
#define _VIEW_PGM_H

  int pgm_grayscale (float **, const char * out, int width, int height);
  int pgm_contour_grayscale (float **, const char * out, int width, int height, int ncontours);
  //int pgm_hillshade (float **, const char * out, int width, int height, float azimuth, float altitude, float cellsize, float z_factor);

#endif
