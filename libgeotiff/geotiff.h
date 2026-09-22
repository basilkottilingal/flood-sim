#ifndef _GEOTIFF_H_
#define _GEOTIFF_H_

  #include "coordinates.h"

  /* APIs */
  void     geotiff_map            ( const char * tiff );
  void     geotiff_map_destroy    ( );
  int      geotiff_elevation      ( CoordL * point_array, double * elevation, int npoints, CoordG c0);   

#endif 
