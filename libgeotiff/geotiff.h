#ifndef _GEOTIFF_H_
#define _GEOTIFF_H_

  #include "coordinates.h"

  typedef struct
  {
    float ** dem;
    double tiepoint [6];
    double scale    [3];
  } GeotiffDEM;

  /* APIs */
  void        geotiff_map         ( const char * tiff );
  void        geotiff_map_destroy ( );
  int         geotiff_elevation   ( CoordL * point_array, float * elevation, int npoints, CoordG c0 );
  //Geotiff_DEM geotiff_dem_window  ( CoordG c0, CoordG Delta );
  //void        geotiff_dem_free    ( GeotiffDEM dem );

#endif 
