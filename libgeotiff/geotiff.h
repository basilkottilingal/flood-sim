#ifndef _GEOTIFF_H_
#define _GEOTIFF_H_

  #include "coordinates.h"

  typedef struct
  {
    /*
    .. A dem raster with the index range : raster [-2 : n+1][-2 : m+1].
    .. i.e ghost layer of 2 cell thickness. A new tiepoint will be 
    .. defined for this raster grid.
    */
    float ** raster;
    double tiepoint [6];
    double scale    [3];
    int n, m;
  } DEM;

  /* APIs */
  void  geotiff_map         ( const char * tiff );
  void  geotiff_map_destroy ( );
  int   geotiff_elevation   ( CoordL * point_array, float * elevation, int npoints, CoordG c0 );
  int   geotiff_dem_window  ( double lon, double lat, double dlon, double dlat, DEM * dem );
  void  geotiff_dem_free    ( DEM dem );

#endif 
