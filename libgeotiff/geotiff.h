#ifndef _GEOTIFF_H_
#define _GEOTIFF_H_

  #include "coordinates.h"

  typedef struct
  {
    /*
    .. to define a square patch whose left-top (or rather west-north) corner
    .. is defined by (lon0, lat0) and width & height of the square patch = Delta
    .. (in meters)
    */
    double lon0, lat0, Delta;
    double tiepoint [6];
  } GTBox;

  /* APIs */
  void   geotiff_map            ( const char * tiff );
  void   geotiff_map_destroy    ( );
  double geotiff_elevation_at   ( CoordG c );
  int    geotiff_box            ( GTBox * box );
  //double geotiff_elevation_grid ( CoordL * point_array, int npoints, GTBox * box );   

#endif 
