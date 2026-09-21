#ifndef _GEOTIFF_H_
#define _GEOTIFF_H_

  #include "coordinates.h"

  /* APIs */
  void   geotiff_map          ( const char * tiff );
  void   geotiff_map_destroy  ( );
  double geotiff_elevation_at ( CoordG c );

#endif 
