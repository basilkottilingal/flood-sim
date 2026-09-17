#ifndef _GEOTIFF_H_
#define _GEOTIFF_H_

  struct Image;

  /* APIs */
  void   geotiff_map          ( const char * tiff );
  void   geotiff_map_destroy  ( );
  double geotiff_elevation_at ( double geo_coord [/* longi, lati */] );

#endif 
