#ifndef _GEOTIFF_WSG84_CRS_H_
#define _GEOTIFF_WSG84_CRS_H_

  #include "coordinates.h"

  /* WGS84 ellipsoid constants : semi-major axis, flattening, square of ecc  */
  #define WGS84_A   6378137.0
  #define WGS84_F  (1.0 / 298.257223563)
  #define WGS84_E2 (2 * WGS84_F - WGS84_F * WGS84_F)

  typedef struct
  {
    double ** lat;
    double ** lon;
    int       n, m;
  } WGS84_GRID;

  int        wgs84_geodetic_to_ecef   ( GeoCoord gc, Coord3 * c);
  int        wgs84_geodetic_local_enu ( GeoCoord gc, GeoCoord gc0, Coord3 * vec);
  WGS84_GRID wgs84_geodetic_grid      ( double * x, double * y, int n, int m); 
  void       wgs84_geodetic_grid_free ( WGS84_GRID );

#endif
