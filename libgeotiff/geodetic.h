#ifndef _GEOTIFF_GEODETIC_H_
#define _GEOTIFF_GEODETIC_H_

  #include "coordinates.h"
  #include "geodetic.h"

  /*
  .. APIs to convert coordinates b/w different Coordinate systems
  .. a. geodetic : (lon, lat, alt) assumes lat is w.r.t to WSG64 datum (PrimeMeridian = 0
  ..    is Greenwich and lon = 0 is obviously the equatorial plane. Do not confuse
  ..    alt with elevation, it rather means the vertical distance from the ellipsoid
  ..    and thus alt = 0 may not fall on the surface of the earth)
  .. b. ecef means Earth Centered Earth Focused coordinate system 
  ..    whose origin is centered at the centroid of the WSG64 ellipsoid.
  .. c. local enu means Local East-North-Up coordinate system centered around
  ..    the provided reference point c0. (If alt = 0, then local East-North
  ..    plane is (ofcourse) a tangent plane of the WSG ellipsoid.)
  */

  Coord3   geodetic_to_ecef      (GeographicCRS * gcrs, CoordG c);
  CoordL   ecef_to_local_enu     (GeographicCRS * gcrs, Coord3 x, CoordG c0);
  CoordL   geodetic_to_local_enu (GeographicCRS * gcrs, CoordG c, CoordG c0);
  Coord3   local_enu_to_ecef     (GeographicCRS * gcrs, CoordL l, CoordG c0);
  CoordG   ecef_to_geodetic      (GeographicCRS * gcrs, Coord3 x);
  CoordG   local_enu_to_geodetic (GeographicCRS * gcrs, CoordL l, CoordG c0);
#endif 
