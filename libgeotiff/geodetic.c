#include <stdio.h>
#include <math.h>

#include "geodetic.h"
#include "coordinates.h"

#define _ec2(gcrs)     ( gcrs->EccentricitySquared )
#define _lat(gcrs,l)   ( l*gcrs->AngularUnitSize )
#define _lon(gcrs,l)   ( (l + gcrs->PrimeMeridianLong) * gcrs->AngularUnitSize )

Coord3 geodetic_to_ecf (GeographicCRS * gcrs, CoordG c)
{
  double lat     = _lat (gcrs, c.lat);
  double lon     = _lon (gcrs, c.lon);
  double sin_lat = sin (lat);
  double cos_lat = cos (lat);
  double e2      = _ec2 (gcrs);
  /* Radius of curvature in the prime vertical */
  double N       = gcrs->SemiMajorAxis / sqrt (1.0 - e2 * sin_lat * sin_lat);

  return
    (Coord3)
      {
        .x = (N + alt) * cos_lat * cos (lon),
        .y = (N + alt) * cos_lat * sin (lon),
        .z = (N * (1.0 - e2) + alt) * sin_lat
      };
}

CoordL ecef_to_local_enu (GeographicCRS * gcrs, Coord3 x, CoordG c0)
{

  Coord3 x0 = geodetic_to_ecef (gcrs, c0),
    dx = (Coord3) { x.x - x0.x, x.y - x0.y, x.z - x0.z};

  double lat0     = _lat (gcrs, c0.lat);
  double lon0     = _lon (gcrs, c0.lon);
  double sin_lat0 = sin (lat0);
  double cos_lat0 = cos (lat0);
  double sin_lon0 = sin (lon0);
  double cos_lon0 = cos (lon0);

  return
    (CoordL)
      {
        .e = -sin_lon0 * dx.x + cos_lon0 * dx.y;
        .n = -sin_lat0 * cos_lon0 * dx.x - sin_lat0 * sin_lon0 * dx.y + cos_lat0 * dx.z;
        .u =  cos_lat0 * cos_lon0 * dx.x + cos_lat0 * sin_lon0 * dx.y + sin_lat0 * dx.z;
      }
}

CoordL geodetic_to_local_enu (GeographicCRS * gcrs, CoordG c, CoordG c0)
{
  Coord3 x = geodetic_to_ecef (gcrs, c);
  return geodetic_to_enu (gcrs, x, c0);
}

Coord3 local_enu_to_ecef (GeographicCRS * gcrs, CoordL l, CoordG c0)
{
  Coord3 x0 = geodetic_to_ecef (gcrs, c0);

  double lat0     = _lat (gcrs, c0.lat);
  double lon0     = _lon (gcrs, c0.lon);
  double sin_lat0 = sin (lat0);
  double cos_lat0 = cos (lat0);
  double sin_lon0 = sin (lon0);
  double cos_lon0 = cos (lon0);

  return
    (Coord3)
      {
        .x = x0.x + -sin_lon0 * l.e
                  +  cos_lon0 * l.n;
        .y = x0.y + -sin_lat0 * cos_lon0 * l.e
                  -  sin_lat0 * sin_lon0 * l.n
                  +  cos_lat0 * l.u;
        .z = x0.z +  cos_lat0 * cos_lon0 * l.e
                  +  cos_lat0 * sin_lon0 * l.n
                  +  sin_lat0 * l.u;
      }
}

/* ECEF -> geodetic, Bowring's method (closed-form, no iteration needed) */
CoordG ecef_to_geodetic (GeographicCRS *gcrs, Coord3 x)
{
  double a       = gcrs->SemiMajorAxis;
  double e2      = gcrs->EccentricitySquared;
  double b       = gcrs->SemiMinorAxis;
  /* second eccentricity squared */
  double ep2     = (a * a) / (b * b) - 1.;
  double p       = sqrt (x.x * x.x + x.y * x.y);
  double theta   = atan2 (x.z * a, p * b);
  double sin_t   = sin (theta);
  double cos_t   = cos (theta);
  double lat     = atan2 ( x.z + ep2 * b * sin_t * sin_t * sin_t,
                           p - e2 * a * cos_t * cos_t * cos_t );
  double lon     = atan2 (x.y, x.x);
  double sin_lat = sin (lat);
  double N       = a / sqrt (1.0 - e2 * sin_lat * sin_lat);
  double h       = p / cos (lat) - N;

  return
    (CoordG)
      {
        .lat = lat / gcrs->AngularUnitSize,
        .lon = lon / gcrs->AngularUnitSize,
        .alt = h
      }
}

CoordG local_enu_to_geodetic (GeographicCRS * gcrs, CoordL l, CoordG c0)
{
  Coord3 x = enu_to_ecef (gcrs, l, c0);
  return ecef_to_geodetic (gcrs, x);
}
