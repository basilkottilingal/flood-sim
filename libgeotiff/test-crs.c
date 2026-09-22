#include "coordinates.h"
#include "geodetic.h"

#include <stdio.h>

int main ()
{
  GeographicCRS gcrs = EPSG4326_CRS;

  CoordG Paris = (CoordG) {2.3522, 48.8566, 0.};
  Coord3 ecef = geodetic_to_ecef (&gcrs, Paris);
  printf ("Paris {%gE %gN %g}\n", Paris.lon, Paris.lat, Paris.alt);
  printf ("Paris ECEF {%g %g %g}\n", ecef.x, ecef.y, ecef.z);
  CoordG Paris2 = ecef_to_geodetic (&gcrs, ecef);
  printf ("Paris {%gE %gN %g}\n", Paris2.lon, Paris2.lat, Paris2.alt);
  CoordL local = (CoordL) {.e = 2000, .n = 0, .u = 0}; //2 km east of Paris;
  CoordG Off = local_enu_to_geodetic (&gcrs, local, Paris);
  printf ("2Km East of Paris {%gE %gN %g}\n", Off.lon, Off.lat, Off.alt);
  local = (CoordL) {.e = -2000, .n = 0, .u = 0};
  Off.alt *= -1.;
  CoordG Paris3 = local_enu_to_geodetic (&gcrs, local, Off);
  printf ("Paris {%gE %gN %g}\n", Paris3.lon, Paris3.lat, Paris3.alt);
}
