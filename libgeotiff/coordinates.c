#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <stdlib.h>
#include <errno.h>
#include <assert.h>
#include <math.h>

#include "coordinates.h"
#include "filemap.h"

/*
References :
http://geotiff.maptools.org/spec/contents.html
https://docs.ogc.org/is/19-008r4/19-008r4.pdf
http://www.opengis.net/doc/IS/GeoTIFF/1.1
https://www.iogp.org/wp-content/uploads/2019/09/373-07-02.pdf
https://gis.stackexchange.com/questions/120636/math-formula-for-transforming-from-epsg4326-to-epsg3857

Model Type
  |_______Projected
  |
  |_______Geographic
          |___ WGS 84
  ??            |_ EPSG:4326 – 2D coordinate reference system (CRS)
              |_ EPSG:4979 – 3D CRS
              |_ EPSG:4978 – geocentric 3D CRS
              |_ EPSG:7030 – reference ellipsoid
              |_ EPSG:6326 – horizontal datum

Raster space to model space mapping
  |
  |_______ tiepoints + scaling
  |       |__ ModelTiepointTag 33922 ( DOUBLE X 6 X N, N >= 1 )
  |       | 
  |       |__ ModelPixelScaleTag 33350  ( DOUBLE X 3 )
  |
  |________  
          |
          |__ ModelTransformationTag 34264 ( DOUBLE X 16 )

(Standard) Coordinate systems used for Geospatial coordinates
  |
  |_ Projected Coordinates    ,ModelTypeProjected   = 1,
  |       |
  |       |_ 
  |
  |_ Geographic Coordinates   ,ModelTypeGeographic  = 2,
  |       |
  |       |_
  |
  |_ Geocentric Coordinates   ,ModelTypeGeocentric  = 3,
  |       |
  |       |_
  |
  |_ Vertical Coordinates (I don't think it is not used in geotiff. May use userdefine option.


Pixel Type
  |_  PixelIsPoint (the data corresponds to the point represented by pixel in the Coord Sys)
  |
  |_  PixelIsArea (the data corresponsd to the center of the [I,I+1]x[J,J+1] patch. More like "cell centered" data)
*/

static GeoCoordSys map = (GeoCoordSys) {0};

#define error(e)                              \
  do {                                        \
    fprintf (stderr, "implementation error"); \
    filemap_close_all (e);                    \
  } while (0) 

void geotiff_key (GeoKey key)
{

  switch (key.KeyID)
  {
    case GTModelTypeGeoKey :
      map.model = key.Value_Offset;
      printf ("\tGTModel %s\n",
        map.model == ModelTypeProjected  ? "projected coordinate system" :
        map.model == ModelTypeGeographic ? "geographic coordinate syatem" : 
        map.model == ModelTypeGeocentric ? "geocentric coordinate syatem" : 
                                           "unknown coordinate system");
      break;

    case GTRasterTypeGeoKey :
      map.pixel = key.Value_Offset;
      printf ("\tGTPixel %s\n",
        map.pixel == PixelIsArea  ? "PixelIsArea"  :
        map.pixel == PixelIsPoint ? "PixelIsPoint" : "PixelIsUnknown");
      if (map.pixel != PixelIsPoint)
        error ("unknown pixel type");
      break;

    case GeographicTypeGeoKey :
      map.gCRS = key.Value_Offset;
      printf ("\tGeographicTypeGeoKey %6u\n", map.gCRS);
      if (map.gCRS != EPSG4326)
        error ("unknown EPSG CRS reference");
      break;

    case GeogCitationGeoKey :
      printf ("\tGeogCitationGeoKey\n");
      break;

    case GeogGeodeticDatumGeoKey :
      printf ("\tGeogGeodeticDatumGeoKey\n");
      break;

    case GeogPrimeMeridianGeoKey :
      printf ("\tGeogPrimeMeridianGeoKey\n");
      break;

    case GeogLinearUnitsGeoKey :
      printf ("\tGeogLinearUnitsGeoKey\n");
      break;

    case GeogLinearUnitSizeGeoKey :
      printf ("\tGeogLinearUnitSizeGeoKey\n");
      break;

    case GeogAngularUnitsGeoKey :
      printf ("\tGeogAngularUnitsGeoKey\n");
      break;

    case GeogAngularUnitSizeGeoKey :
      printf ("\tGeogAngularUnitSizeGeoKey\n");
      break;

    case GeogEllipsoidGeoKey :
      printf ("\tGeogEllipsoidGeoKey\n");
      break;

    case GeogSemiMajorAxisGeoKey :
      printf ("\tGeogSemiMajorAxisGeoKey\n");
      break;

    case GeogSemiMinorAxisGeoKey :
      printf ("\tGeogSemiMinorAxisGeoKey\n");
      break;

    case GeogInvFlatteningGeoKey :
      printf ("\tGeogInvFlatteningGeoKey\n");
      break;

    case GeogPrimeMeridianLongGeoKey :
      printf ("\tGeogPrimeMeridianLongGeoKey\n");
      break;
      
    default :
      printf ("\tGT key : KeyID %6u, TIFFTagLocation %6u, Count %6u, Value_Offset %6u\n",
        key.KeyID, key.TIFFTagLocation, key.Count, key.Value_Offset);
      
  }
}

#undef error
