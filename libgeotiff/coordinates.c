#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <stdlib.h>
#include <errno.h>
#include <assert.h>
#include <math.h>

#include "coordinates.h"
#include "parse.h"
#include "filemap.h"
#include "tiff.h"

/*
References :
http://geotiff.maptools.org/spec/contents.html
https://docs.ogc.org/is/19-008r4/19-008r4.pdf
http://www.opengis.net/doc/IS/GeoTIFF/1.1
https://www.iogp.org/wp-content/uploads/2019/09/373-07-02.pdf
https://gis.stackexchange.com/questions/120636/math-formula-for-transforming-from-epsg4326-to-epsg3857

(Standard) Coordinate Reference Systems (CRS) used for Geospatial coordinates
(Geotiff Key : GTModelTypeGeoKey = 1024)
  |
  |_ Projected Coordinates, ModelTypeProjected   = 1,
  |   |_ Raster to model space
  |   |   |_ Geotiff Tag : ModelTransformationTag 34264 ( DOUBLE X 16 )
  |   |_ Examples
  |       |_ UTM (Universal Traverse Mercator)
  |       |_ Web Mercator
  |       |_ 
  |        
  |_ Geographic Coordinates, ModelTypeGeographic  = 2,
  |   |_ Raster to model space  (tiepoints + scaling).
  |   |   |_ Geotiff Tag : ModelTiepointTag   33922 ( DOUBLE X 6 X N, N >= 1 )
  |   |   |_ Geotiff Tag : ModelPixelScaleTag 33350 ( DOUBLE X 3 )
  |   |_ Examples
  |       |_ WGS 84
  |           |_EPSG:4326 – 2D coordinate reference system (CRS)
  |           |_EPSG:4979 – 3D CRS
  |           |_EPSG:4978 – geocentric 3D CRS
  |           |_EPSG:7030 – reference ellipsoid
  |           |_EPSG:6326 – horizontal datum
  |    
  |_ Geocentric Coordinates, ModelTypeGeocentric  = 3,
  |_ Vertical Coordinates


Pixel Type
  |_  PixelIsPoint (the data corresponds to the point represented by pixel in the CRS)
  |_  PixelIsArea (the data corresponsd to the center of the [I,I+1]x[J,J+1] patch. More like "cell centered" data)
*/

static CRS map = (CRS) {0};

#define error(e)         filemap_close_all (e)
#define not_implemented(e)                      \
  do {                                          \
    fprintf (stderr, "implementation error\n"); \
    filemap_close_all (e);                      \
  } while (0)

static const char * entry_offset_address (uint32_t loc, uint32_t size)
{
  const char * start = filemap_address (FILEMAP_TIFF);
  assert (start != NULL);
  if (loc + size > filemap_size (FILEMAP_TIFF))
    error ("tiff file : insufficient size");
  return start + loc;
}

struct GTTags
{
  TIFFEntry KeyDirectory;
  double *  dparams;
  char   *  aparams;
  uint32_t  ndparams; 
};

void geotiff_keys (struct GTTags gt, CRS * crs)
{

  uint32_t size = gt.KeyDirectory.count * 2;
  const char * address = entry_offset_address (gt.KeyDirectory.value, size);
  uint16_t
    KeyDirectoryVersion = u16 (& address),
    KeyRevision         = u16 (& address),
    MinorRevision       = u16 (& address),
    NumberOfKeys        = u16 (& address);

  printf ("KeyDirectoryVersion %u, KeyRevision %u, MinorRevision %u, NumberOfKeys %u\n",
    KeyDirectoryVersion, KeyRevision, MinorRevision, NumberOfKeys);

  GTModelType gtmodel = ModelTypeUndefined;
  GTPixelType gtpixel = PixelIsUndefined;
  uint16_t    gcrs    = 0;
  uint16_t    pcrs    = 0;
  uint16_t    vcrs    = 0;

  for (int i=0; i<NumberOfKeys; ++i)
  {
    uint16_t KeyID        = u16 (& address);
    uint16_t TIFFTagLoc   = u16 (& address);
    uint16_t Count        = u16 (& address);
    uint16_t Value_Offset = u16 (& address);
        
    switch (KeyID)
    {
      case GTModelTypeGeoKey :
        gtmodel = Value_Offset;
        printf ("\tGTModel %s\n",
          gtmodel == ModelTypeProjected  ? "projected  coordinate system" :
          gtmodel == ModelTypeGeographic ? "geographic coordinate system" : 
          gtmodel == ModelTypeGeocentric ? "geocentric coordinate system" : 
                                             "unknown  coordinate system");
        if (gtmodel != ModelTypeGeographic)
          not_implemented ("Geotiff CRS mode");
        break;
  
      case GTRasterTypeGeoKey :
        gtpixel = Value_Offset;
        printf ("\tGTPixel %s\n",
          gtpixel == PixelIsArea  ? "PixelIsArea"  :
          gtpixel == PixelIsPoint ? "PixelIsPoint" : "PixelIsUnknown");
        if (gtpixel != PixelIsPoint)
          not_implemented ("Geotiff Pixel Type");
        break;
  
      case GeographicTypeGeoKey :
        printf ("\tGeographicTypeGeoKey %6u\n", Value_Offset);
        if (Value_Offset != EPSG4326)
          not_implemented ("Only EPSG4326 available");
        gcrs = EPSG4326;
        break;
  
      case GeogCitationGeoKey :
        printf ("\tGeogCitationGeoKey");
        assert (gt.aparams != NULL);
        printf (" [pos %u] %s\n", Value_Offset, gt.aparams); 
        break;
  
      case GeogGeodeticDatumGeoKey :
        printf ("\tGeogGeodeticDatumGeoKey %u\n", Value_Offset);
        break;
  
      case GeogPrimeMeridianGeoKey :
        printf ("\tGeogPrimeMeridianGeoKey %u\n", Value_Offset);
        break;
  
      case GeogLinearUnitsGeoKey :
        printf ("\tGeogLinearUnitsGeoKey %s\n",
          Value_Offset == 9001 ? "meter" : "unknown");
        break;
  
      case GeogLinearUnitSizeGeoKey :
        printf ("\tGeogLinearUnitSizeGeoKey\n");
        break;
  
      case GeogAngularUnitsGeoKey :
        printf ("\tGeogAngularUnitsGeoKey %s\n",
          Value_Offset == 9101 ? "radian" :
          Value_Offset == 9102 ? "degree" : "unknown");
        break;
  
      case GeogAngularUnitSizeGeoKey :
        printf ("\tGeogAngularUnitSizeGeoKey\n");
        break;
  
      case GeogEllipsoidGeoKey :
        printf ("\tGeogEllipsoidGeoKey\n");
        break;
  
      case GeogSemiMajorAxisGeoKey :
        printf ("\tGeogSemiMajorAxisGeoKey");
        assert (gt.dparams != NULL && Value_Offset < gt.ndparams);
        printf (" %g\n", gt.dparams [Value_Offset]);
        break;
  
      case GeogSemiMinorAxisGeoKey :
        printf ("\tGeogSemiMinorAxisGeoKey");
        assert (gt.dparams != NULL && Value_Offset < gt.ndparams);
        printf (" %10g\n", gt.dparams [Value_Offset]);
        break;
  
      case GeogInvFlatteningGeoKey :
        printf ("\tGeogInvFlatteningGeoKey");
        assert (gt.dparams != NULL && Value_Offset < gt.ndparams);
        printf (" %10g\n", gt.dparams [Value_Offset]);
        break;
  
      case GeogPrimeMeridianLongGeoKey :
        printf ("\tGeogPrimeMeridianLongGeoKey");
        assert (gt.dparams != NULL && Value_Offset < gt.ndparams);
        printf (" %10g\n", gt.dparams [Value_Offset]);
        break;
        
      default :
        printf ("\tGT key : KeyID %6u, TIFFTagLocation %6u, Count %6u, Value_Offset %6u\n",
          KeyID, TIFFTagLoc, Count, Value_Offset);
        
    }
  }

  if (gcrs == EPSG4326)
  {
    GeographicCRS g = EPSG4326_CRS;
    printf ("\t==============default values of EPSG 4326=============\n");
    printf ("\tSemiMajorAxis %g\n", g.SemiMajorAxis);
    printf ("\tSemiMinorAxis %g\n", g.SemiMinorAxis);
    printf ("\tInvFlattening %g\n", g.InvFlattening);
    printf ("\tPrimeMeridian %g degrees\n", g.PrimeMeridianLong);
    printf ("\tLinearUnits Meters\n");
    printf ("\tAngularUnits Degree\n");
  }
}

/*
.. define the coordinate system and mapping parameters using
.. the array of IFD entries. Any inconsistency/non-implemenation error,
.. will call exit (-1).
*/

CRS
geotiff_tags (TIFFEntry * entries)
{
  struct GTTags gt =
    (struct GTTags)
      {
        .KeyDirectory = (TIFFEntry) {.tag = 0},
        .dparams  = NULL,
        .ndparams = 0,
        .aparams  = NULL
      };
  CRS crs = (CRS) {0};
  int exclusive = 0;
  TIFFEntry * ptr = entries;

  while (ptr->tag)
  {

    TIFFEntry entry = *ptr++;
    uint32_t size = entry.count * tiff_datasize (entry.type),
      val_offset = entry.value;
    const char * address =
      size > 4 ? entry_offset_address (val_offset, size) : NULL;

    switch (entry.tag)
    {
      /* geotiff specific tags*/
      case GeoKeyDirectoryTag:
        assert (entry.type == SHORT && gt.KeyDirectory.tag == 0);
        gt.KeyDirectory = entry;
        break;
      case GeoDoubleParamsTag:
        assert (entry.type == DOUBLE && gt.dparams == NULL);
        gt.dparams = malloc (entry.count * sizeof (double));
        gt.ndparams = entry.count;
        for (uint16_t count = 0; count < entry.count; ++count)
          gt.dparams [count] = d64 (&address);
        break;
      case GeoAsciiParamsTag:
        assert (entry.type == ASCII && gt.aparams == NULL);
        gt.aparams = malloc (entry.count);
        memcpy (gt.aparams, address, entry.count);
        assert (gt.aparams [entry.count-1] == '\0');
        break;
      /*
      .. there are two models for raster->coord mapping and
      .. this combination is exclusive and exhaustive.
      .. (a) pixel scale + tie point
      .. (b) matrix transformation model
      */
      case ModelPixelScaleTag :
        assert (
          entry.count     == 3      &&
          entry.type      == DOUBLE &&
          (exclusive & 1) == 0
        );
        for (int i=0; i<3; ++i)
          crs.scale [i] = d64 (&address);
        exclusive |= 1;
        break;
      case ModelTiepointTag :
        assert (
          entry.count % 6 == 0 &&
          entry.type      == DOUBLE &&
          (exclusive & 2) == 0
        );
        for (int i=0; i<6; ++i)
          crs.tiepoint [i] = d64 (&address);
        exclusive |= 2;
        break;
      case ModelTransformationTag :
        assert (
          entry.count     == 16 &&
          entry.type      == DOUBLE &&
          (exclusive & 4) == 0
        );
        for (int i=0; i<4; ++i)
          for (int j=0; j<4; ++j)
            crs.matrix [i][j] = d64 (&address);
        exclusive |= 4;
        not_implemented ("Projection CRS");
        break;

      /* unexpected */
      default :
    }
  }

  if (! (exclusive == 3 || exclusive == 4) )
    error ("inconsistent mapping model params");

  if (gt.KeyDirectory.tag == 0)
    error ("key directory not found");

  geotiff_keys (gt, &crs);

  if (gt.aparams)
    free (gt.aparams);
  if (gt.dparams)
    free (gt.dparams);

  return crs;
}

#if 0


static void read_geo_key_dir (Image * img, TIFFEntry * entry)
{
  const char * const start = filemap_address (FILEMAP_TIFF);
  assert (start != NULL);
  const char * const end   = start + filemap_size (FILEMAP_TIFF);

  assert (entry->tag == GeoKeyDirectoryTag);
  assert (entry->count % 4 == 0);
  assert (entry->type == SHORT);

  is_available (start, end, entry->value + entry->count * 2);

  const char * array = start + entry->value;
  uint16_t
    KeyDirectoryVersion = u16 (&array),
    KeyRevision         = u16 (&array),
    MinorRevision       = u16 (&array),
    NumberOfKeys        = u16 (&array);
  printf ("KeyDirectoryVersion %u, KeyRevision %u, MinorRevision %u, NumberOfKeys %u\n",
    KeyDirectoryVersion, KeyRevision, MinorRevision, NumberOfKeys);

  for (int i=0; i<NumberOfKeys; ++i)
  {
    geotiff_key ( (GeoKey)
      {
        u16 (&array),
        u16 (&array),
        u16 (&array),
        u16 (&array)
      }
    );

  }
}
        
static
void read_geo_pixel_scale (Image * img, TIFFEntry * entry)
{
  /*
  .. scaling of longitude and latitude per pixel.
  */

  const char * const start = filemap_address (FILEMAP_TIFF);
  assert (start != NULL);
  const char * const end   = start + filemap_size (FILEMAP_TIFF);

  assert (entry->tag == ModelPixelScaleTag);
  assert (entry->count == 3);
  assert (entry->type == DOUBLE);
   
  is_available (start, end, entry->value + entry->count * 8);

  const char * array = start + entry->value;
  printf ("pixel scale\n");
  double scale [3] = {d64 (&array), d64 (&array), d64 (&array)};
  printf ("\t(%g, %g, %g) degrees per pixel\n", scale [0], scale [1], scale [2]);
  memcpy (img->coordMap.scale, scale, sizeof (scale));

}

static
void read_geo_tie_point (Image * img, TIFFEntry * entry)
{

  /*
  .. Tie point(s) is a tuple of 6 double numbers.
  .. (I, J, K) represent which pixel corresponds to reference coordinate's origin
  .. (X, Y, Z) represent longitude (-180 deg W, 180 deg E],
  .. latitude [-90 deg N, 90 deg N], and elevation of the reference pixel.
  .. NOTE :
  .. (a) the reference pixel (I, J, K) can be fractional. 
  .. (b) There can be multiple tie points. 

  */

  const char * const start = filemap_address (FILEMAP_TIFF);
  assert (start != NULL);
  const char * const end   = start + filemap_size (FILEMAP_TIFF);

  assert (entry->tag == ModelTiepointTag);
  assert (entry->count % 6 == 0);
  assert (entry->type == DOUBLE);
   
  is_available (start, end, entry->value + entry->count * 8);

  const char * array = start + entry->value;
  printf ("tie point\n");

  double tiepoint [6] = 
    { 
      d64 (&array), d64 (&array), d64 (&array),
      d64 (&array), d64 (&array), d64 (&array)
    };

  printf ("\ttiepoint pixel [%g, %g, %g]\n", tiepoint [0], tiepoint [1], tiepoint [2]);
  printf ("\ttiepoint coord (%g, %g, %g)\n", tiepoint [3], tiepoint [4], tiepoint [5]);
  printf ("\tInference (longitude %g latitude %g)\n", tiepoint [3], tiepoint [4]);

  if (entry->count > 6)
    error ("warning : library not designed for multiple tie points");

  memcpy (img->coordMap.tiepoint, tiepoint, sizeof (tiepoint));
}

static
void read_geo_ascii_params (Image * img, TIFFEntry * entry)
{

  const char * const start = filemap_address (FILEMAP_TIFF);
  assert (start != NULL);
  const char * const end   = start + filemap_size (FILEMAP_TIFF);

  assert (entry->tag  == GeoAsciiParamsTag);
  assert (entry->type == ASCII);
   
  is_available (start, end, entry->value + entry->count);

  const char * params = start + entry->value;
  printf ("geo ascii params\n\t%s\n", params);

}
#endif

#undef error
