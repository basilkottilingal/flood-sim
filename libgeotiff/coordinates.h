#ifndef _GEOTIFF_COORDINATES_H_
#define _GEOTIFF_COORDINATES_H_

  #include <math.h>
  #include "tiff.h"

  typedef enum
  {
    GeoKeyDirectoryTag            = 34735,
    GeoDoubleParamsTag            = 34736,
    GeoAsciiParamsTag             = 34737,
    ModelPixelScaleTag            = 33550,
    ModelTiepointTag              = 33922,
    ModelTransformationTag        = 34264,
  } GTTagsType;

  typedef enum
  {
    GTModelTypeGeoKey             = 1024,
    GTRasterTypeGeoKey            = 1025,
 
    /* related to Geographic CRS */ 
    GeographicTypeGeoKey          = 2048,
    GeogCitationGeoKey            = 2049,
    GeogGeodeticDatumGeoKey       = 2050,

    GeogPrimeMeridianGeoKey       = 2051,
    GeogLinearUnitsGeoKey         = 2052,
    GeogLinearUnitSizeGeoKey      = 2053,
    GeogAngularUnitsGeoKey        = 2054,
    GeogAngularUnitSizeGeoKey     = 2055,
    GeogEllipsoidGeoKey           = 2056,
    GeogSemiMajorAxisGeoKey       = 2057,
    GeogSemiMinorAxisGeoKey       = 2058,
    GeogInvFlatteningGeoKey       = 2059,

    GeogPrimeMeridianLongGeoKey   = 2061, 

    /* Projected CRS */ 
    GeogAzimuthUnitsGeoKey        = 2060,

    ProjectedCSTypeGeoKey         = 3072,
    PCSCitationGeoKey             = 3073,
    ProjectionGeoKey              = 3074,
    ProjCoordTransGeoKey          = 3075,
    ProjLinearUnitsGeoKey         = 3076,
    ProjLinearUnitSizeGeoKey      = 3077,
    ProjStdParallel1GeoKey        = 3078,
    ProjStdParallel2GeoKey        = 3079,
    ProjNatOriginLongGeoKey       = 3080,

    ProjNatOriginLatGeoKey        = 3081,
    ProjFalseEastingGeoKey        = 3082,
    ProjFalseNorthingGeoKey       = 3083,
    ProjFalseOriginLongGeoKey     = 3084,
    ProjFalseOriginLatGeoKey      = 3085,
    ProjFalseOriginEastingGeoKey  = 3086,
    ProjFalseOriginNorthingGeoKey = 3087,
    ProjCenterLongGeoKey          = 3088,
    ProjCenterLatGeoKey           = 3089,
    ProjCenterEastingGeoKey       = 3090,

    ProjCenterNorthingGeoKey      = 3091,
    ProjScaleAtNatOriginGeoKey    = 3092,
    ProjScaleAtCenterGeoKey       = 3093,
    ProjAzimuthAngleGeoKey        = 3094,
    ProjStraightVertPoleLongGeoKey = 3095,

    /* vertical CRS */
    VerticalCSTypeGeoKey          = 4096,
    VerticalCitationGeoKey        = 4097,
    VerticalDatumGeoKey           = 4098,
    VerticalUnitsGeoKey           = 4099
  } GTKeysType;

  typedef enum
  {
    /*
    .. 0              = undefined
    .. [    1,  1000] = Obsolete EPSG/POSC Geographic Codes
    .. [ 1001,  3999] = Reserved by GeoTIFF
    .. [ 4000,  4199] = EPSG GCS Based on Ellipsoid only
    .. [ 4200,  4999] = EPSG GCS Based on EPSG Datum
    .. [ 5000, 32766] = Reserved by GeoTIFF
    .. 32767          = user-defined GCS
    .. [32768, 65535] = Private User Implementations
    */
    EPSG4326                      = 4326,
    /*
    .. note : every other Geographic CRS not implemented.
    */
  } GTGeographicCRS;
  
  typedef enum 
  {
    ModelTypeUndefined            = 0,
    ModelTypeProjected            = 1,
    ModelTypeGeographic           = 2,
    ModelTypeGeocentric           = 3,
    ModelTypeUserDefined          = 32767
  } GTModelType;
  
  typedef enum 
  {
    PixelIsUndefined              = 0,
    PixelIsArea                   = 1,
    PixelIsPoint                  = 2,
  } GTPixelType;

  /*
  .. Coordinates used in
  .. (a) Geographic Coordinate Reference System
  .. (b) Geocentric/Cartesian/ECEF Coordinates System.
  .. (c) Local East-North-Up Coordinates;
  .. Units expected are respectively
  .. (a) (degrees, degrees, meter),
  .. (b) (meter, meter, meter)
  .. (c) (meter, meter, meter)
  */
  typedef struct
  {
    double lon, lat, alt;
  } CoordG;

  typedef struct
  {
    double x, y, z;
  } Coord3;

  typedef struct
  {
    double e, n, u;
  } CoordL;

  typedef struct
  {
    uint16_t ModelType;
    /*
    .. The following are redundant information if ModelType is some
    .. standard models like EPSG 4326 (i.e WSG 84 - 2D ), in which
    .. case radii, flattening, reference datum, ellipsoid model,
    .. are already known.
    */
    uint16_t PrimeMeridian;
    uint16_t AngularUnits;
    uint16_t LinearUnits;
    uint16_t EllipsoidReference;
    uint16_t Dimension;

    double   PrimeMeridianLong; 
    double   AngularUnitSize;
    double   LinearUnitsSize;
    double   SemiMajorAxis;
    double   SemiMinorAxis;
    double   InvFlattening;
    double   EccentricitySquared;
  } GeographicCRS;

  #if 0
  typedef struct
  {
    uint16_t ModelType;
  } ProjectedCRS;

  typedef struct
  {
    uint16_t ModelType;
  } VerticalCRS;
  #endif

  #define EPSG4326_CRS                                     \
    (GeographicCRS)                                        \
      {                                                    \
        .ModelType            = 4326,                      \
        .Dimension            = 2,                         \
        .PrimeMeridian        = 8901,                      \
        .PrimeMeridianLong    = 0.0,                       \
        .AngularUnits         = 9102,                      \
        .AngularUnitSize      = M_PI/180.,                 \
        .LinearUnits          = 9001,                      \
        .LinearUnitsSize      = 1.0,                       \
        .EllipsoidReference   = 7030,                      \
        .SemiMajorAxis        = 6378137.0,                 \
        .SemiMinorAxis        = 6356752.3142,              \
        .InvFlattening        = 298.257223563,             \
        .EccentricitySquared  = 6.69437999014E-3           \
      }

  typedef struct
  {
    GTModelType     model;
    GTPixelType     pixel;

    /* other CRS models are not implemented */
    GeographicCRS   gcrs;

    double matrix   [4][4];
    double scale    [3];
    double tiepoint [6];
   
  } CRS;

  /*
  .. APIs
  .. 1. geotiff_tags ()    : reads all geotiff tags and create a Coordinate Ref.
  .. 2. coordinate_map ()  : maps a geodetic coordinate (long, lat) to the
  ..    index (i, j) of the raster.
  */
  CRS     geotiff_tags          (TIFFEntry  * entries);
  void    coordinate_map        (CRS * crs, CoordG c, double * index);
  void    tiepoint_new          (CRS * crs, CoordG c, double * new_tiepoint);

#endif
