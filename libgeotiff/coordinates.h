#ifndef _GEOTIFF_COORDINATES_H_
#define _GEOTIFF_COORDINATES_H_


  typedef enum
  {
    GeoKeyDirectoryTag            = 34735,
    GeoDoubleParamsTag            = 34736,
    GeoAsciiParamsTag             = 34737,
    ModelPixelScaleTag            = 33550,
    ModelTiepointTag              = 33922,
    ModelTransformationTag        = 34264,
  } GeotiffTags;

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
  } GTKeys;

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

    EPSG4326                = 4326,

    /*
    .. note : every other Geographic CRS not implemented.
    */
  } GTGeographicCRS;
  
  typedef enum 
  {
    ModelTypeUndefined      = 0,
    ModelTypeProjected      = 1,
    ModelTypeGeographic     = 2,
    ModelTypeGeocentric     = 3,
    ModelTypeUserDefined    = 32767
  } GTModel;
  
  typedef enum 
  {
    PixelIsUndefined        = 0,
    PixelIsArea             = 1,
    PixelIsPoint            = 2,
  } GTPixel;

  typedef struct
  {
    GTModel         model;
    GTPixel         pixel;
    GTGeographicCRS gCRS; /* only applies to Geog CRS */
  
    double matrix   [4][4];
    double scale    [3];
    double tiepoint [6];
   
  } GeoCoordSys;

  typedef struct
  {
    uint16_t KeyID;
    uint16_t TIFFTagLocation;
    uint16_t Count;
    uint16_t Value_Offset;
  } GeoKey;

  /* api */
  void geotiff_key (GeoKey key);

#endif
