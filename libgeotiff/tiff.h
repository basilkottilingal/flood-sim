#ifndef _GEOTIFF_TIFF_H_
#define _GEOTIFF_TIFF_H_

  #include <inttypes.h>
  #include <stdlib.h>

  typedef struct
  {
    uint16_t tag;
    uint16_t type; 
    uint32_t count;
    uint32_t value;
  } TIFFEntry;

  typedef struct
  {
    uint16_t nsamples;
    uint16_t bits;
    uint16_t format;
    uint16_t planar;
    void     ** data;
  } TIFFRaster;

  typedef struct
  {
    uint16_t pred;
    uint16_t type;
  } TIFFCompression;

  typedef enum
  {
    T_UINT = 1,
    T_INT,
    T_FLOAT,
    T_UNDEF
  } TIFFFormat;

  typedef enum
  {
    BYTE = 1,
    ASCII,
    SHORT,
    LONG,
    RATIONAL,
    SBYTE,
    UNDEFINED,
    SSHORT,
    SLONG,
    SRATIONAL,
    FLOAT,
    DOUBLE
  } TIFFType;
  
  typedef enum
  {
    TIFF_TAG_IMAGE_WIDTH        = 256,
    TIFF_TAG_IMAGE_LENGTH       = 257,
    TIFF_TAG_BITS_PER_SAMPLE    = 258,
    TIFF_TAG_COMPRESSION        = 259,
    TIFF_TAG_PHOTOMETRIC_INTERP = 262,
  
    TIFF_TAG_STRIP_OFFSETS      = 273,
    TIFF_TAG_SAMPLES_PER_PIXEL  = 277,
    TIFF_TAG_ROWS_PER_STRIP     = 278,
    TIFF_TAG_STRIP_BYTE_COUNTS  = 279,
  
    TIFF_TAG_PLANAR_CONFIG      = 284,
    TIFF_TAG_PREDICTOR          = 317,
  
    TIFF_TAG_TILE_WIDTH         = 322,
    TIFF_TAG_TILE_LENGTH        = 323,
    TIFF_TAG_TILE_OFFSETS       = 324,
    TIFF_TAG_TILE_BYTE_COUNTS   = 325,
  
    TIFF_TAG_SAMPLE_FORMAT      = 339,

    /* this is not a TIFF tag */
    TIFF_TAG_NOT_A_TAG          = 0,

  } TIFFTag;
  
  typedef enum
  {
  
    TIFF_COMP_NONE              = 1,
    TIFF_COMP_LZW               = 5,
    TIFF_COMP_PACKBITS          = 32773,
    TIFF_COMP_DEFLATE           = 8
  } TIFFCompressionType;

  typedef enum
  {
    TIFF_PREDICTOR_DEFAULT      = 1,
    TIFF_PREDICTOR_HORIZONTAL_DIFFERENCING
  } TIFFCompressionPredictor;

  /* api */ 
  uint32_t tiff_datasize (TIFFType type);
  
#endif 
