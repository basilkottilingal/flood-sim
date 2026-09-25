#ifndef _GEOTIFF_PARSE_H_
#define _GEOTIFF_PARSE_H_

  #include <inttypes.h>

  typedef enum
  {
    BIG_ENDIAN    = 0,
    LITTLE_ENDIAN = 1,
  } ENDIANNESS;

  void     parser_endianness (ENDIANNESS);
  uint16_t u16 (const char ** m);
  uint32_t u32 (const char ** m);
  float    f32 (const char ** m);
  double   d64 (const char ** m);

#endif
