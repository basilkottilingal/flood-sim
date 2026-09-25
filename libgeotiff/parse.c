#include <assert.h>
#include <stdlib.h>
#include <string.h>

#include "parse.h"

ENDIANNESS endianness = LITTLE_ENDIAN;

void parser_endianness (ENDIANNESS e)
{
  assert (e == LITTLE_ENDIAN || e == BIG_ENDIAN);
  endianness = e;
}

uint16_t u16 (const char ** m)
{
  uint8_t * b = (uint8_t *) *m;
  *m += 2;
  return endianness ==  LITTLE_ENDIAN ?
    (((uint16_t) b[1] << 8) | (uint16_t) b[0]) :
    (((uint16_t) b[0] << 8) | (uint16_t) b[1]);
}

uint32_t u32 (const char ** m)
{
  uint8_t * b = (uint8_t *) *m;
  *m += 4;
  return endianness ==  LITTLE_ENDIAN ?
    ( ((uint32_t) b[3] << 24) |
      ((uint32_t) b[2] << 16) |
      ((uint32_t) b[1] << 8 ) |
       (uint32_t) b[0]      )
    :
    ( ((uint32_t) b[0] << 24) |
      ((uint32_t) b[1] << 16) |
      ((uint32_t) b[2] << 8 ) |
       (uint32_t) b[3]      );
}

float f32 (const char ** m)
{
  //some_static_assert (sizeof(float) == sizeof (uint32_t));

  uint8_t * b = (uint8_t *) *m;
  *m += 4;
  uint32_t val = endianness ==  LITTLE_ENDIAN ?
    ( ((uint32_t) b[3] << 24) |
      ((uint32_t) b[2] << 16) |
      ((uint32_t) b[1] << 8 ) |
       (uint32_t) b[0]      )
    :
    ( ((uint32_t) b[0] << 24) |
      ((uint32_t) b[1] << 16) |
      ((uint32_t) b[2] << 8 ) |
       (uint32_t) b[3]      );
  float f;
  memcpy (&f, &val, sizeof(float));
  return f;
}

double d64 (const char ** m)
{
  const uint8_t *b = (const uint8_t *) *m;
  uint64_t r;

  *m += 8;

  if (endianness == LITTLE_ENDIAN)
    r =
      ((uint64_t)b[7] << 56) |
      ((uint64_t)b[6] << 48) |
      ((uint64_t)b[5] << 40) |
      ((uint64_t)b[4] << 32) |
      ((uint64_t)b[3] << 24) |
      ((uint64_t)b[2] << 16) |
      ((uint64_t)b[1] <<  8) |
      ((uint64_t)b[0]);
  else
    r =
      ((uint64_t)b[0] << 56) |
      ((uint64_t)b[1] << 48) |
      ((uint64_t)b[2] << 40) |
      ((uint64_t)b[3] << 32) |
      ((uint64_t)b[4] << 24) |
      ((uint64_t)b[5] << 16) |
      ((uint64_t)b[6] <<  8) |
      ((uint64_t)b[7]);

  double val;
  memcpy (&val, &r, sizeof val);
  return val;
}
