#include <assert.h>

#include "tiff.h"

uint32_t tiff_datasize (TIFFType type)
{
  switch ( type )
  {
    case BYTE     :
    case SBYTE    :
    case ASCII    :
      return 1u;

    case SHORT    :
    case SSHORT   :
      return 2u;

    case FLOAT    :
    case LONG     :
    case SLONG    :
      return 4u;

    case DOUBLE   :
    case SRATIONAL:
    case RATIONAL :
      return 8u;

    case UNDEFINED:
      return 0u;

    default       :
  }
  assert (0);
  return 0u;
}
