#include <assert.h>

#include "tiff.h"

size_t tiff_datasize (TIFFType type)
{
  switch ( type )
  {
    case BYTE     :
    case SBYTE    :
    case ASCII    :
      return 1;

    case SHORT    :
    case SSHORT   :
      return 2;

    case FLOAT    :
    case LONG     :
    case SLONG    :
      return 4;

    case DOUBLE   :
    case SRATIONAL:
    case RATIONAL :
      return 8;

    case UNDEFINED:
      return 0;

    default       :
  }
  assert (0);
  return 0;
}
