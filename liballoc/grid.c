#include "grid.h"

#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <stdint.h>

void ** grid_general (size_t s, int w, int h, int ng)
{
  /* warning : use this only for uint8_t, float and double */
  assert (w > 0 && h > 0 && ng >= 0);
  size_t size = (w + 2*ng) * (h + 2*ng) * s;

  #define _grid_(type)                                          \
    if (s == sizeof (type))                                     \
    {                                                           \
      type * mem = malloc (size);                               \
      if (mem == NULL)                                          \
        return NULL;                                            \
      memset (mem, 0, size);                                    \
      type ** memptr = malloc ((h + 2*ng) * sizeof (type *));   \
      if (memptr == NULL)                                       \
        return NULL;                                            \
      for (int y = -ng; y < h + ng; ++y)                        \
        memptr [y + ng] = & mem [(y + ng) * (w + 2*ng) + ng];   \
      return (void **) (memptr + ng);                           \
    }

  _grid_ (uint8_t)
  _grid_ (float)
  _grid_ (double)

  #undef _grid_
  
  return NULL;
}
