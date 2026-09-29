#include "grid.h"

#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <stdint.h>

#include <stdio.h>

/* fixme : use memory pooling */
void ** grid_general (size_t s, int w, int h, int ng)
{
  /* warning : use this only for uint8_t, int and double */
  assert (w > 0 && h > 0 && ng >= 0);
  size_t size = (w + 2*ng) * (h + 2*ng) * s;

  #define _grid_(type)                                          \
    if (s == sizeof (type))                                     \
    {                                                           \
      char * address = malloc (size + sizeof (GridData));       \
      if (address == NULL)                                      \
        return NULL;                                            \
      *(GridData *)address = (GridData){w,h,ng,s};              \
      address += sizeof (GridData);                             \
      type * mem = (type *) (address);                          \
      memset (mem, 0, size);                                    \
      type ** memptr = malloc ((h + 2*ng) * sizeof (type *));   \
      if (memptr == NULL)                                       \
        return NULL;                                            \
      for (int y = -ng; y < h + ng; ++y)                        \
        memptr [y + ng] = & mem [(y + ng) * (w + 2*ng) + ng];   \
      return (void **) (memptr + ng);                           \
    }

  _grid_ (uint8_t)
  _grid_ (double)
  _grid_ (int)

  #undef _grid_
  assert ("unknown type" && 0);
 
  return NULL;
}
