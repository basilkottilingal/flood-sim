#include "stencil.h"
#include <cwrap/event.h>
scalar s;
@event step (i++)
{
  foreach ()
    if (s [] < 0)
      return 3;
  return 0;
}
