#include <assert.h>
#include <string.h>

#include "array.h"

Array * array_new ()
{
  Array * a = malloc (sizeof (Array));
  a->p = NULL;
  a->max = a->len = 0;
  return a;
}

void array_free (Array * a)
{
  free (a->p);
  free (a);
}

void * array_append (Array * a, const void * elem, size_t size)
{
  if (a->len + size >= a->max)
  {
    size_t mask = (1 << 12) - 1;
    a->max += (size + mask) & ~mask;
    a->p = realloc (a->p, a->max);
    assert (a->p);
  }
  char * dest = ((char *)a->p) + a->len;
  memcpy (dest, elem, size);
  a->len += size;
  return (void *) dest;
}
