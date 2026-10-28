#ifndef _ALLOC_ARRAY_H_
#define _ALLOC_ARRAY_H_

  #include <stdlib.h>

  typedef struct
  {
    void * p;
    size_t max, len;
  } Array;

  Array * array_new    ( );
  void    array_free   ( Array * );
  void *  array_append ( Array *, const void *, size_t );

#endif
