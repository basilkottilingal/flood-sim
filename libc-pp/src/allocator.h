#ifndef _ALLOC_ALLOCATOR_H_
#define _ALLOC_ALLOCATOR_H_

  #ifndef POOL
  # define POOL (1<<20)
  #endif

  #if POOL <  (1<<15)
  # error ("insufficient pool size") 
  #endif

  #include <stdint.h>
  #include <stdlib.h>

  typedef struct
  {
    uint16_t n, max;
  } array;

  /* api */
  void * allocate (size_t);
  char * allocate_str (const char *, int);
  void free_arena ();
  /*
  array * array_insert (array *, void *);
  #define foreach_array_elem(type,a) do {                   \
    type ** ARR = (type **) (((char *)a) + sizeof (array)); \
    for (int i = 0; i<a->n; ++i) {                          \
      type * elem = (type *) ARR [i];
  #define foreach_array_elem_end()                          \
    }                                                       \
  } while (0)
  */

#endif
