#ifndef _ALLOC_ALLOCATOR_H_
#define _ALLOC_ALLOCATOR_H_

  #ifndef _ALLOC_PAGE_SIZE_
  #define _ALLOC_PAGE_SIZE_ (1<<16)
  #endif

  #include <stdlib.h>

  /* suitable for arena allocator */

  typedef struct Allocator
  {
    char * mem;
    struct Allocator * next;
    size_t left;
  } Allocator;

  /* api */
  Allocator * allocator ();
  void        allocator_free (Allocator *);
  void      * allocate  (Allocator *, size_t s);

#endif
