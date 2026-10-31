#ifndef _ALLOC_REUSE_ALLOCATOR_H_
#define _ALLOC_REUSE_ALLOCATOR_H_

  #include <stdlib.h>

  typedef struct Freemem
  {
    struct Freemem * next;
  } Freemem;

  typedef struct RAllocator
  {
    struct RAllocator * next;
    Freemem * reuse;
    char * mem;
    size_t size;
    size_t left;
  } RAllocator;

  /* api */
  RAllocator * r_allocator (const size_t s);
  void       * r_allocate  ();
  void         r_allocate_free (RAllocator *, void *);
  void         r_allocator_free (RAllocator *);

#endif
