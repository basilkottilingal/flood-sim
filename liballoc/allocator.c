#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "allocator.h"

Allocator * allocator ()
{
  char * m = malloc (_ALLOC_PAGE_SIZE_);
  assert ("allocator () : malloc ()" && m);
  memset (m, 0, _ALLOC_PAGE_SIZE_);
  Allocator * alloc = (Allocator *) m;
    
  size_t sptr = sizeof (void *);
  size_t mask = sptr - 1;
  size_t salloc = (sizeof (Allocator) + mask) & ~mask;
  *alloc = (Allocator) 
    {
      .mem  = m + salloc,
      .next = NULL,
      .left = (size_t) _ALLOC_PAGE_SIZE_ - salloc
    };
  return alloc;
}

void allocator_free (Allocator * alloc)
{
  while (alloc)
  {
    Allocator * next = alloc->next;
    free (alloc);
    alloc = next; 
  }
}

void * allocate (Allocator * alloc, size_t s)
{
  size_t sptr = sizeof (void *), mask = sptr - 1;
  s = (s + mask) & ~mask;
  if (alloc->left < s)
  {
    assert ("allocate () : heavy demand" && s <= _ALLOC_PAGE_SIZE_ / 2);
    Allocator * newalloc = allocator ();
    newalloc->next = alloc->next;
    alloc->next = newalloc;
    alloc->mem  = newalloc->mem;
    alloc->left = newalloc->left;
  }
  char * m = alloc->mem;
  alloc->left -= s;
  alloc->mem  += s;
  return m;
}
