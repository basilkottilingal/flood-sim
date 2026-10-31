#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

#include "reuse-allocator.h"

RAllocator * r_allocator (size_t s)
{
  size_t sptr  = sizeof (void *);
  size_t mask  = sptr - 1;
  size_t MB    = 1 << 20;
  size_t headsize = (sizeof (RAllocator) + mask) & ~mask;

  s = (s + mask)  & ~mask;
  if (s > MB/2)
  {
    fprintf (stderr, "\n%zd > .5 MB", s);
    fflush (stderr);
  }
  assert (s <= MB/2 ); /* cannot meet huge demands > .5MB */
  size_t blocksize = s * 1024 + headsize;
  if (blocksize > MB)
    blocksize = MB;

  char * m = malloc (blocksize);
  assert ("reuse_allocator () : malloc ()" && m);
  RAllocator * ra = (RAllocator *) m;

  ra->next = NULL;
  ra->size = s;
  ra->left = blocksize - headsize;
  ra->mem  = m + headsize;
  ra->reuse = NULL;

  memset (ra->mem, 0, ra->left);
  return ra;
}

void r_allocator_free (RAllocator * ra)
{
  while (ra)
  {
    RAllocator * next = ra->next;
    free (ra);
    ra = next; 
  }
}

void * r_allocate (RAllocator * ra)
{
  size_t size = ra->size;
  if (ra->left < size)
  {
    /* reuse available? */
    Freemem *m = ra->reuse;
    if (m)
    {
      ra->reuse = ra->reuse->next;
      memset (m, 0, size);
      return m;
    }
    /* createnew block */
    RAllocator * newra = r_allocator (size);
    newra->next = ra->next;
    ra->next = newra;
    ra->mem  = newra->mem;
    ra->left = newra->left;
  }
  char * m = ra->mem;
  ra->left -= size;
  ra->mem  += size;
  return m;
}

void r_allocate_free (RAllocator * ra, void * m)
{
  Freemem * reuse = (Freemem *)m;
  reuse->next = ra->reuse;
  ra->reuse = reuse;
}
