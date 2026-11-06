#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <unistd.h>
#include <sys/mman.h>
#ifndef MAP_ANONYMOUS
#define MAP_ANONYMOUS MAP_ANON
#endif

#include "pool.h"

struct _Pool
{
  size_t size;          /* bytes per object, rounded                  */
  size_t chunk_bytes;   /* bytes per chunk, a multiple of the page    */
  int    large;         /* objects are whole pages, given back on free */
  void * free;          /* free list, linked through each object      */
  char * next, * end;   /* the unused part of the newest chunk        */
  void ** chunks;       /* every chunk, for pool_destroy ()           */
  size_t nchunks, maxchunks;
  size_t used;
};

static size_t round_up (size_t n, size_t a)
{
  return (n + a - 1) / a * a;
}

static size_t page_size (void)
{
  static size_t page = 0;
  if (!page)
  {
    long p = sysconf (_SC_PAGESIZE);
    page = p > 0 ? (size_t) p : 4096;
  }
  return page;
}

static void out_of_memory (const char * what, size_t bytes)
{
  fprintf (stderr, "pool: %s of %zu bytes failed\n", what, bytes);
  abort ();
}

Pool * pool_new (size_t size)
{
  Pool * p = calloc (1, sizeof (Pool));
  if (!p)
    out_of_memory ("calloc", sizeof (Pool));
  if (size < sizeof (void *))
    size = sizeof (void *);
  size_t page = page_size ();
  if (size >= POOL_LARGE)
  {
    p->large = 1;
    p->size = round_up (size, page);
  }
  else
    p->size = round_up (size, size >= 64 ? 64 : 16);
  size_t n = POOL_CHUNK / p->size;
  p->chunk_bytes = round_up ((n ? n : 1) * p->size, page);
  return p;
}

static void new_chunk (Pool * p)
{
  char * c = mmap (NULL, p->chunk_bytes, PROT_READ | PROT_WRITE,
                   MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
  if (c == MAP_FAILED)
    out_of_memory ("mmap", p->chunk_bytes);
  if (p->nchunks == p->maxchunks)
  {
    size_t m = p->maxchunks ? 2 * p->maxchunks : 16;
    void ** a = realloc (p->chunks, m * sizeof (void *));
    if (!a)
      out_of_memory ("realloc", m * sizeof (void *));
    p->chunks = a;
    p->maxchunks = m;
  }
  p->chunks [p->nchunks++] = c;
  p->next = c;
  /* the tail past the last whole object stays unused */
  p->end = c + p->chunk_bytes / p->size * p->size;
}

void * pool_alloc (Pool * p)
{
  void * o;
  if (p->free)
  {
    o = p->free;
    p->free = * (void **) o;
    if (p->large)
      * (void **) o = NULL;     /* the rest was zeroed by madvise () */
    else
      memset (o, 0, p->size);
  }
  else
  {
    if (p->next == p->end)
      new_chunk (p);
    o = p->next;                /* fresh pages read as zero */
    p->next += p->size;
  }
  p->used++;
  return o;
}

void pool_free (Pool * p, void * o)
{
  if (!o)
    return;
  /*
  .. A large object starts on a page and covers whole pages: drop them, so
  .. they read as zero and use no memory until touched again.
  */
  if (p->large && madvise (o, p->size, MADV_DONTNEED))
    memset (o, 0, p->size);
  * (void **) o = p->free;
  p->free = o;
  p->used--;
}

void pool_destroy (Pool * p)
{
  if (!p)
    return;
  for (size_t i = 0; i < p->nchunks; i++)
    munmap (p->chunks [i], p->chunk_bytes);
  free (p->chunks);
  free (p);
}

size_t pool_used (const Pool * p)
{
  return p ? p->used : 0;
}

size_t pool_size (const Pool * p)
{
  return p ? p->size : 0;
}
