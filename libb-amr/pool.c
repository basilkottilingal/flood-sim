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

/* no MADV_FREE in this system's headers: give memory back at once */
#if !_BAMR_MADV_DONTNEED_RECOMMENDED_ && !defined (MADV_FREE)
# undef  _BAMR_MADV_DONTNEED_RECOMMENDED_
# define _BAMR_MADV_DONTNEED_RECOMMENDED_  1
#endif

/* a list of freed objects, linked through the first word of each */
typedef struct
{
  void * top;
  size_t n;
} FreeList;

struct _Pool
{
  size_t size;          /* bytes per object, rounded                  */
  size_t chunk_bytes;   /* bytes per chunk, a multiple of the page    */
  int    large;         /* objects are whole pages                    */
  FreeList hot, cold;   /* see pool.h                                 */
  FreeList muzzy;       /* with MADV_FREE only: see pool.h            */
  size_t hot_max;       /* large objects the hot list may hold        */
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

static void push (FreeList * l, void * o)
{
  * (void **) o = l->top;
  l->top = o;
  l->n++;
}

static void * pop (FreeList * l)
{
  void * o = l->top;
  l->top = * (void **) o;
  l->n--;
  return o;
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
  p->hot_max = POOL_HOT_ALL;
  return p;
}

static void * fresh (Pool * p)
{
  if (p->next == p->end)
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
  void * o = p->next;           /* fresh pages read as zero */
  p->next += p->size;
  return o;
}

void * pool_alloc (Pool * p)
{
  void * o;
  if (p->cold.n)
  {
    o = pop (&p->cold);
    * (void **) o = NULL;       /* the rest is zero since madvise () */
  }
  else if (p->hot.n)
  {
    o = pop (&p->hot);
    memset (o, 0, p->size);
  }
#if !_BAMR_MADV_DONTNEED_RECOMMENDED_
  else if (p->muzzy.n)
  {
    /* each page old or zero, which cannot be told apart cheaply */
    o = pop (&p->muzzy);
    memset (o, 0, p->size);
  }
#endif
  else
    o = fresh (p);
  p->used++;
  return o;
}

void * pool_alloc_uninit (Pool * p)
{
#if !_BAMR_MADV_DONTNEED_RECOMMENDED_
  /* muzzy before cold: no fault, unless the system took its pages */
  void * o = p->hot.n ? pop (&p->hot) : p->muzzy.n ? pop (&p->muzzy) :
    p->cold.n ? pop (&p->cold) : fresh (p);
#else
  void * o = p->hot.n ? pop (&p->hot) : p->cold.n ? pop (&p->cold) : fresh (p);
#endif
  p->used++;
  return o;
}

/* a large object, its pages given back: it reads as zero from now on */
static void make_cold (Pool * p, void * o)
{
  if (madvise (o, p->size, MADV_DONTNEED))
    memset (o, 0, p->size);     /* still zero, if not given back */
  push (&p->cold, o);
}

#if !_BAMR_MADV_DONTNEED_RECOMMENDED_
/*
.. A large object, its pages given back lazily: all but the first, which
.. holds the list's link, so that writing the link costs no fault. An
.. object of one page has nothing else to give. A kernel without
.. MADV_FREE (before Linux 4.5) refuses it: the object goes cold.
*/
static void make_muzzy (Pool * p, void * o)
{
  size_t page = page_size ();
  if (p->size > page &&
      madvise ((char *) o + page, p->size - page, MADV_FREE))
  {
    make_cold (p, o);
    return;
  }
  push (&p->muzzy, o);
}

/* a large object over the hot limit */
static void release (Pool * p, void * o)
{
  make_muzzy (p, o);
}
#else
/* a large object over the hot limit */
static void release (Pool * p, void * o)
{
  make_cold (p, o);
}
#endif

void pool_free (Pool * p, void * o)
{
  if (!o)
    return;
  if (p->large && p->hot.n >= p->hot_max)
    release (p, o);
  else
    push (&p->hot, o);
  p->used--;
}

void pool_set_hot_max (Pool * p, size_t n)
{
  p->hot_max = n;
  if (!p->large)
    return;                     /* small objects share pages */
  while (p->hot.n > n)
    release (p, pop (&p->hot));
}

void pool_trim (Pool * p, size_t keep)
{
  if (!p->large)
    return;                     /* small objects share pages */
  while (p->hot.n > keep)
    make_cold (p, pop (&p->hot));
#if !_BAMR_MADV_DONTNEED_RECOMMENDED_
  /* lazily given back is not given back: really, and out of the RSS */
  while (p->muzzy.n)
    make_cold (p, pop (&p->muzzy));
#endif
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

size_t pool_hot (const Pool * p)
{
  return p ? p->hot.n : 0;
}

size_t pool_cold (const Pool * p)
{
  return p ? p->cold.n : 0;
}

size_t pool_muzzy (const Pool * p)
{
  return p ? p->muzzy.n : 0;
}

size_t pool_size (const Pool * p)
{
  return p ? p->size : 0;
}
