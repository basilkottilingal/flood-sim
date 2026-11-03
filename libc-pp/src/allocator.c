#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

#include "allocator.h"

static struct allocator
{
  char * mem;
  struct allocator * next;
  size_t left;
} * alloc = NULL;

#define mask (size_t) (sizeof (void *) - 1)

void allocator_new ()
{
  char * m = malloc (POOL);
  if (!m)
  {
    perror ("malloc pool failed");
    exit (-1);
  }
  memset (m, 0, POOL);
  struct allocator * a = (struct allocator *) m;
  size_t salloc = (sizeof (struct allocator) + mask) & ~mask;
  a->mem  = m + salloc;
  a->next = alloc;
  a->left = (size_t) POOL - salloc;
  alloc = a;
}

void free_arena ()
{
  while (alloc)
  {
    struct allocator * next = alloc->next;
    free (alloc);
    alloc = next; 
  }
}

void * allocate (size_t s)
{
  s = (s + mask) & ~mask;
  if (!alloc || alloc->left < s)
  {
    if (s > POOL/2) {
      perror ("non practical memory demand");
      exit (-1);
    }
    allocator_new ();
  }
  char * m = alloc->mem;
  alloc->left -= s;
  alloc->mem  += s;
  return m;
}

char * allocate_str (const char * yytext, int yylen)
{
  char * str = allocate (yylen + 1);
  memcpy (str, yytext, yylen);
  str [yylen] = '\0';
  return str;
}
