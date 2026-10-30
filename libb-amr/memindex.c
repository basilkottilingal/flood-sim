/*
 * Portions of this file are derived from basilik <https://basilisk.fr/>
 * Copyright (C) [2019-2026] [Stephane Popinet et al]
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 *
 * Modifications Copyright (C) [2026] [Ahmed Basil KOTTILINGAL]
 */

/**
# Multi-dimensional arrays using virtual memory

This implementation uses the Unix/Linux virtual memory subsystem to
more efficiently store sparse multidimensional arrays.

The `Memindex` structure defines multi-dimensional arrays. The
interface is that used by [/src/grid/tree.h](). */

#include <stdlib.h>
#include <assert.h>
#include <stdio.h>

#include <sys/mman.h>
#ifndef MAP_ANONYMOUS
#define MAP_ANONYMOUS MAP_ANON
#endif

#include "memindex.h"

/**
The `mem_data()` macros return the data stored at a specific
(multidimensional) index. It assumes that the index is allocated. This
can be checked with `mem_allocated()`. */

#if dimension == 1
extern inline
int mem_allocated (Memindex m, int i);
#elif dimension == 2
extern inline
int mem_allocated (Memindex m, int i, int j);
#else // dimension == 3
extern inline
int mem_allocated (Memindex m, int i, int j, int k);
#endif // dimension == 3

/**
The `mem_new()` function returns a new (empty) `Memindex`. */

Memindex mem_new (int len)
{
  (void) len; //??
  Memindex m = calloc (1, sizeof (struct _Memindex));
  return m;
}

/**
The `mem_destroy()` function frees all the memory allocated by a given
`Memindex`. */

void mem_destroy (Memindex m, int len)
{
#if dimension > 1
  for (int i = 0; i < len; i++)
    if (m->b[i]) {
  #if dimension > 2
      for (int j = 0; j < len; j++)
	if (m->b[i][j])
	  munmap (m->b[i][j], m->len*sizeof (char *));
  #endif // dimension > 2
      munmap (m->b[i], m->len*sizeof (char *));
    }
#endif // dimension > 1
  if (m->b)
    munmap (m->b, m->len*sizeof (char *));
  free (m);
}

/**
The `mem_assign()` function assigns a (pointer) value to a given index. */
#define mymap(len) \
  mmap (NULL, len*sizeof (char *), PROT_READ|PROT_WRITE, MAP_PRIVATE|MAP_ANONYMOUS, -1, 0)

#if dimension == 1
void mem_assign (Memindex m, int i, int len, void * b)
{
  if (!m->b) {
    m->b = mymap (len);
    m->len = len;
  }
  mem_data(m,i) = b;
}
#elif dimension == 2
void mem_assign (Memindex m, int i, int j, int len, void * b)
{
  assert (b != NULL);
  if (!m->b) {
    m->b = mymap (len);
    m->len = len;
  }
  if (!m->b[i])
    m->b[i] = mymap (len);
  mem_data(m,i,j) = b;
}
#else // dimension == 3
void mem_assign (Memindex m, int i, int j, int k, int len, void * b)
{
  if (!m->b) {
    m->b = mymap (len);
    m->len = len;
  }
  if (!m->b[i])
    m->b[i] = mymap (len);
  if (!m->b[i][j])
    m->b[i][j] = mymap (len);  
  mem_data(m,i,j,k) = b;
}
#endif // dimension == 3

/**
The `mem_free()` function frees a given index. */

static int mem_free1d (char ** const b, const int i, int len)
{
 // Assumes 4k pages = 512*sizeof(char *) = 512*8 = 4096
  char ** const start = b + (i/512)*512;
  if (len > 512) len = 512;
  for (char ** j = start; j < start + len; j++)
    if (*j)
      return 0;
  if (madvise (start, 4096, MADV_DONTNEED))
    perror ("virtual.h: madvise: ");
  return 1;
}

#if dimension == 1
void mem_free (Memindex m, int i, int len, void * b)
{
  (void) b;
  if (m->b[i]) {
    m->b[i] = NULL;
    mem_free1d (m->b, i, len);
  }
}
#elif dimension == 2
void mem_free (Memindex m, int i, int j, int len, void * b)
{
  (void) b;
  if (m->b[i][j]) {
    m->b[i][j] = NULL;
    if (mem_free1d (m->b[i], j, len))
      mem_free1d ((char **)m->b, i, len);
  }
}
#else // dimension == 3
void mem_free (Memindex m, int i, int j, int k, int len, void * b)
{
  (void) b;
  if (m->b[i][j][k]) {
    m->b[i][j][k] = NULL;
    if (mem_free1d (m->b[i][j], k, len))
      if (mem_free1d ((char **)m->b[i], j, len))
	mem_free1d ((char **)m->b, i, len);
  }
}
#endif // dimension == 3

#if 0
/**
The `foreach_mem()` macro traverses every `_i` allocated elements of
array `_m` taking into account a periodicity of `_len` (and ghost
cells). */

macro foreach_mem (Memindex index, int len, int _i) {
  Memindex _m = index;
  int _len = len;
  Point point = {0};
  for (point.i = max(Period.x*GHOSTS, 0);
       point.i < min(_len - Period.x*GHOSTS, _len);
       point.i += _i)
    if (_m->b[point.i])
#if dimension > 1
      for (point.j = max(Period.y*GHOSTS, 0);
	   point.j < min(_len - Period.y*GHOSTS, _len);
	   point.j += _i)
	if (_m->b[point.i][point.j])
#if dimension > 2
	  for (point.k = max(Period.z*GHOSTS, 0);
	       point.k < min(_len - Period.z*GHOSTS, _len);
	       point.k += _i)
	    if (_m->b[point.i][point.j][point.k])
#endif // dimension > 2
#endif // dimension > 1
	      {...}
}
#endif
