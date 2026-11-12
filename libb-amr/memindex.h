/*
 * Portions of this header and corresponding definition are derived from
 * basilik <https://basilisk.fr/>
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


#ifndef _MEMINDEX_H_
#define _MEMINDEX_H_

  #include <stddef.h>                     /* NULL */

  #ifndef dimension
  #  error "dimension must be defined as 1, 2 or 3"
  #elif (dimension < 1) || (dimension > 3)
  #  error "dimension must be 1, 2 or 3"
  #endif

  struct _Memindex {
    int len;
    #if dimension == 1
    char ** b;
    #elif dimension == 2  
    char *** b;
    #else
    char **** b;
    #endif
  };
  #define Memindex struct _Memindex *

  #if dimension == 1
  void mem_free (Memindex m, int i, int len, void * b);
  void mem_assign (Memindex m, int i, int len, void * b);
  inline int mem_allocated (Memindex m, int i)
  {
    return (i >= 0 && i < m->len && m->b[i]);
  }
  #define mem_data(m,i) ((m)->b[i])
  #elif dimension == 2
  void mem_free (Memindex m, int i, int j, int len, void * b);
  void mem_assign (Memindex m, int i, int j, int len, void * b);
  inline int mem_allocated (Memindex m, int i, int j)
  {
    return (i >= 0 && i < m->len && m->b[i] &&
	    j >= 0 && j < m->len && m->b[i][j]);
  }
  #define mem_data(m,i,j) ((m)->b[i][j])
  #else // dimension == 3
  void mem_free (Memindex m, int i, int j, int k, int len, void * b);
  void mem_assign (Memindex m, int i, int j, int k, int len, void * b);
  inline int mem_allocated (Memindex m, int i, int j, int k)
  {
    return (i >= 0 && i < m->len && m->b[i] &&
	    j >= 0 && j < m->len && m->b[i][j] &&
	    k >= 0 && k < m->len && m->b[i][j][k]);
  }
  #define mem_data(m,i,j,k) ((m)->b[i][j][k])
  #endif // dimension == 3

  Memindex mem_new (int len);
  void mem_destroy (Memindex m, int len);     /* an empty index is fine */

  /*
  .. The interface shared with voxel-memindex.h and nano-vdb.h (see
  .. level-index.h): mem_get () is the pointer at an index, or NULL, out of
  .. range too. memindex finds an index by arithmetic alone, so its cursor
  .. holds nothing and mem_get_cursor () is mem_get ().
  */
  typedef struct { int unused; } MemCursor;

  #if dimension == 1
  static inline void * mem_get (Memindex m, int i)
  {
    return mem_allocated (m, i) ? m->b[i] : NULL;
  }
  static inline void * mem_get_cursor (Memindex m, MemCursor * c, int i)
  {
    (void) c;
    return mem_get (m, i);
  }
  #elif dimension == 2
  static inline void * mem_get (Memindex m, int i, int j)
  {
    return mem_allocated (m, i, j) ? m->b[i][j] : NULL;
  }
  static inline void * mem_get_cursor (Memindex m, MemCursor * c, int i, int j)
  {
    (void) c;
    return mem_get (m, i, j);
  }
  #else
  static inline void * mem_get (Memindex m, int i, int j, int k)
  {
    return mem_allocated (m, i, j, k) ? m->b[i][j][k] : NULL;
  }
  static inline void * mem_get_cursor (Memindex m, MemCursor * c, int i, int j,
                                       int k)
  {
    (void) c;
    return mem_get (m, i, j, k);
  }
  #endif

#endif
