#ifndef _LEVEL_INDEX_H_
#define _LEVEL_INDEX_H_

  /*
  .. The index of one level of the tree: a sparse (1 << level)^dimension
  .. array of node pointers. Three implementations share one interface, and
  .. BAMR_INDEX chooses the one the tree uses, when compiling:
  ..
  ..   -DBAMR_INDEX=1  memindex.h        per-level virtual arrays (Basilisk):
  ..                                     d loads, pages committed by the OS
  ..                                     (the default)
  ..   -DBAMR_INDEX=2  voxel-memindex.h  a tree of 64-way nodes
  ..   -DBAMR_INDEX=3  nano-vdb.h        a VDB tree: root hash, 32^d, 16^d
  ..                                     and 8^d nodes; -DNVDB_COMPACT for
  ..                                     masks and popcount rank
  ..
  .. The interface, for an index m of type LevelIndex (IJK: i, j, k as the
  .. dimension has them; len: 1 << level):
  ..   lidx_get (m, IJK)                 the pointer, or NULL (out of range too)
  ..   lidx_allocated (m, IJK)           lidx_get () != NULL
  ..   lidx_assign (m, IJK, len, p)      p, not NULL, at IJK
  ..   lidx_free (m, IJK, len, p)        IJK empty again
  ..   lidx_destroy (m)                  everything freed; m empty, reusable
  ..   lidx_get_cursor (m, &c, IJK)      lidx_get () through a cursor c, which
  ..                                     makes nearby lookups cheaper. A
  ..                                     LevelCursor starts as {0}, belongs to
  ..                                     one index, and must be reset ({0})
  ..                                     after any lidx_free () or
  ..                                     lidx_destroy () of that index.
  .. BAMR_INDEX_NAME names the implementation.
  ..
  .. The three can also be used side by side, by their own names (mem_,
  .. vmem_, nvdb_): the benchmarks in tests/ do.
  */

  #define BAMR_INDEX_MEMINDEX  1
  #define BAMR_INDEX_VOXEL     2
  #define BAMR_INDEX_NANOVDB   3

  #ifndef BAMR_INDEX
  #  define BAMR_INDEX  BAMR_INDEX_MEMINDEX
  #endif

  #if BAMR_INDEX == BAMR_INDEX_MEMINDEX
    #include "memindex.h"
    typedef struct _Memindex  LevelIndex;
    typedef MemCursor         LevelCursor;
    #define lidx_get          mem_get
    #define lidx_get_cursor   mem_get_cursor
    #define lidx_assign       mem_assign
    #define lidx_free         mem_free
    #define lidx_destroy_     mem_destroy
    #define BAMR_INDEX_NAME   "memindex"
  #elif BAMR_INDEX == BAMR_INDEX_VOXEL
    #include "voxel-memindex.h"
    typedef struct _VMemindex LevelIndex;
    typedef VMemCursor        LevelCursor;
    #define lidx_get          vmem_get
    #define lidx_get_cursor   vmem_get_cursor
    #define lidx_assign       vmem_assign
    #define lidx_free         vmem_free
    #define lidx_destroy_     vmem_destroy
    #define BAMR_INDEX_NAME   "voxel-memindex"
  #elif BAMR_INDEX == BAMR_INDEX_NANOVDB
    #include "nano-vdb.h"
    typedef struct _NanoVDB   LevelIndex;
    typedef NvdbCursor        LevelCursor;
    #define lidx_get          nvdb_get
    #define lidx_get_cursor   nvdb_get_cursor
    #define lidx_assign       nvdb_assign
    #define lidx_free         nvdb_free
    #define lidx_destroy_     nvdb_destroy
    #ifdef NVDB_COMPACT
    #  define BAMR_INDEX_NAME "nano-vdb (compact)"
    #else
    #  define BAMR_INDEX_NAME "nano-vdb"
    #endif
  #else
  #  error "BAMR_INDEX must be 1 (memindex), 2 (voxel-memindex) or 3 (nano-vdb)"
  #endif

  #define lidx_allocated(m, ...)  (lidx_get (m, __VA_ARGS__) != NULL)

  static inline void lidx_destroy (LevelIndex * m)
  {
    if (m->len)
      lidx_destroy_ (m, m->len);
    *m = (LevelIndex) {0};
  }

#endif
