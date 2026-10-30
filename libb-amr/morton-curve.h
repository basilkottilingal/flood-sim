#ifndef _QUADTREE_MORTON_CURVE_H_
#define _QUADTREE_MORTON_CURVE_H_

  #include <assert.h>
#include<stdio.h>
  #include "tree.h"
  #include "liballoc/array.h"

  #if dimension == 2

  #define morton_start(_i_,_j_,_l_) do {                           \
    int _ii_ = _i_, _jj_ = _j_, _ll_ = _l_;                        \
    do {                                                           \
      int i = _ii_, j = _jj_, level = _ll_;                        \
      BAMR_NODE_ENV ();                                            \
      if (!node) break;                                            \
      struct _index_ {int i, j, level, visited;}                   \
        stack [(1 << dimension) * BAMR_MAX_DEPTH + 1];             \
      int depth = 0, visited = 0;                                  \
      do                                                           \
      {                                                            \
        while (!visited)                                           \
        {                                                          \
          BAMR_NODE_ENV ();                                        \
          assert (node);                                           \
          if (!tree_has_children (node))                           \
            break;                                                 \
          stack [depth++] = (struct _index_) {i,   j,   level, 1}; \
          i *= 2, j *= 2, level++;                                 \
          assert (level <= BAMR_MAX_DEPTH);                        \
          stack [depth++] = (struct _index_) {i+1, j+1, level, 0}; \
          stack [depth++] = (struct _index_) {i+1, j,   level, 0}; \
          stack [depth++] = (struct _index_) {i,   j+1, level, 0}; \
        }                                                          \
        do {                                                       \
          BAMR_NODE_ENV ();                                        \
          /* do something here */
    #define morton_end()                                           \
        } while (0);                                               \
        if (!depth)                                                \
          break;                                                   \
        depth--;                                                   \
        i = stack [depth].i;                                       \
        j = stack [depth].j;                                       \
        level   = stack [depth].level;                             \
        visited = stack [depth].visited;                           \
      } while (1);                                                 \
    } while (0);                                                   \
  } while (0)

  #endif

  /* API */
  Array * tree_block_indices ();
  
#endif
