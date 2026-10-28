#ifndef _QUADTREE_MORTON_CURVE_H_
#define _QUADTREE_MORTON_CURVE_H_

  #include <assert.h>
  #include "quadtree.h"
  #include "array.h"

  #define quadtree_morton_start(i,j,l) do {                 \
    struct _index_ {int i, j, level, visited;}              \
      stack [ (1 << dimension) * QAMR_MAX_DEPTH + 1];       \
    int depth = 0, v = 0;                                   \
    stack [depth++] = (struct _index_) {i, j, l, 0};        \
    do                                                      \
    {                                                       \
      while (!v && quadtree_has_children (i, j, l))         \
      {                                                     \
        i *= 2, j *= 2, l++;                                \
        assert (l <= QAMR_MAX_DEPTH);                       \
        stack [depth-1].visited = 1;                        \
        stack [depth++] = (struct _index_) {i+1, j+1, l, 0};\
        stack [depth++] = (struct _index_) {i+1, j,   l, 0};\
        stack [depth++] = (struct _index_) {i,   j+1, l, 0};\
        stack [depth++] = (struct _index_) {i,   j,   l, 0};\
      }                                                     \
      /* do something here */
  #define quadtree_morton_end(i,j,l)                        \
      if (depth) {                                          \
        depth--;                                            \
        i = stack [depth].i;                                \
        j = stack [depth].j;                                \
        l = stack [depth].level;                            \
        v = stack [depth].visited;                          \
      }                                                     \
    } while (depth);                                        \
  } while (0)

  /* API */
  Array * quadtree_block_indices ();
  
#endif
