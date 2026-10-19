#ifndef _QUADTREE_MORTON_H_
#define _QUADTREE_MORTON_H_

  #include <assert.h>
  #include "qamr.h"

  #ifndef QAMR_MAX_DEPTH
  #define QAMR_MAX_DEPTH 20
  #endif

  #define QAMR_MAX_STACK_SIZE ( (1 << DIMENSION) * QAMR_MAX_DEPTH )

  #define qamr_morton_start(i,j,l) do {                     \
    struct q { int i, j, l; } stack [QAMR_MAX_STACK_SIZE ]; \
    depth = 0;                                              \
    do                                                      \
    {                                                       \
      assert (depth < Q_AMR_MAX_STACK_SIZE);                \
      while (qamr_child (i, j, l))                          \
      {                                                     \
        i *= 2, j *= 2, l++;                                \
        stack [depth++] = (struct q) {i+1, j+1, l};         \
        stack [depth++] = (struct q) {i+1, j,   l};         \
        stack [depth++] = (struct q) {i,   j+1, l};         \
        stack [depth++] = (struct q) {i,   j,   l};         \
      }                                                     \
      /* do something here */
  #define qamr_morton_end(i,j,l)                            \
                                                            \
      if (depth) {                                          \
        depth--;                                            \
        i = stack [depth].i;                                \
        j = stack [depth].j;                                \
        l = stack [depth].l;                                \
      }                                                     \
    } while (depth)                                         \
  } while (0)

  /* API */
  IndexStack qamr_block_indices ();
  
#endif
