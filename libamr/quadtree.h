#ifndef _QUADTREE_H_
#define _QUADTREE_H_

  /*
  .. 2:1 balanced, [block based] Adaptive Mesh Refinement (AMR).
  .. Note that the AMR grid is a global grid, and you can't
  .. create multiple quadtree in a single run 
  */

  #define dimension    2
  #define BLOCK_AMR    1
  #define _2_1_BALANCE 1

  #ifndef QAMR_MAX_DEPTH
  #define QAMR_MAX_DEPTH 20
  #endif

  typedef struct
  {
    int i, j, l;
  } Index;

  typedef struct
  {
    Index * indices;
    int n;
  } IndexStack;

  /* API (s)
  .. assumes depthmax in [1, 20), blocksize in 2 ^ {1, 2,.., 10},
  .. scalar_size in {sizeof (double), sizeof (int), sizeof (uint8_t)}
  */
  
  int quadtree (int depthmax, int blocksize);
  int scalar   (size_t scalar_size);
  
#endif
