#ifndef _QUADTREE_H_
#define _QUADTREE_H_

  /*
  .. block based AMR. note that the AMR grid is global grid, and you can't
  .. create multiple quadtree in a single run 
  */

  #define DIMENSION 2
  #define BLOCK_AMR 1

  typedef struct
  {
    struct index {int i, j, l;} * indices;
    int n;
  } IndexStack;

  /* API
  .. assumes depthmax in [1, 20), blocksize in 2 ^ {1, 2,.., 10},
  .. scalar_size in {sizeof (double), sizeof (int), sizeof (uint8_t)}
  */
  
  int quadtree (int depthmax, int blocksize);
  int scalar   (size_t scalar_size);
  
#endif
