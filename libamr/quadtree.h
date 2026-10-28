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
  # define QAMR_MAX_DEPTH 20
  #endif

  typedef enum
  {
    QUADTREE_NODE_LEAF   = 1,
    QUADTREE_NODE_PRUNED = 2,
    QUADTREE_NODE_CHILDREN = 4,
  } QuadtreeNodeFlags;

  typedef struct
  {
    int flags;
    int pid;
  } QuadtreeNode;

  #if dimension == 1
  # define quadtree_node(i,l)  
  #endif
  #define quadtree_has_children(n)           \
    ((n->flags & QUADTREE_NODE_CHILDREN) &&  \
     !(n->flags & QUADTREE_NODE_PRUNED))


  /* API (s)
  .. assumes depthmax in [1, 20), blocksize in 2 ^ {1, 2,.., 10},
  .. scalar_size in {sizeof (double), sizeof (int), sizeof (uint8_t)}
  */

  int  quadtree (int depthmax, int blocksize);
  void quadtree_shift_translate (double x0, double y0, double h);
  QuadtreeNode quadtree_node (int i, int j, int level); 
  
#endif
