#ifndef _BAMR_H_
#define _BAMR_H_

  #include "memindex.h"

  /*
  .. 2:1 balanced, [block based] Adaptive Mesh Refinement (AMR).
  .. Note that the AMR grid is a global grid, and you can't
  .. create multiple tree in a single run 
  */

  #define BLOCK_AMR        1
  #define BALANCE_2_1      1
  #ifndef BAMR_MAX_DEPTH
  #define BAMR_MAX_DEPTH   20
  #endif

  #if dimension == 1
    # define IJK           i
    # define IJKargs       i
  #elif dimension == 2
    # define IJK           i, j
    # define IJKargs       i, int j
  #else
    # define IJK           i, j, k
    # define IJKargs       i, int j, int k
  #endif
  #define IJKL             IJK, level
  #define IJKLargs         IJKargs, int level

  extern struct _Memindex grid_info_layers [BAMR_MAX_DEPTH + 1];

  typedef struct
  {
    int IJKL;
  } Index;

  typedef struct
  {
    Index * indices;
    int n;
  } IndexStack;

  typedef enum
  {
    BAMR_NODE_LEAF   = 1,
    BAMR_NODE_IS_PRUNED = 2,
    BAMR_NODE_HAS_CHILDREN = 4,
  } BAMRNodeFlags;

  typedef struct
  {
    int flags;
    int pid;
  } BAMRNode;

  #if dimension == 1
    #define tree_node(i, level)                                   \
      ((BAMRNode *) mem_data (&grid_info_layers [level], i))
  #elif dimension == 2
    #define tree_node(i, j, level)                                \
      ((BAMRNode *) mem_data (&grid_info_layers [level], i, j))
  #else
    #define tree_node(i, j, k, level)                             \
      ((BAMRNode *) mem_data (&grid_info_layers [level], i, j, k))
  #endif

  #define tree_has_children(n)                                    \
    ((n->flags & BAMR_NODE_HAS_CHILDREN) &&                       \
    !(n->flags & BAMR_NODE_IS_PRUNED))
  
  #if dimension == 1
    #define tree_is_allocated(i, level)                          \
      (level <= BAMR_MAX_DEPTH &&                                \
        mem_allocated (&grid_info_layers [level], i))
  #elif dimension == 2
    #define tree_is_allocated(i, j, level)                       \
      (level <= BAMR_MAX_DEPTH &&                                \
        mem_allocated (&grid_info_layers [level], i, j))
  #else
    #define tree_is_allocated(i, j, k, level)                    \
      (level <= BAMR_MAX_DEPTH &&                                \
        mem_allocated (&grid_info_layers [level], i, j, k))
  #endif

  #if dimension == 1
    #define BAMR_NODE_ENV()                                      \
        BAMRNode * node = tree_is_allocated (i, level) ?         \
          tree_node (i, level) : NULL;                           \
        (void) node
  #elif dimension == 2
    #define BAMR_NODE_ENV()                                      \
        BAMRNode * node = tree_is_allocated (i, j, level) ?      \
          tree_node (i, j, level) : NULL;                        \
        (void) node
  #else
    #define BAMR_NODE_ENV()                                      \
        BAMRNode * node = tree_is_allocated (i, j, k, level) ?   \
          tree_node (i, j, k, level) : NULL;                     \
        (void) node
  #endif

  /*
  .. API (s)
  */
  int  tree (int depthmax, int blocksize);
  void tree_shift_resize (double x0, double y0, double z0, double h0);
  void tree_destroy ();
  
#endif
