#ifndef _WEBGL_TREE_H_
#define _WEBGL_TREE_H_

  /*
  .. libb-amr's tree in a scene (webgl.h). Built with cwrap, as tree.h is;
  .. 2D for now.
  ..
  .. webgl_tree (S, name, n, s): a blocks layer of the tree's leaves
  .. (tree_caches ()->leaves: under MPI the halo too), each a block of m x m
  .. cells, and on it a field a scalar s[0 ... n - 1], named after the
  .. scalar (scalar_name (), else s0, s1, ...), plus two of each block:
  .. its level and its owner's rank (pid), drawn categorical. The cells'
  .. values are read in the blocks, ghosts not read. The layer's number, or
  .. -1 (not 2D, or n < 0, or s NULL with n > 0).
  */

  #include "tree.h"
  #include "webgl.h"

  int webgl_tree (WebglScene * S, const char * name, int n, const scalar * s);

#endif
