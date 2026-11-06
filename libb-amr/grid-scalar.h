#ifndef _BAMR_SCALAR_H_
#define _BAMR_SCALAR_H_

  #include <stddef.h>

  /*
  .. The blocks of scalars: each node has a table of NSCALAR_MAX / w groups,
  .. each of w block pointers, so node->s [id/w][id%w] is the block of
  .. scalar id. A group is allocated when the first of its scalars needs one.
  .. Blocks, groups and tables come from pools (pool.h).
  ..
  .. grid_scalar_init () sets scalar_new_hook and scalar_free_hook: a scalar
  .. made while the tree exists gets a block in every node, and freeing a
  .. scalar frees its blocks. Scalars that exist when tree () runs get theirs
  .. as the nodes are made.
  */

  /* 64 scalars : a practical limit (atleast for now) */
  # define SCALAR_ARRAY_WIDTH   8

  /* api */
  void    grid_scalar_init (size_t blocksize);  /* Reals per block */
  void    grid_scalar_destroy (void);           /* frees every block */
  void *** grid_scalar_node_new (void);  /* a node's table, with a block
                                            for each scalar in use */
  void    grid_scalar_node_free (void *** s);
  void    grid_scalar_alloc_blocks (scalar s);  /* in every node */
  void    grid_scalar_free_blocks (scalar s);   /* in every node */

  /*
  .. Debug: counts what the tree's nodes hold and compares it with what the
  .. pools have handed out and not had back (pool_used ()):
  ..   blocks of scalars in use      == blocks out of the block pool
  ..   blocks in every slot          == the same (none left by a freed scalar)
  ..   groups, tables                == out of their pools
  .. Prints each mismatch to stderr, then returns 0; returns 1 if all agree.
  .. Its cost is a full traversal: for checks, not for every step.
  */
  long    grid_scalar_blocks_attached (void);   /* blocks of scalars in use */
  int     grid_scalar_check (void);

#endif
