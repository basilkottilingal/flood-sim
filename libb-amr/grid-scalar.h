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
  void *** grid_scalar_table_new (void);    /* a node's table, empty */
  void    grid_scalar_root_blocks (void *** t);  /* the root: a block of
                                                    its own per scalar */
  /*
  .. A family of TREE_NCHILDREN children, tables t [0..]: for each scalar in
  .. use, one object holding the blocks of all of them, child c's the c-th.
  .. uninit: the caller writes every cell, so no zeroing. Freed together,
  .. tables and all, when the family is coarsened.
  */
  void    grid_scalar_family_new (void *** t [], int uninit);
  void    grid_scalar_family_free (void *** t []);
  void    grid_scalar_alloc_blocks (scalar s);  /* in every node */
  void    grid_scalar_free_blocks (scalar s);   /* in every node */

  /*
  .. Debug: counts what the tree's nodes hold and compares it with what the
  .. pools have handed out and not had back (pool_used ()):
  ..   blocks of scalars in use      == the root's, plus TREE_NCHILDREN per
  ..                                    family out of the family pool
  ..   blocks in every slot          == the same (none left by a freed scalar)
  ..   groups, tables                == out of their pools
  .. Prints each mismatch to stderr, then returns 0; returns 1 if all agree.
  .. Its cost is a full traversal: for checks, not for every step.
  */
  long    grid_scalar_blocks_attached (void);   /* blocks of scalars in use */

  /*
  .. Freed blocks stay in memory, so that the next refinement reuses them
  .. with no page fault (pool.h). After a coarsening, this gives back the
  .. pages of all but 'keep' of them.
  */
  void    grid_scalar_trim (size_t keep);
  int     grid_scalar_check (void);

#endif
