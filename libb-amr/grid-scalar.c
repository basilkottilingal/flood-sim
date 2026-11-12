#include <assert.h>
#include <stdio.h>

#include "tree.h"
#include "grid-scalar.h"
#include "morton-curve.h"
#include "pool.h"
#include "block-ghosts.h"

/* the precision this library is built with: see real.h */
const int BAMR_REAL_TAG = 1;

#define w        SCALAR_ARRAY_WIDTH
#define NGROUPS  ((NSCALAR_MAX + w - 1) / w)

static Pool * pool_table  = NULL;  /* a node's NGROUPS group pointers  */
static Pool * pool_group  = NULL;  /* w block pointers                 */
static Pool * pool_block  = NULL;  /* one block of Reals: the root's   */
static Pool * pool_family = NULL;  /* TREE_NCHILDREN blocks, the blocks
                                      of a family of children, at once */
static size_t block_reals = 0;     /* Reals in a block                 */

/*
.. used as the constructor for creating blocks for scalar 's'. It will be
.. hooked to the constructor that allocate 's' (libc-pp/srs/scalar.c)
*/
static void grid_scalar_new (scalar s)
{
  scalar_ops_reset (s);
  grid_scalar_alloc_blocks (s);
}

/*
.. used as the destructor that frees blocks of scalar 's'. It will be
.. hooked to the destructor that deallocate 's' (libc-pp/srs/scalar.c)
*/
static void grid_scalar_free (scalar s)
{
  grid_scalar_free_blocks (s);
  scalar_ops_reset (s);
}

/*
.. (a) attaches constructor/destructor hooks.
.. (b) create pool for blocks, scalar block arrays
*/
void grid_scalar_init (size_t blocksize)
{
  assert (!pool_block);
  block_reals = blocksize;
  pool_table  = pool_new (NGROUPS * sizeof (void **));
  pool_group  = pool_new (w * sizeof (void *));
  pool_block  = pool_new (blocksize * sizeof (Real));
  pool_family = pool_new (TREE_NCHILDREN * blocksize * sizeof (Real));
  scalar_new_hook  = grid_scalar_new;
  scalar_free_hook = grid_scalar_free;
}

/*
.. (a) unhooks constructor/destructor
.. (b) free 3 pools
*/
void grid_scalar_destroy (void)
{
  if (scalar_new_hook == grid_scalar_new)
    scalar_new_hook = NULL;
  if (scalar_free_hook == grid_scalar_free)
    scalar_free_hook = NULL;
  pool_destroy (pool_family);
  pool_destroy (pool_block);
  pool_destroy (pool_group);
  pool_destroy (pool_table);
  pool_family = pool_block = pool_group = pool_table = NULL;
}

/* the group of scalar s in table t, made if missing */
static void ** group_of (void *** t, int s)
{
  if (!t [s/w])
    t [s/w] = pool_alloc (pool_group);
  return t [s/w];
}

void *** grid_scalar_table_new (void)
{
  return pool_alloc (pool_table);
}

/* the root has no siblings: a block of its own for scalar s */
static void root_block (void *** t, int s)
{
  void ** g = group_of (t, s);
  if (!g [s%w])
    g [s%w] = pool_alloc (pool_block);
}

void grid_scalar_root_blocks (void *** t)
{
  foreach_scalar (s)
    root_block (t, s);
}

/*
.. The blocks of scalar s for a family of children, tables t [0..], as one
.. object of the family pool: child c's block is the c-th of it. uninit:
.. the caller writes every cell (prolongation), so no zeroing.
*/
static void family_block (void *** t [TREE_NCHILDREN], int s, int uninit)
{
  Real * base = uninit ? pool_alloc_uninit (pool_family)
                       : pool_alloc (pool_family);
  for (int c = 0; c < TREE_NCHILDREN; c++)
    group_of (t [c], s) [s%w] = base + (size_t) c * block_reals;
}

void grid_scalar_family_new (void *** t [], int uninit)
{
  foreach_scalar (s)
    family_block (t, s, uninit);
}

/* the family's block of scalar s, back to the pool, and its slots emptied */
static void family_block_free (void *** t [TREE_NCHILDREN], int s)
{
  if (!t [0][s/w] || !t [0][s/w][s%w])
    return;
  pool_free (pool_family, t [0][s/w][s%w]);   /* child 0's is the base */
  for (int c = 0; c < TREE_NCHILDREN; c++)
    t [c][s/w][s%w] = NULL;
}

void grid_scalar_family_free (void *** t [])
{
  /*
  .. Only scalars in use have blocks: the free hook takes a scalar's blocks
  .. back from every family (grid_scalar_free_blocks ()).
  */
  foreach_scalar (s)
    family_block_free (t, s);
  for (int c = 0; c < TREE_NCHILDREN; c++)
  {
    for (int g = 0; g < NGROUPS; g++)
      pool_free (pool_group, t [c][g]);
    pool_free (pool_table, t [c]);
  }
}

/* the tables of the children of the node at p */
static void children_tables (Index p, void *** t [TREE_NCHILDREN])
{
  int c = 0;
  children_of (p)
    t [c++] = node->s;
}

void grid_scalar_alloc_blocks (scalar s)
{
  if (!pool_block)        /* no tree yet: tree () allocates them */
    return;
  root_block (index_node (TREE_ROOT_INDEX)->s, s);
  const TreeCaches * tc = tree_caches ();
  for (int q = 0; q < tc->nall; q++)
    if (tree_has_children (tc->all [q].node)) {
      void *** t [TREE_NCHILDREN];
      children_tables (tc->all [q].p, t);
      family_block (t, s, 0);
    }
}

void grid_scalar_free_blocks (scalar s)
{
  /*
  .. A scalar may have no blocks: made before tree () and freed after
  .. tree_destroy ().
  .. Group tables are kept: fixme : may free sblocks [s/w] once empty.
  */
  if (!pool_block)
    return;
  void *** root = index_node (TREE_ROOT_INDEX)->s;
  if (root [s/w] && root [s/w][s%w]) {
    pool_free (pool_block, root [s/w][s%w]);
    root [s/w][s%w] = NULL;
  }
  const TreeCaches * tc = tree_caches ();
  for (int q = 0; q < tc->nall; q++)
    if (tree_has_children (tc->all [q].node)) {
      void *** t [TREE_NCHILDREN];
      children_tables (tc->all [q].p, t);
      family_block_free (t, s);
    }
}

/*
.. used for debugging. counts number of blocks attached. lazy and
.. not optimized. recommended for debugging runs only
*/
long grid_scalar_blocks_attached (void)
{
  long n = 0;
  if (!pool_block)
    return 0;
  foreach_node () {
    foreach_scalar (s)
      if (sblocks [s/w] && sblocks [s/w][s%w])
        n++;
  }
  return n;
}

/*
.. compare number of blocks, scalar table, scalar_arrays in use are the
.. same as their respective pools have released
*/
int grid_scalar_check (void)
{
  if (!pool_block)
    return 1;
  long attached = 0, slots = 0, groups = 0, tables = 0, missing = 0;
  foreach_node () {
    tables++;
    foreach_scalar (s)
      if (sblocks [s/w] && sblocks [s/w][s%w])
        attached++;
      else
        missing++;
    for (int g = 0; g < NGROUPS; g++)
      if (sblocks [g])
      {
        groups++;
        for (int b = 0; b < w; b++)
          slots += sblocks [g][b] != NULL;
      }
  }

  int ok = 1;
  long out;
  #define MISMATCH(what, have, pool, more)                                \
    if ((out = (long) pool_used (pool) + (more)) != (have)) {             \
      fprintf (stderr, "grid_scalar_check: %ld %s attached, but %ld out " \
               "of their pool\n", (long) (have), what, out);             \
      ok = 0;                                                             \
    }
  /* the root's blocks are single, every other block one of a family */
  MISMATCH ("blocks of scalars in use", attached,
            pool_block, TREE_NCHILDREN * (long) pool_used (pool_family));
  MISMATCH ("tables", tables, pool_table, 0);
  MISMATCH ("groups", groups, pool_group, 0);
  #undef MISMATCH
  if (slots != attached)
  {
    fprintf (stderr, "grid_scalar_check: %ld blocks belong to scalars not "
             "in use\n", slots - attached);
    ok = 0;
  }
  if (missing)
  {
    fprintf (stderr, "grid_scalar_check: %ld blocks of scalars in use are "
             "missing\n", missing);
    ok = 0;
  }
  return ok;
}

void grid_scalar_trim (size_t keep)
{
  if (pool_block) {
    pool_trim (pool_block, keep);
    pool_trim (pool_family, keep / TREE_NCHILDREN);
  }
}
