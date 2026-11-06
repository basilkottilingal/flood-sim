#include <assert.h>
#include <stdio.h>

#include "tree.h"
#include "grid-scalar.h"
#include "morton-curve.h"
#include "pool.h"

#define w        SCALAR_ARRAY_WIDTH
#define NGROUPS  ((NSCALAR_MAX + w - 1) / w)

static Pool * pool_table = NULL;   /* a node's NGROUPS group pointers */
static Pool * pool_group = NULL;   /* w block pointers                */
static Pool * pool_block = NULL;   /* one block of doubles            */

/*
.. used as the constructor for creating blocks for scalar 's'. It will be
.. hooked to the constructor that allocate 's' (libc-pp/srs/scalar.c)
*/
static void grid_scalar_new (scalar s)
{
  grid_scalar_alloc_blocks (s);
}

/*
.. used as the destructor that frees blocks of scalar 's'. It will be
.. hooked to the destructor that deallocate 's' (libc-pp/srs/scalar.c)
*/
static void grid_scalar_free (scalar s)
{
  grid_scalar_free_blocks (s);
}

/*
.. (a) attaches constructor/destructor hooks.
.. (b) create pool for blocks, scalar block arrays
*/
void grid_scalar_init (size_t blocksize)
{
  assert (!pool_block);
  pool_table = pool_new (NGROUPS * sizeof (void **));
  pool_group = pool_new (w * sizeof (void *));
  pool_block = pool_new (blocksize * sizeof (double));
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
  pool_destroy (pool_block);
  pool_destroy (pool_group);
  pool_destroy (pool_table);
  pool_block = pool_group = pool_table = NULL;
}

/*
.. alloc block at a node corresponding for scalar 's'
*/
static void alloc_block (void *** t, int s)
{
  if (!t [s/w])
    t [s/w] = pool_alloc (pool_group);
  if (!t [s/w][s%w])
    t [s/w][s%w] = pool_alloc (pool_block);
}

/*
.. alloc blocks for each scalar which will be later attached to scalar table
.. of a b-AMR node
*/
void *** grid_scalar_node_new (void)
{
  void *** t = pool_alloc (pool_table);
  foreach_scalar (s)
    alloc_block (t, s);
  return t;
}

void grid_scalar_node_free (void *** t)
{
  if (!t)
    return;
  /*
  .. Only scalars in use have blocks: the free hook takes a scalar's blocks
  .. back from every node (grid_scalar_free_blocks ()).
  */
  foreach_scalar (s)
    if (t [s/w])
      pool_free (pool_block, t [s/w][s%w]);
  for (int g = 0; g < NGROUPS; g++)
    pool_free (pool_group, t [g]);
  pool_free (pool_table, t);
}

void grid_scalar_alloc_blocks (scalar s)
{
  if (!pool_block)        /* no tree yet: tree () allocates them */
    return;
  morton (TREE_ROOT_INDEX) {
    alloc_block (sblocks, s);
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
  morton (TREE_ROOT_INDEX) {
    if (sblocks [s/w] && sblocks [s/w][s%w])
    {
      pool_free (pool_block, sblocks [s/w][s%w]);
      sblocks [s/w][s%w] = NULL;
    }
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
  morton (TREE_ROOT_INDEX) {
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
  morton (TREE_ROOT_INDEX) {
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
  #define MISMATCH(what, have, pool)                                      \
    if ((out = (long) pool_used (pool)) != (have)) {                      \
      fprintf (stderr, "grid_scalar_check: %ld %s attached, but %ld out " \
               "of their pool\n", (long) (have), what, out);             \
      ok = 0;                                                             \
    }
  MISMATCH ("blocks of scalars in use", attached, pool_block);
  MISMATCH ("tables", tables, pool_table);
  MISMATCH ("groups", groups, pool_group);
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
