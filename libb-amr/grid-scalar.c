#include <assert.h>

#include "tree.h"
#include "grid-scalar.h"
#include "morton-curve.h"
#include "pool.h"

#define w        SCALAR_ARRAY_WIDTH
#define NGROUPS  ((NSCALAR_MAX + w - 1) / w)

static Pool * pool_table = NULL;   /* a node's NGROUPS group pointers */
static Pool * pool_group = NULL;   /* w block pointers                */
static Pool * pool_block = NULL;   /* one block of doubles            */

static void grid_scalar_new (scalar s)
{
  grid_scalar_alloc_blocks (s);
}

static void grid_scalar_free (scalar s)
{
  grid_scalar_free_blocks (s);
}

void grid_scalar_init (size_t blocksize)
{
  assert (!pool_block);
  pool_table = pool_new (NGROUPS * sizeof (void **));
  pool_group = pool_new (w * sizeof (void *));
  pool_block = pool_new (blocksize * sizeof (double));
  scalar_new_hook  = grid_scalar_new;
  scalar_free_hook = grid_scalar_free;
}

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

static void alloc_block (void *** t, int s)
{
  if (!t [s/w])
    t [s/w] = pool_alloc (pool_group);
  if (!t [s/w][s%w])
    t [s/w][s%w] = pool_alloc (pool_block);
}

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
  for (int g = 0; g < NGROUPS; g++)
    if (t [g])
    {
      for (int b = 0; b < w; b++)
        pool_free (pool_block, t [g][b]);
      pool_free (pool_group, t [g]);
    }
  pool_free (pool_table, t);
}

void grid_scalar_alloc_blocks (scalar s)
{
  if (!pool_block)        /* no tree yet: tree () allocates them */
    return;
  morton_start (TREE_ROOT) {
    alloc_block (sblocks, s);
  } morton_end ();
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
  morton_start (TREE_ROOT) {
    if (sblocks [s/w] && sblocks [s/w][s%w])
    {
      pool_free (pool_block, sblocks [s/w][s%w]);
      sblocks [s/w][s%w] = NULL;
    }
  } morton_end ();
}
