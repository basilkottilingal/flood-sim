/*
.. The tree's leaves as a blocks layer: see webgl-tree.h. The leaf's node
.. and block come from tree_caches (), a cell by block_at (): no lookup in
.. the level index.
*/

#include <stdio.h>
#include <stdlib.h>

#include "tree.h"
#include "block-ghosts.h"
#include "webgl-tree.h"

int webgl_tree (WebglScene * S, const char * name, int n, const scalar * s)
{
#if CWRAP_DIMENSION != 2
  (void) S, (void) name, (void) n, (void) s;
  return -1;
#else
  if (!S || !name || n < 0 || (n > 0 && !s))
    return -1;
  const TreeCaches * tc = tree_caches ();
  int nl = tc->nleaves, m = tree_get_block_m ();
  double o[2] = { tree_get_origin (0), tree_get_origin (1) };
  double size = tree_get_size ();
  size_t mm = (size_t) m * (size_t) m;
  double * box = malloc (4 * (size_t) (nl ? nl : 1) * sizeof (double));
  float * level = malloc ((size_t) (nl ? nl : 1) * sizeof (float));
  float * pid = malloc ((size_t) (nl ? nl : 1) * sizeof (float));
  float * v = malloc ((size_t) (nl ? nl : 1) * mm * sizeof (float));
  if (!box || !level || !pid || !v) {
    fprintf (stderr, "webgl: out of memory\n");
    abort ();
  }
  for (int q = 0; q < nl; q++) {
    Index p = tc->leaves[q].p;
    double h = size / (double) (1L << p.level);
    box[4 * q] = o[0] + p.i * h, box[4 * q + 1] = o[1] + p.j * h;
    box[4 * q + 2] = box[4 * q + 3] = h;
    level[q] = (float) p.level;
    pid[q] = (float) tc->leaves[q].node->pid;
  }
  int layer = webgl_blocks (S, name, nl, m, box);
  BlockShape sh = tree_block_shape ();
  for (int k = 0; layer >= 0 && k < n; k++) {
    for (int q = 0; q < nl; q++) {
      const Real * b = block_data (tc->leaves[q].node->s, s[k]);
      for (int j = 0; j < m; j++)
        for (int i = 0; i < m; i++) {
          int x[3] = { i, j, 0 };
          v[(size_t) q * mm + (size_t) j * (size_t) m + (size_t) i] =
            (float) b[block_at (sh, x)];
        }
    }
    const char * sn = scalar_name (s[k]);
    char fallback[32];
    if (!sn || !*sn) {
      snprintf (fallback, sizeof (fallback), "s%d", k);
      sn = fallback;
    }
    webgl_field (S, layer, sn, WEBGL_CELL, v, 0);
  }
  if (layer >= 0) {
    webgl_field (S, layer, "level", WEBGL_BLOCK, level, 1);
    webgl_field (S, layer, "pid", WEBGL_BLOCK, pid, 1);
  }
  free (box), free (level), free (pid), free (v);
  return layer;
#endif
}
