#include <stdio.h>
#include "stencil.h"

scalar s, t, q;

int main (void)
{
  /* t: its storage index, and 100 per leaf, ghosts included */
  foreach_leaf () {
    Real * a = block_data (sblocks, t);
    for (int k = 0; k < NX * NX; k++)
      a [k] = k + 100 * i;
  }

  /* a level up: coarse (t), its neighbours, child.x */
  foreach_up ()
    s[] = coarse (t) + child.x * (coarse (t, 1) - coarse (t, -1)) / 8.
          + child.y * (coarse (t, 0, 1) - coarse (t, 0, -1)) / 8. + child.z;
  BlockShape shape = tree_block_shape ();
  double err = 0;
  {
    const Real * u = block_data (leaf_blocks (1), t), * v = block_data (leaf_blocks (0), s);
    for (int cj = shape.g; cj < shape.g + shape.m; cj++)
      for (int ci = shape.ox; ci < shape.ox + shape.m; ci++) {
        int c = ci + cj * shape.sy;
        int cc = shape.ox + (ci - shape.ox) / 2 + (shape.g + (cj - shape.g) / 2) * shape.sy;
        int sx = (ci - shape.ox) % 2 ? 1 : -1, sy = (cj - shape.g) % 2 ? 1 : -1;
        double w = u [cc] + sx * (u [cc + 1] - u [cc - 1]) / 8.
                   + sy * (u [cc + shape.sy] - u [cc - shape.sy]) / 8.;
        err += v [c] - w;
      }
  }

  /* a level down: the mean of the children, leaf 2 from leaves 0, 1, 2, 0 */
  foreach_leaf () {
    Real * a = block_data (sblocks, q);
    for (int k = 0; k < NX * NX; k++)
      a [k] = i;
  }
  foreach_down () {
    double sum = 0;
    foreach_child ()
      sum += q[];
    t[] = sum / 4;
  }
  double sum = 0;
  {
    const Real * a = block_data (leaf_blocks (2), t);
    for (int cj = shape.g; cj < shape.g + shape.m; cj++)
      for (int ci = shape.ox; ci < shape.ox + shape.m; ci++)
        sum += a [ci + cj * shape.sy];
  }
  printf ("%g %g | %s\n", err, sum, ghosts_log);   /* quadrants 0 1 2 0 */
  return 0;
}
