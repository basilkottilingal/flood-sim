/* vof_advection () timed: a sphere (circle) in the vortex of test_vof.c,
 * on a uniform tree, LEVEL (default 3 in 2D, 2 in 3D: 256^2, 64^3 cells),
 * STEPS steps (default 50); the time per step and per cell, and the share
 * of interface cells. One JSON line. */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <cwrap/cache.h>
#include "adapt.h"
#include "block-ghosts.h"
#include "vof.h"
#include "vof-ghosts.h"

scalar c;
face_vector uf;

static double now (void)
{
  struct timespec t;
  clock_gettime (CLOCK_MONOTONIC, &t);
  return t.tv_sec + 1e-9 * t.tv_nsec;
}

static double psi (double a, double b)
{
  return sin (M_PI * a) * sin (M_PI * a) * sin (M_PI * b) * sin (M_PI * b)
         / M_PI;
}

static double A_x (double x, double y, double z)
{
  (void) x;
  return dimension > 2 ? psi (y, z) : 0;
}

static double A_z (double x, double y, double z)
{
  (void) z;
  return psi (x, y);
}

static void set_vortex (void)
{
  foreach_face_x () {
    double h = Delta / 2;
    uf.x[] = (Real) ((A_z (x, y + h, z) - A_z (x, y - h, z)) / Delta);
  }
#if dimension > 1
  foreach_face_y () {
    double h = Delta / 2;
    uf.y[] = (Real) (((A_x (x, y, z + h) - A_x (x, y, z - h))
                      - (A_z (x + h, y, z) - A_z (x - h, y, z))) / Delta);
  }
#endif
#if dimension > 2
  foreach_face_z () {
    double h = Delta / 2;
    uf.z[] = (Real) (- (A_x (x, y + h, z) - A_x (x, y - h, z)) / Delta);
  }
#endif
}

int main (void)
{
  const char * e = getenv ("LEVEL");
  int level = e ? atoi (e) : dimension == 3 ? 2 : dimension == 2 ? 3 : 6;
  e = getenv ("STEPS");
  int steps = e ? atoi (e) : 50;
  ghost_ops (c, prolongation (prolong_fraction));
  tree (level + 1, 4);
  for (int l = 0; l < level; l++) {
    cache (Index) todo;
    foreach_leaf ()
      if (level < l + 1)
        cache_push (todo, (Index) {IJKL});
    foreach (p in todo)
      tree_refine (p);
  }
  long cells = 0;
  /* the sphere's fraction of each cell, by 4^dimension points */
  foreach () {
    int in = 0, all = 0;
    for (int a = 0; a < 4; a++)
      for (int b = 0; b < (dimension > 1 ? 4 : 1); b++)
        for (int e = 0; e < (dimension > 2 ? 4 : 1); e++) {
          double q [3] = { x + ((a + 0.5) / 4 - 0.5) * Delta - 0.5,
                           y + ((b + 0.5) / 4 - 0.5) * Delta - 0.7,
                           z + ((e + 0.5) / 4 - 0.5) * Delta - 0.5 }, d = 0;
          for (int k = 0; k < dimension; k++)
            d += q [k] * q [k];
          in += d < 0.15 * 0.15;
          all++;
        }
    c[] = (Real) in / all;
    cells++;
  }
  set_vortex ();
  double dt = 0.2 * tree_get_size () / (tree_block_shape ().m << tree_get_depth ());
  vof_advection (c, uf, dt);              /* warm up */
  double t0 = now ();
  for (int s = 0; s < steps; s++)
    vof_advection (c, uf, dt);
  double t = (now () - t0) / steps;
  long interface = 0;
  foreach ()
    interface += c[] > 1e-6 && c[] < 1 - 1e-6;
  printf ("{\"dimension\": %d, \"cells\": %ld, \"interface\": %.4f, "
          "\"ms_per_step\": %.3f, \"ns_per_cell\": %.2f}\n", dimension, cells,
          (double) interface / cells, 1e3 * t, 1e9 * t / cells);
  return 0;
}
