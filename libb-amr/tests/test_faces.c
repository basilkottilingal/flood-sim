/* Face vectors (libc-pp's face_vector) on the tree: foreach_face (), each
 * face made once, the fluxes made to agree at level jumps (refluxing) when
 * the ghosts of the face vector are made, lazily, by cwrap; face
 * restriction and prolongation. */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cwrap/cache.h>
#include "adapt.h"
#include "block-ghosts.h"

scalar s, lin;
face_vector F, G, P;

static int checks, fails;
#define CHECK(c) do {                                                      \
    checks++;                                                              \
    if (!(c)) {                                                            \
      fails++;                                                             \
      fprintf (stderr, "FAIL line %d: %s\n", __LINE__, #c);                \
    }                                                                      \
  } while (0)

static const double tol = 1e3 * REAL_EPSILON;

/* refine the leaves coarser than 'to' near (0.3, 0.3, 0.3): level jumps */
static void refine_near (int to)
{
  for (int pass = 0; pass < 3 * to; pass++) {
    cache (Index) todo;
    foreach_leaf () {
      double n = 1 << level, r2 = 0, c = (i + 0.5) / n - 0.3;
      r2 += c * c;
#if dimension > 1
      c = (j + 0.5) / n - 0.3, r2 += c * c;
#endif
#if dimension > 2
      c = (k + 0.5) / n - 0.3, r2 += c * c;
#endif
      if (level < to && sqrt (r2) < 0.2)
        cache_push (todo, (Index) {IJKL});
    }
    if (!todo.n)
      break;
    foreach (p in todo)
      tree_refine (p);
  }
}

/* the net flux out of every leaf, times its face area, summed: what the
 * walls let out (none here) if every face is made once; against the sum of
 * every |flux| times area, the scale of its round-off */
static double net_flux (face_vector f)
{
  double total = 0, scale = 0;
  foreach () {
    double out = 0, area = 1;
    foreach_dimension ()
      out += f.x[1] - f.x[];
    for (int a = 1; a < dimension; a++)
      area *= Delta;
    total += out * area;
    foreach_dimension ()
      scale += fabs (f.x[]) * area;
  }
  return fabs (total) / (scale > 0 ? scale : 1);
}

int main (void)
{
  int deep = dimension == 3 ? 3 : 4;
  tree (deep + 1, 4);
  refine_near (deep);
  CHECK (tree_get_depth () == deep && tree_get_depth_min () < deep);
  double L = tree_get_size ();

  /* 1. a random field, fluxes of its gradient, none through the walls:
   * nothing created or lost, level jumps included */
  srand (7);
  foreach ()
    s[] = (Real) (rand () / (double) RAND_MAX);
  foreach_face ()
    F.x[] = x < 0.5 * Delta || x > L - 0.5 * Delta ? 0 : (s[] - s[-1]) / Delta;
  double net = net_flux (F);
  printf ("  levels %d to %d, %d leaves: net flux %.3g of the sum of fluxes\n",
          tree_get_depth_min (), tree_get_depth (), tree_get_nleaves (), net);
  CHECK (net < tol);

  /* 2. a linear field: every face its gradient, at level jumps too */
  foreach () {
    double v = 1 * x;
#if dimension > 1
    v += 2 * y;
#endif
#if dimension > 2
    v += 3 * z;
#endif
    lin[] = (Real) v;
  }
  foreach_face ()
    G.x[] = (lin[] - lin[-1]) / Delta;
  /* each cell's lower faces, but on a wall, where the ghost is a mirror */
  double err = 0, g_x = 1, g_y = 2, g_z = 3;
  (void) g_y, (void) g_z;
  foreach ()
    foreach_dimension ()
      if (x > Delta)
        err = fmax (err, fabs (G.x[] - g_x));
  printf ("  linear field: largest face gradient error %.3g\n", err);
  CHECK (err < 100 * tol);

  /* 3. the positions in foreach_face (): the face's centre */
  foreach_face ()
    P.x[] = (Real) x;
  double perr = 0;
  foreach ()
    perr = fmax (perr, fabs (P.x[] - (x - 0.5 * Delta)));
  printf ("  face positions: error %.3g\n", perr);
  CHECK (perr < tol);

  /* 3b. read at the cell only, a face vector is made to agree first too
   * (ghosts_update_dirty_located (), written by cwrap): Q.x = the size of
   * the cell that made the face; a coarse cell whose lower face meets a
   * finer region reads the fine cells' size there */
  foreach_face ()
    P.x[] = (Real) Delta;
  int finer = 0;
  foreach ()
    finer += P.x[] < 0.75 * Delta;
  printf ("  read at the cell: %d lower faces from a finer region\n", finer);
  CHECK (finer > 0);

  /* 4. refined further: the fluxes prolonged to the new faces agree with
   * the old, and still balance */
  refine_near (deep + 1);
  net = net_flux (F);
  printf ("  refined to level %d: net flux %.3g\n", tree_get_depth (), net);
  CHECK (tree_get_depth () == deep + 1 && net < tol);

  tree_destroy ();
  printf ("dimension %d: %d checks, %d failed\n", dimension, checks, fails);
  return fails != 0;
}
