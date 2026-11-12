/* Vertex scalars (libc-pp's vertex_scalar) on the tree: foreach_vertex (),
 * each vertex made once by its owner, every other copy its value (ghosts),
 * continuous at level jumps (hanging vertices interpolated), refined and
 * coarsened. */
#include <stdio.h>
#include <math.h>
#include <cwrap/cache.h>
#include "adapt.h"
#include "block-ghosts.h"

vertex_scalar v, w;

static int checks, fails;
#define CHECK(c) do {                                                      \
    checks++;                                                              \
    if (!(c)) {                                                            \
      fails++;                                                             \
      fprintf (stderr, "FAIL line %d: %s\n", __LINE__, #c);                \
    }                                                                      \
  } while (0)

static const double tol = 1e3 * REAL_EPSILON;

static double linear (double x, double y, double z)
{
  (void) y, (void) z;
  double l = 1 + x;
#if dimension > 1
  l += 2 * y;
#endif
#if dimension > 2
  l += 3 * z;
#endif
  return l;
}

static double centre = 0.3;

static void refine_near (int to)
{
  for (int pass = 0; pass < 3 * to; pass++) {
    cache (Index) todo;
    foreach_leaf () {
      double n = 1 << level, r2 = 0, c = (i + 0.5) / n - centre;
      r2 += c * c;
#if dimension > 1
      c = (j + 0.5) / n - centre, r2 += c * c;
#endif
#if dimension > 2
      c = (k + 0.5) / n - centre, r2 += c * c;
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

/* every corner of every leaf cell against the linear field (read into
 * locals first: a comma in v[1,1] would split a macro's argument) */
static double corners_error (void)
{
  double err = 0;
  foreach () {
    double x0 = x - 0.5 * Delta, D = Delta;
#if dimension == 1
    double c [2] = { v[], v[1] };
    double l [2] = { linear (x0, 0, 0), linear (x0 + D, 0, 0) };
#elif dimension == 2
    double y0 = y - 0.5 * Delta;
    double c [4] = { v[], v[1], v[0,1], v[1,1] };
    double l [4] = { linear (x0, y0, 0), linear (x0 + D, y0, 0),
                     linear (x0, y0 + D, 0), linear (x0 + D, y0 + D, 0) };
#else
    double y0 = y - 0.5 * Delta, z0 = z - 0.5 * Delta;
    double c [8] = { v[], v[1], v[0,1], v[1,1],
                     v[0,0,1], v[1,0,1], v[0,1,1], v[1,1,1] };
    double l [8] = { linear (x0, y0, z0), linear (x0 + D, y0, z0),
                     linear (x0, y0 + D, z0), linear (x0 + D, y0 + D, z0),
                     linear (x0, y0, z0 + D), linear (x0 + D, y0, z0 + D),
                     linear (x0, y0 + D, z0 + D),
                     linear (x0 + D, y0 + D, z0 + D) };
#endif
    for (int q = 0; q < (1 << dimension); q++)
      err = fmax (err, fabs (c [q] - l [q]));
  }
  return err;
}

int main (void)
{
  int deep = dimension == 3 ? 3 : 4;
  tree (deep + 1, 4);
  refine_near (deep);
  CHECK (tree_get_depth () == deep && tree_get_depth_min () < deep);

  /* 1. a linear field made at the vertices: every corner of every cell
   * holds it, at level jumps (hanging vertices) and walls too */
  foreach_vertex ()
    v[] = (Real) linear (x, y, z);
  double err = corners_error ();
  printf ("  levels %d to %d, %d leaves: corners' error %.3g\n",
          tree_get_depth_min (), tree_get_depth (), tree_get_nleaves (), err);
  CHECK (err < tol);

  /* 2. positions in foreach_vertex (): the cell's lower corner */
  foreach_vertex ()
    w[] = (Real) x;
  double perr = 0;
  foreach ()
    perr = fmax (perr, fabs (w[] - (x - 0.5 * Delta)));
  printf ("  vertex positions: error %.3g\n", perr);
  CHECK (perr < tol);

  /* 3. refined further: new vertices from the old, still exact */
  refine_near (deep + 1);
  err = corners_error ();
  printf ("  refined to level %d: corners' error %.3g\n", tree_get_depth (), err);
  CHECK (tree_get_depth () == deep + 1 && err < tol);

  /* 4. coarsened: every family of leaves that can be, merged */
  int merged = 0;
  {
    cache (Index) parents;
    foreach_node () {
      if (!tree_has_children (node) || level != tree_get_depth () - 1)
        continue;
      cache_push (parents, (Index) {IJKL});
    }
    foreach (p in parents)
      merged += tree_coarsen (p);
  }
  err = corners_error ();
  printf ("  %d families coarsened: corners' error %.3g\n", merged, err);
  CHECK (merged > 0 && err < tol);

  tree_destroy ();

  /* 5. refined at the upper walls, beside coarser leaves there: the
   * vertices on the walls owned by blocks of every level */
  centre = 0.85;
  tree (deep + 1, 4);
  refine_near (deep);
  foreach_vertex ()
    v[] = (Real) linear (x, y, z);
  err = corners_error ();
  printf ("  refined at the upper walls, levels %d to %d: corners' error %.3g\n",
          tree_get_depth_min (), tree_get_depth (), err);
  CHECK (tree_get_depth_min () < deep && err < tol);
  tree_destroy ();
  printf ("dimension %d: %d checks, %d failed\n", dimension, checks, fails);
  return fails != 0;
}
