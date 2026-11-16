/* reduction () with OpenMP: its clauses on the loop directives of the two
 * @loops (a parallel for over blocks, a simd over cells, a target one gets
 * map (tofrom: ...) too), not on the body's own. Many blocks and threads:
 * exact integer sums, a max. Build with -dimension=2 -fopenmp. */
#include <stdio.h>
#include "stencil.h"

scalar s;

#define NB 64

@loop omp_blocks () : sblocks, i, j, level
{
  _Pragma ("omp parallel for schedule (dynamic)")
  for (int b = 0; b < NB; b++) {
    int i = b % NLEAF, j = 0, level = 2;
    (void) i, (void) j, (void) level;
    void *** sblocks = groups [i];
    {...}
  }
}

@loop omp_cells (BlockShape shape) : c, shape, ci, cj
{
  for (int cj = shape.g; cj < shape.g + shape.m; cj++) {
    _Pragma ("omp simd")
    for (int ci = shape.ox; ci < shape.ox + shape.m; ci++) {
      int c = ci + cj * shape.sy;
      {...}
    }
  }
}

@loop omp_device (BlockShape shape) : c, shape, ci, cj
{
  _Pragma ("omp target teams distribute parallel for")
  for (int cj = shape.g; cj < shape.g + shape.m; cj++)
    for (int ci = shape.ox; ci < shape.ox + shape.m; ci++) {
      int c = ci + cj * shape.sy;
      {...}
    }
}

@foreach omp_foreach () = omp_blocks () omp_cells (tree_block_shape ());
@foreach omp_target () = omp_blocks () omp_device (tree_block_shape ());

int main (void)
{
  for (int l = 0; l < NLEAF; l++)
    leaf_blocks (l);
  foreach ()
    s[] = 1 + (int) (16 * x);           /* 1 .. 12 */
  long n = 0;
  double hi = 0;
  omp_foreach (reduction (+:n) reduction (max:hi)) {
    int t = 0;
    _Pragma ("omp simd")
    for (int q = 0; q < 4; q++)
      t += q;                           /* the body's own: no clause */
    n += (long) s[] + t - 6;
    hi = s[] > hi ? s[] : hi;
  }
  printf ("%ld %g\n", n, hi);
#ifdef TARGET_CHECK
  /* translated only: the target directive's clauses */
  long m = 0;
  omp_target (reduction (+:m))
    m += (long) s[];
  printf ("%ld\n", m);
#endif
  return 0;
}
