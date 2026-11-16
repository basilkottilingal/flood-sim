/* reduction (): +, *, min, max over the 48 cells of stencil.h's three
 * leaves, exact; the value before the loop counted once, on rank 0 only
 * (tree_rank () and tree_reduce () mocked: each call logged); a loop with
 * a parameter; a loop by axis, reduced after each axis. Build with
 * -dimension=2. */
#include <stdio.h>
#include "stencil.h"

scalar s;

/* the processes, mocked: this one's rank, and each tree_reduce () logged
 * as type:op=value */
static int rank;
static char log_ [512];
static inline int tree_rank (void) { return rank; }
static inline void tree_reduce (void * v, int type, int op)
{
  char * end = log_;
  while (*end)
    end++;
  if (type == 7)
    sprintf (end, " %d:%d=%g", type, op, *(double *) v);
  else if (type == 1)
    sprintf (end, " %d:%d=%ld", type, op, *(long *) v);
  else
    sprintf (end, " %d:%d=%d", type, op, *(int *) v);
}

@foreach mock_face_x () = foreach_leaf () foreach_cell (tree_block_shape (), 0);
@foreach mock_face_y () = foreach_leaf () foreach_cell (tree_block_shape (), 0);

static void run (void)
{
  double sum = 5, lo = 1, hi = -1;
  long n = 0;
  int prod = 3, k = 0, f = 0;
  log_ [0] = '\0';
  foreach (reduction (+:sum, n) reduction (*:prod) reduction (min:lo)
           reduction (max:hi)) {
    sum += s[];
    n += s[] > 0;
    prod *= 1 + (x < 0.0625);
    lo = x < lo ? x : lo;
    hi = x > hi ? x : hi;
  }
  foreach_mock_level (2, reduction (+:k))
    k += 2 * (s[] > 0);
  mock_face (reduction (+:f))
    f += s[] > 0;
  printf ("rank %d: %g %ld %d %g %g %d %d;%s\n", rank, sum, n, prod, lo, hi,
          k, f, log_);
}

int main (void)
{
  foreach ()
    s[] = x;
  run ();
  rank = 1;
  run ();
  return 0;
}
