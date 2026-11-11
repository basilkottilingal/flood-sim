#include <stdio.h>
#include "stencil.h"

scalar p, q, r;           /* ids 0, 1, 2 */
scalar a0, a1;            /* 3, 4 */
scalar list [2];          /* an array: a0 and a1 below */

int main (void)
{
  list [0] = a0;
  list [1] = a1;
  foreach ()
    p[] = 1;                          /* p stale */
  foreach ()
    q[] = p[1] - p[-1];               /* p updated first; q stale */
  foreach ()
    r[] = p[0,1];                     /* p up to date: nothing */
  foreach ()
    r[] = q[] + q[0,0];               /* the cell only: nothing */
  foreach ()
    for (int k = 0; k < 2; k++)
      list[k][] = k;                  /* a0 and a1 stale */
  foreach ()
    r[] = a1[1];                      /* a1 updated */
  printf ("updated %s| stale %llu\n", ghosts_log, ghosts_dirty);
  return 0;
}
