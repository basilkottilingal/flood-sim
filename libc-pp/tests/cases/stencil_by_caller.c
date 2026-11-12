#include <stdio.h>
#include "stencil.h"

scalar p, q;              /* ids 0, 1 */

int main (void)
{
  foreach_mock_level (2)
    p[] = 1;                          /* not marked: the caller's */
  foreach_mock_level (2)
    q[] = p[1] - p[-1];               /* not updated */
  printf ("level: updated %s| stale %llu; ", ghosts_log, ghosts_dirty);
  foreach ()
    p[] = 2;                          /* foreach () marks it */
  foreach_mock_level (2)
    q[] = p[1];                       /* and still no update here */
  printf ("then: updated %s| stale %llu\n", ghosts_log, ghosts_dirty);
  return 0;
}
