/* No event: events_run () has nothing to run, and returns 0. */
#include <stdio.h>
#include <cwrap/event.h>

int main (void)
{
  int rc = events_run ();

  events_print (stdout);
  printf ("rc %d i %d\n", rc, events_i);
  return 0;
}
