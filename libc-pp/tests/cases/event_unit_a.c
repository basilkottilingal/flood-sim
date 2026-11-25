/* Two files with a slot each that nothing orders: an error at startup.
 * NEXT orders them; TWICE includes event_init.h here and in event_unit_b.c. */
#include <stdio.h>
#include <cwrap/event.h>
#if TWICE
#include "event_init.h"
#endif

@event x (i = 0) { printf ("x\n"); }

int main (void)
{
  printf ("rc %d\n", events_run ());
  return 0;
}
