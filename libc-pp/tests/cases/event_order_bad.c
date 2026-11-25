/* Orders events_run () refuses, one per macro: it reports the error and
 * returns -1, with nothing run. */
#include <stdio.h>
#include <cwrap/event.h>

#if UNKNOWN
@event a (i++, next.nothere) {}
@event b (i++, next.a) {}
#elif RESERVED
@event a (i++, prev.end) {}
#elif CYCLE
@event a (i++, next.c) {}
@event b (i++, next.a) {}
@event c (i++, next.b) {}
@event d (i++) {}
#elif SELF
@event a (i++, next.a) {}
#endif
@event init { printf ("ran\n"); }

int main (void)
{
  printf ("rc %d\n", events_run ());
  return 0;
}
