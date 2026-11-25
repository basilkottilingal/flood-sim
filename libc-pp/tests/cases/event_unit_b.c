/* The second file of event_unit_a.c */
#include <stdio.h>
#include <cwrap/event.h>
#if TWICE
#include "event_init.h"
#endif

#if NEXT
@event y (i = 0, next.x) { printf ("y\n"); }
#else
@event y (i = 0) { printf ("y\n"); }
#endif
