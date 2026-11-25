/* An overload of init, with an include guard: fine in one file, linked
 * twice when two files include it. */
#ifndef EVENT_INIT_H
#define EVENT_INIT_H
#include <stdio.h>
#include <cwrap/event.h>
@event init { printf ("init\n"); }
#endif
