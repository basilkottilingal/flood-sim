/* A solver as a header: its events, compiled in the one file including it. */
#ifndef EVENT_SOLVER_H
#define EVENT_SOLVER_H
#include "event_log.h"

@event init { LOG ("init.solver "); }
@event flux (i++, i <= 2) { LOG ("flux%d ", i); }
#endif
