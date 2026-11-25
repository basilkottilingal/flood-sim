/* An event defined in a file of its own, as a library would: other files
 * declare it (@event adv;). next.flux names an event of another file. Its
 * init is first: two files' init bodies at one rank have no order (TIE).
 * With SOLVER, it includes the solver's events too: compiled twice. */
#ifdef SOLVER
#include "event_solver.h"
#else
#include "event_log.h"
#endif

@event adv (i++, next.flux) { LOG ("adv.lib "); }
#ifndef TIE
@event init (first) { LOG ("init.lib "); }
#else
@event init { LOG ("init.lib "); }
#endif
