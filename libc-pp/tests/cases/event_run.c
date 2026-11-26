/* Events in one file: the order, overloads, the phases, landing on times. */
#include <stdio.h>
#include <cwrap/event.h>

static char log_ [4096];
static int n_;
#define LOG(...) (n_ += snprintf (log_ + n_, sizeof (log_) - n_, __VA_ARGS__))

@event defaults { LOG ("D "); }
@event init () { LOG ("I "); }
@event stability (i++, t <= 1) { events_dtmax (0.3); LOG ("s%d ", i); }
@event vof (i++, next.stability) { LOG ("v "); }
@event vof (first) { LOG ("V "); }
@event output (t += 0.5,
               prev.vof)
{
  LOG ("o%g@%d ", t, __LINE__);
}
@event vof { LOG ("w "); }
@event init (first) { LOG ("J "); }
@event end () { LOG ("E%d", event_status ()); }

int main (void)
{
  int rc;

  events_print (stdout);
  rc = events_run ();
  printf ("%s\nrc %d i %d t %g\n", log_, rc, events_i, events_t);
  return 0;
}
