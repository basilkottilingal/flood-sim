/* dt limits: the lowest of three bodies' wins, the step lands on 0.65, dt
 * read twice is the same. MODE: 1 a limit after dt is fixed, 2 a dt of 0,
 * 3 a limit in init, 4 dt read with no limit and no time ahead, 5 a limit
 * in an event after the one that fixed dt; each an error that fails its
 * body. */
#include <stdio.h>
#include <stdlib.h>
#include <cwrap/event.h>

static char log_ [4096];
static int n_, mode;
#define LOG(...) (n_ += snprintf (log_ + n_, sizeof (log_) - n_, __VA_ARGS__))

@event init
{
  mode = atoi (getenv ("MODE") ? getenv ("MODE") : "0");
  if (mode == 3)
    events_dtmax (1);
}
@event limits (i++, i <= 3)
{
  if (!(mode == 4 && i == 3))
    events_dtmax (mode == 2 && i == 2 ? 0 : 0.5);
}
@event limits { if (!(mode == 4 && i == 3)) events_dtmax (0.2 + 0.1 * i); }
@event limits (last) { if (!(mode == 4 && i == 3)) events_dtmax (0.45); }
@event out (t = 0.65) { LOG ("out%.17g ", t); }
@event advance (i++, next.limits)
{
  double a = events_dtnext (), b = events_dtnext ();

  LOG ("i%d t%g dt%g%s ", i, t, a, a == b ? "" : "!");
  if (mode == 1 && i == 1)
    events_dtmax (0.01);
}
@event late (i++, next.advance) { if (mode == 5 && i == 2) events_dtmax (0.1); }
@event end { LOG ("E%d", event_status ()); }

int main (void)
{
  int rc = events_run ();

  printf ("%s\nrc %d i %d t %g dt %g\n", log_, rc, events_i, events_t,
          events_dt);
  return 0;
}
