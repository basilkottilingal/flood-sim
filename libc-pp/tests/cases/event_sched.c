/* Schedules by step and by time, bounds, conditions; times landed on. */
#include <stdio.h>
#include <cwrap/event.h>

static char log_ [4096];
static int n_;
#define LOG(...) (n_ += snprintf (log_ + n_, sizeof (log_) - n_, __VA_ARGS__))

static int every, nout, bad;
static double dtmin = 1, dtmax;

/* read when the steps start: after init */
@event init { every = 3; }
@event a (i = 2) { LOG ("a%d ", i); }
@event b (i = 1, i += every, i <= 10) { LOG ("b%d ", i); }
@event c (i = 1; i += 3; i < 10) { LOG ("c%d ", i); }
@event d (i++, i % 4 == 0, i <= 12) { LOG ("d%d ", i); }
@event e (i % 5 == 0 && i > 0, i <= 12) { LOG ("e%d ", i); }

/* t lands on 0.1 k exactly, on 0.35, and no step is longer than 0.03 */
@event step (i++)
{
  double dt = events_dtnext (0.03);

  dtmin = dt < dtmin ? dt : dtmin;
  dtmax = dt > dtmax ? dt : dtmax;
}
@event out (t += 0.1, t <= 1) { bad += t != nout * 0.1; nout++; }
@event once (t = 0.35) { LOG ("once%.17g ", t); }
@event early (i++, t < 0.1) { LOG ("x"); }

int main (void)
{
  int rc = events_run ();

  printf ("%s\nout %d bad %d dt %.3g %.3g rc %d t %.17g\n", log_, nout, bad,
          dtmin, dtmax, rc, events_t);
  return 0;
}
