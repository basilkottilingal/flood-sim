/* Long steps that land on times exactly: t + (tn - t) is not tn for t =
 * 0.03, tn = 0.29, so tn must be taken as it is; a bound t <= 0.97 that no
 * event's time falls on. Slots first and last against their file's order. */
#include <stdio.h>
#include <cwrap/event.h>

static char log_ [4096];
static int n_;
#define LOG(...) (n_ += snprintf (log_ + n_, sizeof (log_) - n_, __VA_ARGS__))

static double a, b, c;

@event tail (i++, last) { LOG ("tail "); }
@event step (i++) { events_dtnext (1); LOG ("s%d ", i); }
@event at3 (t = 0.03) { a = t; LOG ("a "); }
@event at29 (t = 0.29) { b = t; LOG ("b "); }
@event upto (i++, t <= 0.97) { c = t; }
@event head (i++, first) { LOG ("head "); }

int main (void)
{
  int rc = events_run ();

  printf ("%s\n%d %d %d rc %d t %g\n", log_, a == 0.03, b == 0.29, c == 0.97,
          rc, events_t);
  return 0;
}
