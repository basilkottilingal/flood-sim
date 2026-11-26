/* event_fail (code): the body goes on to its end (its collectives, here a
 * log), then the run stops as if it had returned the code; the first code
 * is kept over the one returned; event_fail (0) and a call outside a body
 * do nothing; the reduce hook sees the code. */
#include <stdio.h>
#include <cwrap/event.h>

static char log_ [4096];
static int n_;
#define LOG(...) (n_ += snprintf (log_ + n_, sizeof (log_) - n_, __VA_ARGS__))

static int reduce (int status)
{
  LOG ("r%d ", status);
  return status;
}

@event a (i++, i <= 5)
{
  event_fail (0);
  if (i == 2)
    event_fail (4);
  LOG ("a%d ", i);                  /* the collective the others enter */
  return i == 2 ? 5 : 0;
}
@event a { LOG ("A%d ", i); }
@event end { LOG ("E%d", event_status ()); }

int main (void)
{
  int rc;

  event_fail (9);                   /* no body running: nothing */
  events_set_reduce (reduce);
  rc = events_run ();
  printf ("%s\nrc %d i %d\n", log_, rc, events_i);
  return 0;
}
