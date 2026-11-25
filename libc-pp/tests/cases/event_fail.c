/* A failure stops the run, the slot's later bodies and steps skipped; end
 * runs and sees the code. events_stop (), event_cancel (), and the reduce
 * hook called after every body. FAIL_AT=-1: no failure. */
#include <stdio.h>
#include <stdlib.h>
#include <cwrap/event.h>

static char log_ [4096];
static int n_, nreduce, fail_at;
#define LOG(...) (n_ += snprintf (log_ + n_, sizeof (log_) - n_, __VA_ARGS__))

/* as an MPI_Allreduce (MAX) would: here, one process */
static int reduce (int status)
{
  nreduce++;
  return status;
}

@event init { fail_at = atoi (getenv ("FAIL_AT") ? getenv ("FAIL_AT") : "-1"); }
@event a (i++, i <= 6) { LOG ("a%d ", i); return i == fail_at ? 3 : 0; }
@event a { LOG ("A%d ", i); }
@event b (i++, next.a) { LOG ("b%d ", i); if (i == 4) events_stop (); }
@event c (i++, next.b)
{
  LOG ("c%d ", i);
  if (i == 1)
    event_cancel ();
}
@event c { LOG ("C%d ", i); }
@event end { LOG ("E%d ", event_status ()); return 7; }
@event end (last) { LOG ("F%d", event_status ()); }

int main (void)
{
  int rc;

  events_set_reduce (reduce);
  rc = events_run ();
  printf ("%s\nrc %d i %d reduce %d\n", log_, rc, events_i, nreduce);
  return 0;
}
