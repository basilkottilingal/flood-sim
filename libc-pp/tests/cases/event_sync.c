/* The sync, on one process, with a sync of its own that logs each word:
 * sent at the end of each step and read at the end of the next; a local
 * failure seen one step late; the time limit given to each start; stop ()
 * before end. MODE: 1 a failure at step 2; 2 a failure in init; 3 a word
 * read back changed (out of sync): abort (); 4 CWRAP_EVENTS_TIMEOUT. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <cwrap/event.h>

static char log_ [8192];
static int n_, mode;
#define LOG(...) (n_ += snprintf (log_ + n_, sizeof (log_) - n_, __VA_ARGS__))

static double word [16];

static void start (const double * v, int n, double timeout)
{
  memcpy (word, v, n * sizeof (double));
  LOG ("[s %g %g %g %g %g %g /%g] ", v [0], v [1], v [2], v [3], v [4], v [5],
       timeout);
}

static void finish (double * v, int n)
{
  memcpy (v, word, n * sizeof (double));
  if (mode == 3 && v [2] == 2)
    v [5] = -0.4;                   /* another process at t = 0.4 */
  LOG ("[f] ");
}

static void stop (void) { LOG ("[stop] "); }
static void abort_ (int code)
{
  printf ("%s\nabort %d\n", log_, code);
  exit (0);
}
static int rank (void) { return 0; }

@event init
{
  mode = atoi (getenv ("MODE") ? getenv ("MODE") : "0");
  return mode == 2 ? 6 : 0;
}
@event step (i++, i <= 3)
{
  events_dtmax (0.25);
  LOG ("%d ", i);
  return mode == 1 && i == 2 ? 3 : 0;
}
@event end { LOG ("end%d", event_status ()); }

int main (void)
{
  static const cwrap_event_sync s = { start, finish, stop, abort_, rank };
  int rc;

  events_set_sync (&s);
  events_set_sync_timeout (60);
  rc = events_run ();
  printf ("%s\nrc %d i %d\n", log_, rc, events_i);
  return 0;
}
