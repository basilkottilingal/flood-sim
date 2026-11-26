/*
.. The sync of events over MPI processes. Include it in one file of the
.. program, after <mpi.h>, compiled with MPI; it gives events_run () the
.. sync of include/cwrap/event.h, at start-up:
..
..   - at the end of each step, a status word of 6 doubles goes out by
..     MPI_Iallreduce (MPI_MAX), on a communicator of its own (MPI_COMM_WORLD
..     duplicated, so that it never meets the program's collectives), and
..     is waited for (MPI_Wait) at the end of the next step: the time it
..     takes is hidden behind that step;
..   - each step has at most CWRAP_EVENTS_TIMEOUT seconds (600 by default, 0
..     for no limit; events_set_sync_timeout () in the program), from the
..     end of the step before: a timer (setitimer, SIGALRM) that fires
..     wherever the process is, in an event, in a collective waiting for a
..     process stuck or gone, or in the sync. It writes why, then ends the
..     process (_exit (1): MPI_Abort is not safe in a signal handler), and
..     the launcher ends the others. The time limit is off from the start
..     of end. A program that uses SIGALRM or ITIMER_REAL itself sets 0;
..   - abort () is MPI_Abort (MPI_COMM_WORLD, code).
..
.. MPI must be started (MPI_Init, or libb-amr's tree ()) before
.. events_run (). On one process, no message is sent.
*/
#ifndef CWRAP_EVENT_MPI_H
#define CWRAP_EVENT_MPI_H

#include <mpi.h>
#include <signal.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <unistd.h>
#include <cwrap/event.h>

  enum { CWRAP_EVENT_SYNC_MAX = 16 };

  static MPI_Comm __cwrap_event_comm = MPI_COMM_NULL;
  static MPI_Request __cwrap_event_req = MPI_REQUEST_NULL;
  static double __cwrap_event_out [CWRAP_EVENT_SYNC_MAX],
    __cwrap_event_in [CWRAP_EVENT_SYNC_MAX];
  static int __cwrap_event_npe, __cwrap_event_me;
  static double __cwrap_event_limit;

  static int __cwrap_event_rank (void)
  {
    int r = 0, on = 0;

    MPI_Initialized (&on);
    if (on)
      MPI_Comm_rank (MPI_COMM_WORLD, &r);
    return r;
  }

  /* the step ran out of time: say so, as a signal handler may, and end */
  static void __cwrap_event_alarm (int sig)
  {
    char m [640];
    int n = snprintf (m, sizeof (m), "cwrap: events (): process %d: step %d "
                      "(t = %g) has run more than %g s "
                      "(CWRAP_EVENTS_TIMEOUT). Possible reasons: (a) a "
                      "process crashed or exited, and the others wait for "
                      "it in a collective; (b) a process is stuck in an "
                      "event; (c) the processes are out of sync, one gone "
                      "on to a collective the others never reach. "
                      "Aborting\n", __cwrap_event_me, events_i, events_t,
                      __cwrap_event_limit);

    (void) sig;
    if (n > 0 && write (2, m, n < (int) sizeof (m) ? n : (int) sizeof (m))) {}
    _exit (1);
  }

  static void __cwrap_event_arm (double s)
  {
    struct itimerval it;

    memset (&it, 0, sizeof (it));
    it.it_value.tv_sec = (time_t) s;
    it.it_value.tv_usec = (suseconds_t) ((s - (double) (time_t) s) * 1e6);
    setitimer (ITIMER_REAL, &it, NULL);
  }

  static void __cwrap_event_start (const double * v, int n, double timeout)
  {
    if (__cwrap_event_comm == MPI_COMM_NULL) {
      int on = 0;

      MPI_Initialized (&on);
      if (!on) {
        fprintf (stderr, "cwrap: events (): <cwrap/event-mpi.h> needs MPI "
                 "started (MPI_Init (), or tree ()) before events_run ()\n");
        exit (1);
      }
      MPI_Comm_dup (MPI_COMM_WORLD, &__cwrap_event_comm);
      MPI_Comm_size (__cwrap_event_comm, &__cwrap_event_npe);
      MPI_Comm_rank (__cwrap_event_comm, &__cwrap_event_me);
      if (timeout > 0) {
        struct sigaction sa;

        memset (&sa, 0, sizeof (sa));
        sa.sa_handler = __cwrap_event_alarm;
        sigaction (SIGALRM, &sa, NULL);
      }
    }
    __cwrap_event_limit = timeout;
    if (timeout > 0)
      __cwrap_event_arm (timeout);  /* the next step's time, from now */
    memcpy (__cwrap_event_out, v, n * sizeof (double));
    if (__cwrap_event_npe > 1)
      MPI_Iallreduce (__cwrap_event_out, __cwrap_event_in, n, MPI_DOUBLE,
                      MPI_MAX, __cwrap_event_comm, &__cwrap_event_req);
    else
      memcpy (__cwrap_event_in, v, n * sizeof (double));
  }

  static void __cwrap_event_finish (double * v, int n)
  {
    if (__cwrap_event_npe > 1)
      MPI_Wait (&__cwrap_event_req, MPI_STATUS_IGNORE);
    memcpy (v, __cwrap_event_in, n * sizeof (double));
  }

  static void __cwrap_event_stop (void)
  {
    if (__cwrap_event_limit > 0)
      __cwrap_event_arm (0);        /* off */
  }

  static void __cwrap_event_abort (int code)
  {
    MPI_Abort (MPI_COMM_WORLD, code);
  }

  __attribute__ ((constructor))
  static void __cwrap_event_mpi (void)
  {
    static const cwrap_event_sync s = {
      __cwrap_event_start, __cwrap_event_finish, __cwrap_event_stop,
      __cwrap_event_abort, __cwrap_event_rank,
    };

    events_set_sync (&s);
  }

#endif
