/* Events on 3 MPI processes, with <cwrap/event-mpi.h>. MODE: 0 a run in
 * step; 1 process 1 fails at step 2, every process stops after step 3
 * (the sync one step late) and end sees 3; 2 process 1 gives a dt limit of
 * its own: out of sync, MPI_Abort; 3 process 1 is stuck at step 2: the
 * others wait CWRAP_EVENTS_TIMEOUT s, then MPI_Abort; 4 process 2 fails in
 * init: no step runs. Each process prints its log, once. */
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <cwrap/event.h>
#include <cwrap/event-mpi.h>

static char log_ [4096];
static int n_, mode, rank;
#define LOG(...) (n_ += snprintf (log_ + n_, sizeof (log_) - n_, __VA_ARGS__))

@event init
{
  mode = atoi (getenv ("MODE") ? getenv ("MODE") : "0");
  MPI_Comm_rank (MPI_COMM_WORLD, &rank);
  return mode == 4 && rank == 2 ? 4 : 0;
}
@event limit (i++, i <= 5)
{
  double dt = mode == 2 && rank == 1 ? 0.1 : 0.2;
  double all;

  /* what a foreach () reduction does: a value, the same on every process */
  MPI_Allreduce (&dt, &all, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  events_dtmax (mode == 2 ? dt : all);
  LOG ("%d ", i);
  if (mode == 3 && rank == 1 && i == 2)
    sleep (30);
  return mode == 1 && rank == 1 && i == 2 ? 3 : 0;
}
@event end { LOG ("end%d", event_status ()); }

int main (int argc, char ** argv)
{
  int rc;

  MPI_Init (&argc, &argv);
  rc = events_run ();
  printf ("p%d: %s rc %d t %g\n", rank, log_, rc, events_t);
  MPI_Finalize ();
  return 0;
}
