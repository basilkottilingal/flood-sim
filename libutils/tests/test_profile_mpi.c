/* profile.h on MPI processes (real MPI, MPI-3): rank r busy for (r + 1)
 * ms a step, 20 steps, a snapshot each; rank 2 enters a region of its own
 * from step 5 on, a path new after five snapshots; at step 3 every rank
 * enters a region named after it and a shared one, many paths new at
 * once, some on every process; a snapshot begun while rank 0 is 100 ms
 * late returns at once on the others. On rank 0, after the end: the
 * paths, calls, mean times and imbalance against the work. */
#include <mpi.h>
#include <stdio.h>
#include <time.h>
#include "profile.h"

static int checks, fails, rank, npe;
#define CHECK(c) do { checks++; if (!(c)) { fails++; \
  printf ("rank %d: FAIL line %d: %s\n", rank, __LINE__, #c); } } while (0)

static double clock_now (void)
{
  struct timespec t;
  clock_gettime (CLOCK_MONOTONIC, &t);
  return (double) t.tv_sec + 1e-9 * (double) t.tv_nsec;
}

static void spin (double s)
{
  double end = clock_now () + s;
  while (clock_now () < end)
    ;
}

static void only2 (void)
{
  PROFILE ();
  spin (0.0005);
}

static void step (int n)
{
  PROFILE ();
  spin (0.001 * (rank + 1));
  if (rank == 2 && n >= 5)
    only2 ();
  if (n == 3)
  {
    char name [32];
    snprintf (name, sizeof name, "rank %d", rank);
    { PROFILE_SCOPE ("shared"); }
    { PROFILE_SCOPE (name); }
  }
}

static int near (double a, double b, double rel)
{
  double d = a > b ? a - b : b - a;
  return d <= rel * b || d <= 5e-4;
}

int main (int argc, char ** argv)
{
  MPI_Init (&argc, &argv);
  MPI_Comm_rank (MPI_COMM_WORLD, &rank);
  MPI_Comm_size (MPI_COMM_WORLD, &npe);
  enum { STEPS = 20 };
  step (-1);                          /* a path, before the first snapshot */

  /* a snapshot while rank 0 is late: no wait on the others */
  if (rank == 0)
    spin (0.1);
  double t = clock_now ();
  profile_snapshot ();
  t = clock_now () - t;
  if (rank != 0)
    CHECK (t < 0.02);

  for (int n = 0; n < STEPS; n++)
  {
    step (n);
    profile_snapshot ();
  }
  profile_finish ();

  if (rank == 0)
  {
    /* step: STEPS + 1 calls a process; rank r (r + 1) ms each: mean
     * (P + 1) / 2 ms */
    double calls, incl, self;
    double mean = 0.001 * (npe + 1) / 2 * (STEPS + 1);
    double imbal, only = npe > 2 ? 15 * 0.0005 : 0;
    CHECK (profile_query ("step", &calls, &incl, &self, &imbal));
    CHECK (calls == STEPS + 1);
    CHECK (near (incl, mean, 0.15));
    /* the slowest: rank P - 1, and rank 2's own region when P is 3 */
    double most = 0.001 * npe * (STEPS + 1) + (npe == 3 ? only : 0);
    CHECK (near (imbal, most / (mean + only / npe), 0.15));
    /* one call on each process, and each process's own: 1 / P */
    CHECK (profile_query ("step;shared", &calls, NULL, NULL, NULL)
           && calls == 1);
    for (int r = 0; r < npe; r++)
    {
      char path [64];
      snprintf (path, sizeof path, "step;rank %d", r);
      CHECK (profile_query (path, &calls, NULL, NULL, NULL)
             && near (calls, 1. / npe, 1e-9));
    }
    if (npe > 2)
    {
      /* only rank 2, from step 5: 15 calls, the mean over P processes */
      CHECK (profile_query ("step;only2", &calls, &incl, NULL, NULL));
      CHECK (near (calls, 15. / npe, 1e-9));
      CHECK (near (incl, 15 * 0.0005 / npe, 0.25));
    }
    else
      CHECK (!profile_query ("step;only2", NULL, NULL, NULL, NULL));
    printf ("%d processes: %d checks, %d failed (rank 0)\n", npe, checks,
            fails);
  }
  int f = fails, all;
  MPI_Allreduce (&f, &all, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  MPI_Finalize ();
  return all != 0;
}
