/* The repartition's collectives on real MPI ranks, before and after
 * 2026-11-24: plain MPI, no libb-amr. Build and run from libb-amr:
 *   mpicc -O2 -std=gnu11 tests/bench_table.c -o bench_table
 *   mpirun -np P ./bench_table
 * Each case, timed over many calls with a barrier after each, the slowest
 * rank's mean in microseconds, the barrier included (its time alone is
 * printed first: on one machine the calls cost about as much as it):
 *   small, before   the weight check, 2 MPI_Exscan, 2 MPI_Allreduce
 *   small, now      1 MPI_Allreduce of 3 doubles, 1 MPI_Exscan of 2
 *   table, before   MPI_Allreduce of P entries of 24 B, a user op
 *   table, now      MPI_Allgather of an int, MPI_Allgatherv of P entries:
 *                   one each (spread), or all from rank 0 (concentrated)
 * The tables of both ways are checked equal. */
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>

typedef struct { uint64_t key; double load; long count; } Entry;

static void merge (void * in, void * inout, int * len, MPI_Datatype * type)
{
  (void) type;
  const Entry * a = in;
  Entry * b = inout;
  for (int k = 0; k < *len; k++)
  {
    b [k].key = a [k].key < b [k].key ? a [k].key : b [k].key;
    b [k].load += a [k].load;
    b [k].count += a [k].count;
  }
}

static int rank, npe;

/* the time of one call of f and a barrier, the slowest rank's mean */
static double timed (void (* f) (void), int reps)
{
  for (int k = 0; k < 20; k++)
    f (), MPI_Barrier (MPI_COMM_WORLD);
  double t = MPI_Wtime ();
  for (int k = 0; k < reps; k++)
    f (), MPI_Barrier (MPI_COMM_WORLD);
  t = (MPI_Wtime () - t) / reps;
  double most;
  MPI_Allreduce (&t, &most, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  return most * 1e6;
}

static void nothing (void) {}

static void small_before (void)
{
  int bad = 0, any;
  long cnt = 100, before, total;
  double w = 100, wbefore, wtotal;
  MPI_Allreduce (&bad, &any, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
  MPI_Exscan (&cnt, &before, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
  MPI_Exscan (&w, &wbefore, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce (&cnt, &total, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce (&w, &wtotal, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
}

static void small_now (void)
{
  double s [3] = { 100, 100, 0 }, tot [3], pre [2];
  MPI_Allreduce (s, tot, 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Exscan (s, pre, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
}

/* every rank's leaves are 100, rank r's new piece starts at its own first
 * leaf (spread) or every new piece starts in rank 0's (concentrated) */
static int concentrated;
static Entry * mine, * all;
static MPI_Datatype type;
static MPI_Op op;
static int * lasts, * count, * at;

static uint64_t key_of (int s) { return (uint64_t) s * 1000; }

static void table_before (void)
{
  for (int s = 0; s < npe; s++)
    mine [s] = (Entry) { UINT64_MAX, 0, 0 };
  if (concentrated ? rank == 0 : 1)
    for (int s = concentrated ? 0 : rank; s <= (concentrated ? npe - 1 : rank);
         s++)
      mine [s] = (Entry) { key_of (s), 100, 100 };
  MPI_Allreduce (mine, all, npe, type, op, MPI_COMM_WORLD);
}

static void table_now (void)
{
  int last = concentrated ? npe - 1 : rank, n = 0, prev = -1, given = 0;
  MPI_Allgather (&last, 1, MPI_INT, lasts, 1, MPI_INT, MPI_COMM_WORLD);
  for (int r = 0; r < npe; r++)
  {
    int c = lasts [r] > prev ? lasts [r] - prev : 0;
    if (r == rank)
      given = prev + 1, n = c;
    count [r] = c * (int) sizeof (Entry);
    at [r] = (r ? at [r - 1] + count [r - 1] : 0);
    prev = lasts [r] > prev ? lasts [r] : prev;
  }
  for (int k = 0; k < n; k++)
    mine [k] = (Entry) { key_of (given + k), 100.0 * (given + k),
                         100L * (given + k) };
  MPI_Allgatherv (mine, count [rank], MPI_BYTE, all, count, at, MPI_BYTE,
                  MPI_COMM_WORLD);
}

int main (int argc, char ** argv)
{
  MPI_Init (&argc, &argv);
  MPI_Comm_rank (MPI_COMM_WORLD, &rank);
  MPI_Comm_size (MPI_COMM_WORLD, &npe);
  mine = malloc ((size_t) npe * sizeof (Entry));
  all = malloc ((size_t) npe * sizeof (Entry));
  lasts = malloc ((size_t) npe * sizeof (int));
  count = malloc ((size_t) npe * sizeof (int));
  at = malloc ((size_t) npe * sizeof (int));
  MPI_Type_contiguous ((int) sizeof (Entry), MPI_BYTE, &type);
  MPI_Type_commit (&type);
  MPI_Op_create (merge, 1, &op);
  int reps = argc > 1 ? atoi (argv [1]) : 20000;

  /* the two ways give the same keys (concentrated and spread) */
  int same = 1;
  for (concentrated = 0; concentrated < 2; concentrated++)
  {
    table_before ();
    uint64_t * k0 = malloc ((size_t) npe * sizeof (uint64_t));
    for (int s = 0; s < npe; s++)
      k0 [s] = all [s].key;
    table_now ();
    for (int s = 0; s < npe; s++)
      same &= all [s].key == k0 [s];
    free (k0);
  }

  double b = timed (nothing, reps);
  double sb = timed (small_before, reps), sn = timed (small_now, reps);
  concentrated = 0;
  double tb = timed (table_before, reps), tn = timed (table_now, reps);
  concentrated = 1;
  double tc = timed (table_now, reps);
  if (rank == 0)
    printf ("P %d: barrier %.2f us; small before %.2f, now %.2f; table "
            "before %.2f, now %.2f (spread), %.2f (concentrated); same keys "
            "%s\n", npe, b, sb, sn, tb, tn, tc, same ? "yes" : "NO");
  MPI_Op_free (&op);
  MPI_Type_free (&type);
  free (mine), free (all), free (lasts), free (count), free (at);
  MPI_Finalize ();
  return !same;
}
