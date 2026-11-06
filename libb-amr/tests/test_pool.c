/* pool.c: the hot and cold free lists, with page faults counted */
#include <stdio.h>
#include <string.h>
#include <sys/resource.h>
#include <time.h>
#include "pool.h"
static int fails;
#define CHECK(c) do { if (!(c)) { fails++; printf ("FAIL line %d: %s\n", __LINE__, #c); } } while (0)
static long faults (void) { struct rusage u; getrusage (RUSAGE_SELF, &u); return u.ru_minflt; }
static double now (void) { struct timespec t; clock_gettime (CLOCK_MONOTONIC, &t); return t.tv_sec + 1e-9 * t.tv_nsec; }
enum { N = 64, SIZE = 20 * 20 * 20 * 8 };     /* 3D blocks of 16 + 2 ghosts */

/* free all, then allocate them again with f, writing each in full */
static long cycle (Pool * p, void ** o, void * (*f) (Pool *), double * secs)
{
  for (int i = 0; i < N; i++) pool_free (p, o [i]);
  long f0 = faults (); double t0 = now ();
  for (int i = 0; i < N; i++) { o [i] = f (p); memset (o [i], 1, SIZE); }
  *secs = now () - t0;
  return faults () - f0;
}

int main (void)
{
  void * o [N];
  double secs;
  /* hot list (the default): reuse costs no fault */
  Pool * p = pool_new (SIZE);
  for (int i = 0; i < N; i++) { o [i] = pool_alloc (p); memset (o [i], 1, SIZE); }
  long hot = cycle (p, o, pool_alloc_uninit, &secs);
  printf ("hot list, uninit:  %5ld faults, %.1f us per block\n", hot, 1e6 * secs / N);
  long hotz = cycle (p, o, pool_alloc, &secs);
  printf ("hot list, zeroed:  %5ld faults, %.1f us per block\n", hotz, 1e6 * secs / N);
  CHECK (hot < 8 && hotz < 8 && pool_used (p) == N);

  /* uninit keeps the old content: no clearing was done */
  pool_free (p, o [0]);
  void * q = pool_alloc_uninit (p);
  CHECK (q == o [0] && ((unsigned char *) q) [100] == 1);
  /* zeroed really is zero, from the hot list */
  pool_free (p, q);
  q = pool_alloc (p);
  int zero = 1;
  for (int i = 0; i < SIZE; i++) zero &= ((unsigned char *) q) [i] == 0;
  CHECK (q == o [0] && zero);
  o [0] = q;

  /* no hot list: as before, every reuse faults on every page */
  pool_set_hot_max (p, 0);
  long cold = cycle (p, o, pool_alloc_uninit, &secs);
  printf ("no hot list:       %5ld faults, %.1f us per block\n", cold, 1e6 * secs / N);
  CHECK (cold >= N * (SIZE / 4096));
  /* cold objects read as zero */
  pool_free (p, o [0]);
  CHECK (pool_cold (p) == 1 && pool_hot (p) == 0);
  q = pool_alloc (p);
  zero = 1;
  for (int i = 0; i < SIZE; i++) zero &= ((unsigned char *) q) [i] == 0;
  CHECK (zero);
  o [0] = q;

  /* trim: hot beyond keep goes cold */
  pool_set_hot_max (p, POOL_HOT_ALL);
  for (int i = 0; i < N; i++) pool_free (p, o [i]);
  CHECK (pool_hot (p) == N && pool_cold (p) == 0 && pool_used (p) == 0);
  pool_trim (p, 10);
  CHECK (pool_hot (p) == 10 && pool_cold (p) == N - 10);
  /* zeroed takes cold first, uninit hot first */
  void * a = pool_alloc (p), * b = pool_alloc_uninit (p);
  CHECK (pool_cold (p) == N - 11 && pool_hot (p) == 9);
  (void) a; (void) b;
  pool_destroy (p);

  /* small objects: hot only, zeroed or not */
  Pool * s = pool_new (24);
  int * x = pool_alloc (s);
  x [0] = 7;
  pool_free (s, x);
  int * y = pool_alloc_uninit (s);
  CHECK (y == x);
  pool_free (s, y);
  y = pool_alloc (s);
  CHECK (y == x && y [0] == 0 && y [5] == 0);
  pool_trim (s, 0);
  CHECK (pool_cold (s) == 0);
  pool_destroy (s);
  printf ("%s\n", fails ? "FAILED" : "pool ok");
  return fails != 0;
}
