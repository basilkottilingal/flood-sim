/* profile.h on one process: regions of known lengths (busy waits on the
 * clock), nested, left by return, break and goto, recursion past the
 * depth the profiler keeps, and threads; the paths, calls and times
 * against the lengths; then the outputs written (the run script checks
 * the files). Built with -DPROFILE_ENABLE, and without it (nothing timed,
 * no path found). */
#include <pthread.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include "profile.h"

static int checks, fails;
#define CHECK(c) do { checks++; if (!(c)) { fails++; \
  printf ("FAIL line %d: %s\n", __LINE__, #c); } } while (0)

static double clock_now (void)
{
  struct timespec t;
  clock_gettime (CLOCK_MONOTONIC, &t);
  return (double) t.tv_sec + 1e-9 * (double) t.tv_nsec;
}

/* busy for s seconds: a sleep could last longer than asked */
static void spin (double s)
{
  double end = clock_now () + s;
  while (clock_now () < end)
    ;
}

static void inner (void)
{
  PROFILE ();
  spin (0.001);
}

static void outer (void)
{
  PROFILE ();
  spin (0.002);
  inner ();
  inner ();
}

/* left by an early return, by break and by goto */
static int early (int k)
{
  PROFILE ();
  if (k % 2)
    return 1;
  for (int i = 0; ; i++)
  {
    PROFILE_SCOPE ("loop body");
    if (i == 3)
      break;
  }
  {
    PROFILE_SCOPE ("block");
    goto out;
  }
out:
  return 0;
}

static int depth_reached;
static void recurse (int n)
{
  PROFILE ();
  if (n > depth_reached)
    depth_reached = n;
  if (n < 80)
    recurse (n + 1);
}

static void * worker (void * arg)
{
  (void) arg;
  for (int k = 0; k < 5; k++)
    outer ();
  return NULL;
}

#ifdef PROFILE_ENABLE
/* |a - b| within rel of b, or within 0.5 ms (the clock, the scheduler) */
static int near (double a, double b, double rel)
{
  double d = a > b ? a - b : b - a;
  return d <= rel * b || d <= 5e-4;
}
#endif

int main (void)
{
  for (int k = 0; k < 10; k++)
    outer ();
  for (int k = 0; k < 4; k++)
    early (k);
  recurse (0);
  pthread_t t [3];
  for (int i = 0; i < 3; i++)
    pthread_create (&t [i], NULL, worker, NULL);
  for (int i = 0; i < 3; i++)
    pthread_join (t [i], NULL);
  profile_snapshot ();

  double calls, incl, self;
#ifdef PROFILE_ENABLE
  /* outer: 10 calls here, 5 on each worker (the same path), times the
   * longest thread's (10 here: 40 ms in all, 20 of them its own) */
  CHECK (profile_query ("outer", &calls, &incl, &self, NULL));
  CHECK (calls == 25);
  CHECK (near (incl, 0.040, 0.15) && near (self, 0.020, 0.15));
  CHECK (profile_query ("outer;inner", &calls, &incl, &self, NULL));
  CHECK (calls == 50 && near (incl, 0.020, 0.15) && near (self, incl, 0.01));
  CHECK (!profile_query ("inner", NULL, NULL, NULL, NULL));     /* never a root */
  /* every way out of a block ends its region */
  CHECK (profile_query ("early", &calls, NULL, NULL, NULL) && calls == 4);
  CHECK (profile_query ("early;loop body", &calls, NULL, NULL, NULL)
         && calls == 8);
  CHECK (profile_query ("early;block", &calls, NULL, NULL, NULL) && calls == 2);
  /* recursion: a path a level, 64 of them kept, the rest not timed */
  char path [4096] = "recurse";
  for (int d = 1; d < 64; d++)
    strcat (path, ";recurse");
  CHECK (depth_reached == 80);
  CHECK (profile_query (path, &calls, NULL, NULL, NULL) && calls == 1);
  strcat (path, ";recurse");
  CHECK (!profile_query (path, NULL, NULL, NULL, NULL));
#else
  CHECK (!profile_query ("outer", &calls, &incl, &self, NULL));
#endif
  profile_finish ();
  printf ("%d checks, %d failed\n", checks, fails);
  return fails != 0;
}
