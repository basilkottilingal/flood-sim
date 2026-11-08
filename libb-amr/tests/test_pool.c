/* pool.c: the hot, muzzy and cold free lists, with page faults counted, in
 * either way of giving memory back (_BAMR_MADV_DONTNEED_RECOMMENDED_) */
#include <stdio.h>
#include <string.h>
#include <sys/resource.h>
#include <time.h>
#include "pool.h"
static int fails;
#define CHECK(c) do { if (!(c)) { fails++; printf ("FAIL line %d: %s\n", __LINE__, #c); } } while (0)
static long faults (void) { struct rusage u; getrusage (RUSAGE_SELF, &u); return u.ru_minflt; }
static double now (void) { struct timespec t; clock_gettime (CLOCK_MONOTONIC, &t); return t.tv_sec + 1e-9 * t.tv_nsec; }
enum { N = 64, SIZE = 20 * 20 * 20 * 8, PAGE = 4096 };   /* 3D blocks of 16 + 2 ghosts */

/* free all, then allocate them again with f, writing each in full; the
 * time of the frees (madvise, if any) and of the allocations */
static long cycle (Pool * p, void ** o, void * (*f) (Pool *), double * fsecs, double * secs)
{
  double t0 = now ();
  for (int i = 0; i < N; i++) pool_free (p, o [i]);
  *fsecs = now () - t0;
  long f0 = faults (); t0 = now ();
  for (int i = 0; i < N; i++) { o [i] = f (p); memset (o [i], 1, SIZE); }
  *secs = now () - t0;
  return faults () - f0;
}

static int all (const void * q, size_t from, size_t to, unsigned char v)
{
  for (size_t i = from; i < to; i++)
    if (((const unsigned char *) q) [i] != v)
      return 0;
  return 1;
}

int main (void)
{
  void * o [N];
  double fsecs, secs;
  int dontneed = _BAMR_MADV_DONTNEED_RECOMMENDED_;
  printf ("giving memory back with %s\n", dontneed ? "MADV_DONTNEED" : "MADV_FREE");

  /* hot list (the default): reuse costs no fault */
  Pool * p = pool_new (SIZE);
  for (int i = 0; i < N; i++) { o [i] = pool_alloc (p); memset (o [i], 1, SIZE); }
  long hot = cycle (p, o, pool_alloc_uninit, &fsecs, &secs);
  printf ("hot list, uninit:  %5ld faults, free %.1f us, alloc %.1f us per block\n",
          hot, 1e6 * fsecs / N, 1e6 * secs / N);
  long hotz = cycle (p, o, pool_alloc, &fsecs, &secs);
  printf ("hot list, zeroed:  %5ld faults, free %.1f us, alloc %.1f us per block\n",
          hotz, 1e6 * fsecs / N, 1e6 * secs / N);
  CHECK (hot < 8 && hotz < 8 && pool_used (p) == N);

  /* uninit keeps the old content: no clearing was done */
  pool_free (p, o [0]);
  void * q = pool_alloc_uninit (p);
  CHECK (q == o [0] && ((unsigned char *) q) [100] == 1);
  /* zeroed really is zero, from the hot list */
  pool_free (p, q);
  q = pool_alloc (p);
  CHECK (q == o [0] && all (q, 0, SIZE, 0));
  o [0] = q;

  /* no hot list: objects freed go cold (MADV_DONTNEED), or muzzy
   * (MADV_FREE) */
  pool_set_hot_max (p, 0);
  long none = cycle (p, o, pool_alloc_uninit, &fsecs, &secs);
  printf ("no hot list:       %5ld faults, free %.1f us, alloc %.1f us per block\n",
          none, 1e6 * fsecs / N, 1e6 * secs / N);
  if (dontneed)
    CHECK (none >= N * (SIZE / PAGE));          /* every page, every reuse */
  else
    /* no fault while memory is plentiful: the system took no page */
    CHECK (none < N * (SIZE / PAGE) / 2);

  /* a freed object read as zero, whichever list it is on */
  pool_free (p, o [0]);
  CHECK (dontneed ? pool_cold (p) == 1 && pool_muzzy (p) == 0
                  : pool_muzzy (p) == 1 && pool_cold (p) == 0);
  CHECK (pool_hot (p) == 0);
  q = pool_alloc (p);
  CHECK (q == o [0] && all (q, 0, SIZE, 0));
  o [0] = q;

  /* trim: hot beyond keep, and every muzzy object, go cold */
  pool_set_hot_max (p, POOL_HOT_ALL);
  for (int i = 0; i < N; i++) pool_free (p, o [i]);
  CHECK (pool_hot (p) == N && pool_cold (p) == 0 && pool_used (p) == 0);
  pool_trim (p, 10);
  CHECK (pool_hot (p) == 10 && pool_cold (p) == N - 10 && pool_muzzy (p) == 0);
  /* zeroed takes cold first, uninit hot first */
  void * a = pool_alloc (p), * b = pool_alloc_uninit (p);
  CHECK (pool_cold (p) == N - 11 && pool_hot (p) == 9);
  (void) a; (void) b;
  pool_destroy (p);

  /* MADV_FREE: the muzzy list */
  if (!dontneed) {
    p = pool_new (SIZE);
    for (int i = 0; i < 3; i++) { o [i] = pool_alloc (p); memset (o [i], 7, SIZE); }
    /* over the limit: muzzy, and each page but the link's old or zero */
    pool_set_hot_max (p, 1);
    pool_free (p, o [0]);
    pool_free (p, o [1]);
    CHECK (pool_hot (p) == 1 && pool_muzzy (p) == 1 && pool_cold (p) == 0);
    pool_set_hot_max (p, 0);                    /* the hot one goes muzzy */
    CHECK (pool_hot (p) == 0 && pool_muzzy (p) == 2);
    q = pool_alloc_uninit (p);
    int pages_ok = all (q, sizeof (void *), PAGE, 7);   /* never advised */
    for (size_t pg = PAGE; pg < SIZE; pg += PAGE) {
      size_t end = pg + PAGE < SIZE ? pg + PAGE : SIZE;
      pages_ok &= all (q, pg, end, 7) || all (q, pg, end, 0);
    }
    CHECK (pages_ok);
    /* the order: uninit hot, muzzy, cold; zeroed cold, hot, muzzy */
    pool_set_hot_max (p, POOL_HOT_ALL);
    pool_free (p, q);                           /* hot again */
    pool_trim (p, 1);                           /* the muzzy one: cold */
    CHECK (pool_hot (p) == 1 && pool_cold (p) == 1 && pool_muzzy (p) == 0);
    pool_set_hot_max (p, 1);
    pool_free (p, o [2]);                       /* over the limit: muzzy */
    CHECK (pool_hot (p) == 1 && pool_muzzy (p) == 1 && pool_cold (p) == 1);
    void * h = pool_alloc_uninit (p);
    CHECK (pool_hot (p) == 0 && pool_muzzy (p) == 1);
    void * m = pool_alloc_uninit (p);
    CHECK (pool_muzzy (p) == 0 && pool_cold (p) == 1);
    pool_free (p, h);                           /* hot */
    pool_free (p, m);                           /* over the limit: muzzy */
    void * z = pool_alloc (p);
    CHECK (pool_cold (p) == 0 && all (z, 0, SIZE, 0));
    z = pool_alloc (p);
    CHECK (pool_hot (p) == 0 && pool_muzzy (p) == 1 && all (z, 0, SIZE, 0));
    z = pool_alloc (p);
    CHECK (pool_muzzy (p) == 0 && all (z, 0, SIZE, 0));   /* cleared */
    pool_destroy (p);
  }

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
  pool_set_hot_max (s, 0);
  pool_free (s, y);
  pool_trim (s, 0);
  CHECK (pool_cold (s) == 0 && pool_muzzy (s) == 0 && pool_hot (s) == 1);
  pool_destroy (s);
  printf ("%s\n", fails ? "FAILED" : "pool ok");
  return fails != 0;
}
