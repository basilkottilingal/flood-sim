/*
.. memindex.[ch] against voxel-memindex.[ch]: memory, page faults and the
.. average time of a lookup, for one implementation, pattern and level.
..
..   bench_memindex IMPL PATTERN LEVEL [N]
..     IMPL     mem | voxel
..     PATTERN  dense    every index of the level
..              shell    a sphere (3D) or circle (2D) one index thick
..              river    a meandering band (2D) or tube (3D), N wide
..              random   N indices, uniform
..
.. Points are inserted in Morton order, as a refinement would roughly do.
.. Phases, each timed, with minor page faults and, if the system allows
.. (kernel.perf_event_paranoid <= 2), cache and dTLB load misses:
..   build     insert every point
..   ordered   look every point up, Morton order
..   shuffled  look every point up, random order
..   nbr       the 2*dimension face neighbours of every point, Morton
..             order (hits and misses, out of range too)
..   nbrcur    the same through a cursor (voxel only)
..   uniform   indices drawn uniformly over the level (mostly misses
..             when sparse)
..   free      free every point, random order
.. One line of key=value pairs out; tests/bench_memindex.sh tabulates.
.. Build: gcc -O2 -Ddimension=D -D_GNU_SOURCE -I. tests/bench_memindex.c
..        memindex.c voxel-memindex.c pool.c
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include <time.h>
#include <unistd.h>
#include <sys/resource.h>
#include <sys/ioctl.h>
#include <sys/syscall.h>
#include <linux/perf_event.h>

#include "memindex.h"
#include "voxel-memindex.h"

#define D dimension

typedef struct { int x [D]; } Pt;

#if D == 1
# define P(p)  (p).x [0]
#elif D == 2
# define P(p)  (p).x [0], (p).x [1]
#else
# define P(p)  (p).x [0], (p).x [1], (p).x [2]
#endif

/* f (a, b, ...) after the arguments are expanded, as TREE_CALL in tree.h */
#define CALL(f, ...)  f (__VA_ARGS__)

static int use_voxel, len;
static struct _Memindex mi;
static struct _VMemindex vi;

static inline void * get (Pt p)
{
  if (use_voxel)
    return vmem_get (&vi, P (p));
  /* as tree.h's BAMR_NODE_ENV () does it */
  return mem_allocated (&mi, P (p)) ? CALL (mem_data, &mi, P (p)) : NULL;
}

static inline void assign (Pt p, void * b)
{
  if (use_voxel)
    vmem_assign (&vi, P (p), len, b);
  else
    mem_assign (&mi, P (p), len, b);
}

static inline void release (Pt p)
{
  if (use_voxel)
    vmem_free (&vi, P (p), len, NULL);
  else
    mem_free (&mi, P (p), len, NULL);
}

/* ---- measuring ------------------------------------------------------ */

static double now (void)
{
  struct timespec t;
  clock_gettime (CLOCK_MONOTONIC, &t);
  return t.tv_sec + 1e-9 * t.tv_nsec;
}

static long minflt (void)
{
  struct rusage u;
  getrusage (RUSAGE_SELF, &u);
  return u.ru_minflt;
}

/* kB of a /proc/self/status field */
static long status_kb (const char * field)
{
  char line [256];
  long v = -1;
  size_t n = strlen (field);
  FILE * f = fopen ("/proc/self/status", "r");
  while (f && fgets (line, sizeof line, f))
    if (!strncmp (line, field, n)) {
      v = atol (line + n + 1);
      break;
    }
  if (f)
    fclose (f);
  return v;
}

static int perf_fd [2] = {-1, -1};

static int perf_open (uint32_t type, uint64_t config)
{
  struct perf_event_attr a;
  memset (&a, 0, sizeof a);
  a.size = sizeof a;
  a.type = type;
  a.config = config;
  a.disabled = 1;
  a.exclude_kernel = 1;
  a.exclude_hv = 1;
  return syscall (SYS_perf_event_open, &a, 0, -1, -1, 0);
}

static void perf_init (void)
{
  perf_fd [0] = perf_open (PERF_TYPE_HARDWARE, PERF_COUNT_HW_CACHE_MISSES);
  perf_fd [1] = perf_open (PERF_TYPE_HW_CACHE, PERF_COUNT_HW_CACHE_DTLB |
                           PERF_COUNT_HW_CACHE_OP_READ << 8 |
                           PERF_COUNT_HW_CACHE_RESULT_MISS << 16);
}

typedef struct {
  double t;
  long flt;
} Mark;

static Mark mark (void)
{
  for (int e = 0; e < 2; e++)
    if (perf_fd [e] >= 0) {
      ioctl (perf_fd [e], PERF_EVENT_IOC_RESET, 0);
      ioctl (perf_fd [e], PERF_EVENT_IOC_ENABLE, 0);
    }
  return (Mark) {now (), minflt ()};
}

/* ns per operation, faults, and misses per operation, as key=value */
static void report (const char * name, Mark m, double ops)
{
  double t = now () - m.t;
  long flt = minflt () - m.flt;
  printf (" %s_ns=%.2f %s_flt=%ld", name, 1e9 * t / ops, name, flt);
  const char * ev [2] = {"llc", "tlb"};
  for (int e = 0; e < 2; e++) {
    long long c = -1;
    if (perf_fd [e] >= 0) {
      ioctl (perf_fd [e], PERF_EVENT_IOC_DISABLE, 0);
      if (read (perf_fd [e], &c, sizeof c) != sizeof c)
        c = -1;
    }
    if (c >= 0)
      printf (" %s_%s=%.3f", name, ev [e], c / ops);
    else
      printf (" %s_%s=n/a", name, ev [e]);
  }
}

/* ---- points --------------------------------------------------------- */

static uint64_t rng = 88172645463325252ull;
static uint64_t xorshift (void)
{
  rng ^= rng << 13;
  rng ^= rng >> 7;
  rng ^= rng << 17;
  return rng;
}

static uint64_t morton (Pt p)
{
  uint64_t k = 0;
  for (int b = 0; b < 21; b++)
    for (int d = 0; d < D; d++)
      k |= (uint64_t) ((p.x [d] >> b) & 1) << (b * D + d);
  return k;
}

typedef struct { uint64_t key; Pt p; } Keyed;

static int by_key (const void * a, const void * b)
{
  uint64_t x = ((const Keyed *) a)->key, y = ((const Keyed *) b)->key;
  return x < y ? -1 : x > y;
}

static Keyed * pts;
static size_t npts, maxpts;

static void add (Pt p)
{
  for (int d = 0; d < D; d++)
    if (p.x [d] < 0 || p.x [d] >= len)
      return;
  if (npts == maxpts) {
    maxpts = maxpts ? 2 * maxpts : 1 << 16;
    pts = realloc (pts, maxpts * sizeof (Keyed));
    if (!pts)
      abort ();
  }
  pts [npts++].p = p;
}

static void dense (int d, Pt p)
{
  if (d == D) {
    add (p);
    return;
  }
  for (p.x [d] = 0; p.x [d] < len; p.x [d]++)
    dense (d + 1, p);
}

static void make_points (const char * pattern, long n)
{
  double c = 0.5 * len, R = 0.4 * len;
  Pt p = {{0}};
  if (!strcmp (pattern, "dense"))
    dense (0, p);
  else if (!strcmp (pattern, "random"))
    for (long q = 0; q < n; q++) {
      for (int d = 0; d < D; d++)
        p.x [d] = xorshift () % len;
      add (p);
    }
  else if (!strcmp (pattern, "shell")) {
    /* |x - c| within half an index of R: along the last axis, both sides */
#if D == 1
    p.x [0] = c - R; add (p);
    p.x [0] = c + R; add (p);
#else
    long outer = 1;
    for (int d = 0; d < D - 1; d++)
      outer *= len;
    for (long q = 0; q < outer; q++) {
      long r = q;
      double s = 0;
      for (int d = 0; d < D - 1; d++) {
        p.x [d] = r % len;
        r /= len;
        s += (p.x [d] + 0.5 - c) * (p.x [d] + 0.5 - c);
      }
      double lo2 = (R - 0.5) * (R - 0.5) - s, hi2 = (R + 0.5) * (R + 0.5) - s;
      if (hi2 < 0)
        continue;
      int lo = lo2 > 0 ? (int) ceil (sqrt (lo2) - 0.5) : 0;
      int hi = (int) floor (sqrt (hi2) - 0.5);
      for (int t = lo; t <= hi; t++) {
        p.x [D - 1] = (int) (c + 0.5 + t);
        add (p);
        p.x [D - 1] = (int) (c - 0.5 - t);
        add (p);
      }
    }
#endif
  }
  else if (!strcmp (pattern, "river")) {
    /* along x: a meander n wide (2D), or a tube of diameter n (3D) */
    for (int x = 0; x < len; x++) {
      double u = (double) x / len;
      double y = c + 0.25 * len * sin (6.283185 * 3 * u) +
        0.05 * len * sin (6.283185 * 17 * u);
      p.x [0] = x;
#if D == 2
      for (int w = 0; w < n; w++) {
        p.x [1] = (int) (y - 0.5 * n) + w;
        add (p);
      }
#elif D == 3
      double z = c + 0.1 * len * cos (6.283185 * 5 * u);
      for (int a = 0; a < n; a++)
        for (int b = 0; b < n; b++)
          if ((a - 0.5 * n + 0.5) * (a - 0.5 * n + 0.5) +
              (b - 0.5 * n + 0.5) * (b - 0.5 * n + 0.5) <= 0.25 * n * n) {
            p.x [1] = (int) (y - 0.5 * n) + a;
            p.x [2] = (int) (z - 0.5 * n) + b;
            add (p);
          }
#endif
    }
  }
  else {
    fprintf (stderr, "unknown pattern '%s'\n", pattern);
    exit (1);
  }
  /* Morton order, without duplicates */
  for (size_t q = 0; q < npts; q++)
    pts [q].key = morton (pts [q].p);
  qsort (pts, npts, sizeof (Keyed), by_key);
  size_t u = 0;
  for (size_t q = 0; q < npts; q++)
    if (!u || pts [q].key != pts [u - 1].key)
      pts [u++] = pts [q];
  npts = u;
}

/* ---- the benchmark -------------------------------------------------- */

/* repeats of a loop over n items so that a phase does about 4M lookups */
static int reps (size_t n)
{
  size_t r = (4000000 + n - 1) / n;
  return r < 1 ? 1 : r;
}

#define VALUE(q) ((void *) (uintptr_t) (((q) + 1) << 4))

int main (int argc, char ** argv)
{
  if (argc < 4) {
    fprintf (stderr, "usage: %s mem|voxel PATTERN LEVEL [N]\n", argv [0]);
    return 1;
  }
  use_voxel = !strcmp (argv [1], "voxel");
  const char * pattern = argv [2];
  int level = atoi (argv [3]);
  long n = argc > 4 ? atol (argv [4]) : 1000000;
  len = 1 << level;
  make_points (pattern, n);

  Pt * order = malloc (npts * sizeof (Pt)), * shuf = malloc (npts * sizeof (Pt));
  for (size_t q = 0; q < npts; q++)
    order [q] = shuf [q] = pts [q].p;
  free (pts);
  for (size_t q = npts - 1; q > 0; q--) {
    size_t r = xorshift () % (q + 1);
    Pt t = shuf [q]; shuf [q] = shuf [r]; shuf [r] = t;
  }
  size_t nuni = 4000000;
  Pt * uni = malloc (nuni * sizeof (Pt));
  for (size_t q = 0; q < nuni; q++)
    for (int d = 0; d < D; d++)
      uni [q].x [d] = xorshift () % len;

  perf_init ();
  printf ("impl=%s dim=%d pattern=%s level=%d n=%zu", argv [1], D, pattern,
          level, npts);
  long rss0 = status_kb ("VmRSS"), pte0 = status_kb ("VmPTE");

  /* build */
  Mark m = mark ();
  for (size_t q = 0; q < npts; q++)
    assign (order [q], VALUE (q));
  report ("build", m, npts);
  long rss = status_kb ("VmRSS") - rss0, pte = status_kb ("VmPTE") - pte0;
  printf (" rss_kb=%ld pte_kb=%ld bytes_per_pt=%.1f", rss, pte,
          1024. * (rss + pte) / npts);

  uintptr_t sum = 0, check = 0;
  int r;

  /* lookups of every point: Morton order, then random order */
  r = reps (npts);
  m = mark ();
  for (int t = 0; t < r; t++)
    for (size_t q = 0; q < npts; q++)
      sum += (uintptr_t) get (order [q]);
  report ("ordered", m, (double) r * npts);
  check += sum;

  m = mark ();
  sum = 0;
  for (int t = 0; t < r; t++)
    for (size_t q = 0; q < npts; q++)
      sum += (uintptr_t) get (shuf [q]);
  report ("shuffled", m, (double) r * npts);
  check += sum;

  /* face neighbours, Morton order */
  r = reps (2 * D * npts);
  m = mark ();
  sum = 0;
  for (int t = 0; t < r; t++)
    for (size_t q = 0; q < npts; q++)
      for (int d = 0; d < D; d++)
        for (int s = -1; s <= 1; s += 2) {
          Pt p = order [q];
          p.x [d] += s;
          sum += (uintptr_t) get (p);
        }
  report ("nbr", m, (double) r * 2 * D * npts);
  check += sum;
  uintptr_t nbrsum = sum;

  if (use_voxel) {
    m = mark ();
    sum = 0;
    for (int t = 0; t < r; t++) {
      VMemCursor cur = {0};
      for (size_t q = 0; q < npts; q++)
        for (int d = 0; d < D; d++)
          for (int s = -1; s <= 1; s += 2) {
            Pt p = order [q];
            p.x [d] += s;
            sum += (uintptr_t) vmem_get_cursor (&vi, &cur, P (p));
          }
    }
    report ("nbrcur", m, (double) r * 2 * D * npts);
    if (sum != nbrsum) {
      fprintf (stderr, "cursor lookups differ\n");
      return 1;
    }
  }

  /* uniform over the level */
  m = mark ();
  sum = 0;
  size_t hits = 0;
  for (size_t q = 0; q < nuni; q++) {
    void * v = get (uni [q]);
    sum += (uintptr_t) v;
    hits += v != NULL;
  }
  report ("uniform", m, nuni);
  printf (" uniform_hits=%.3f", (double) hits / nuni);
  check += sum;

  /* free everything */
  m = mark ();
  for (size_t q = 0; q < npts; q++)
    release (shuf [q]);
  report ("free", m, npts);
  printf (" rss_after_free_kb=%ld pte_after_free_kb=%ld",
          status_kb ("VmRSS") - rss0, status_kb ("VmPTE") - pte0);
  for (size_t q = 0; q < npts; q += 97)
    if (get (order [q])) {
      fprintf (stderr, "point %zu still there after free\n", q);
      return 1;
    }

  printf (" checksum=%lx\n", (unsigned long) check);
  return 0;
}
