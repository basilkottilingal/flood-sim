/* profile.h: regions timed a thread and a process, reported as a table, a
.. flame graph or a timeline */

#ifndef _GNU_SOURCE
# define _GNU_SOURCE                  /* clock_gettime, strdup */
#endif

#include <math.h>
#include <pthread.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#ifdef _MPI
# include <mpi.h>
#endif

#include "profile.h"

enum { DEPTH = 64 };                  /* regions open at once a thread */

static double now (void)
{
  struct timespec t;
  clock_gettime (CLOCK_MONOTONIC, &t);
  return (double) t.tv_sec + 1e-9 * (double) t.tv_nsec;
}

/* zeroed memory, or the program stops */
static void * alloc (size_t n)
{
  void * p = calloc (1, n ? n : 1);
  if (!p)
  {
    fprintf (stderr, "profile: out of memory\n");
    abort ();
  }
  return p;
}

/* p grown to hold n items of size, from *cap (doubled) */
static void * grow (void * p, int * cap, int n, size_t size)
{
  if (n < *cap)
    return p;
  int c = *cap ? 2 * *cap : 64;
  while (c <= n)
    c *= 2;
  void * q = realloc (p, (size_t) c * size);
  if (!q)
  {
    fprintf (stderr, "profile: out of memory\n");
    abort ();
  }
  memset ((char *) q + (size_t) *cap * size, 0, (size_t) (c - *cap) * size);
  *cap = c;
  return q;
}

/* ---------------------------------------------------------------------- */
/* the regions: names, ids                                                */
/* ---------------------------------------------------------------------- */

static pthread_mutex_t lock = PTHREAD_MUTEX_INITIALIZER;
static char ** names;                 /* by region id */
static int nnames, cnames;

/* the id of a name, made if new (under the lock) */
static int region_of (const char * name)
{
  for (int i = 0; i < nnames; i++)
    if (!strcmp (names [i], name))
      return i;
  names = grow (names, &cnames, nnames, sizeof (char *));
  names [nnames] = strdup (name);
  if (!names [nnames])
    abort ();
  return nnames++;
}

/* ---------------------------------------------------------------------- */
/* a thread: its tree of call paths, its open regions, its events         */
/* ---------------------------------------------------------------------- */

/*
.. A node of the calling-context tree: a region on a call path. Node 0 is
.. the root (no region); a node's children are a list, short in practice,
.. found by a walk.
*/
typedef struct
{
  int region, parent, child, sibling;
  long calls;
  double incl, self, longest;
} Node;

typedef struct { int node; double t0, inside; } Frame;
typedef struct { int region; double t0, t1; } Event;

typedef struct Thread
{
  Node * node;
  int nnode, cnode;
  Frame frame [DEPTH];
  int depth, too_deep;                /* past DEPTH: not timed, counted */
  Event * event;
  int nevent, cevent;
  long dropped;                       /* events past the buffer */
  int tid;
  struct Thread * next;
} Thread;

static _Thread_local Thread * self;
static Thread * threads;              /* every thread's, under the lock */
static int nthreads;

/* settings, from the environment */
static int on = 1, want_table = 1, want_flame, want_trace;
static int trace_max = 65536;
static char prefix [256] = "profile";
static int started, finished;
static double t_start;                /* the first region entered */
static double overhead;               /* an enter and a leave, seconds */

static Thread * thread_new (int keep)
{
  Thread * t = alloc (sizeof (Thread));
  t->node = grow (NULL, &t->cnode, 0, sizeof (Node));
  t->node [0] = (Node) { .region = -1, .parent = -1, .child = -1,
                         .sibling = -1 };
  t->nnode = 1;
  if (keep)
  {
    pthread_mutex_lock (&lock);
    t->tid = nthreads++;
    t->next = threads;
    threads = t;
    pthread_mutex_unlock (&lock);
  }
  return t;
}

static void thread_free (Thread * t)
{
  free (t->node), free (t->event), free (t);
}

/*
.. the child of node p for region r, made if new: under the lock, as a
.. snapshot may be reading the nodes (and a realloc would move them); a
.. path is new only the first time, so the lock is rare
*/
static int child_of (Thread * t, int p, int r)
{
  for (int c = t->node [p].child; c >= 0; c = t->node [c].sibling)
    if (t->node [c].region == r)
      return c;
  pthread_mutex_lock (&lock);
  t->node = grow (t->node, &t->cnode, t->nnode, sizeof (Node));
  int c = t->nnode;
  t->node [c] = (Node) { .region = r, .parent = p, .child = -1,
                         .sibling = t->node [p].child };
  t->node [p].child = c;
  t->nnode++;
  pthread_mutex_unlock (&lock);
  return c;
}

/* enter region r on thread t: 0 if too deep (not timed) */
static int enter (Thread * t, int r)
{
  if (t->depth == DEPTH)
  {
    t->too_deep++;
    return 0;
  }
  int p = t->depth ? t->frame [t->depth - 1].node : 0;
  t->frame [t->depth++] = (Frame) { child_of (t, p, r), now (), 0 };
  return 1;
}

static void leave (Thread * t, int record)
{
  double t1 = now ();
  Frame * f = &t->frame [--t->depth];
  Node * n = &t->node [f->node];
  double dt = t1 - f->t0;
  n->calls++;
  n->incl += dt;
  n->self += dt - f->inside;
  n->longest = dt > n->longest ? dt : n->longest;
  if (t->depth)
    t->frame [t->depth - 1].inside += dt;
  if (record && want_trace)
  {
    if (t->nevent < trace_max)
    {
      if (t->nevent == t->cevent)     /* a snapshot may be reading */
      {
        pthread_mutex_lock (&lock);
        t->event = grow (t->event, &t->cevent, t->nevent, sizeof (Event));
        pthread_mutex_unlock (&lock);
      }
      t->event [t->nevent++] = (Event) { n->region, f->t0, t1 };
    }
    else
      t->dropped++;
  }
}

/* ---------------------------------------------------------------------- */
/* start: the settings, the cost of a region                              */
/* ---------------------------------------------------------------------- */

static void finish_at_exit (void);

/* the settings from PROFILE and PROFILE_OUT; the overhead, on a thread of
.. its own not kept */
static void start (void)
{
  const char * p = getenv ("PROFILE");
  if (p)
  {
    on = !strstr (p, "off");
    want_table = strstr (p, "table") != NULL;
    want_flame = strstr (p, "flame") != NULL;
    want_trace = strstr (p, "trace") != NULL;
  }
  const char * o = getenv ("PROFILE_OUT");
  if (o && o [0])
    snprintf (prefix, sizeof prefix, "%s", o);
  const char * m = getenv ("PROFILE_TRACE_EVENTS");
  if (m && atoi (m) > 0)
    trace_max = atoi (m);
  /* an enter and a leave, the fastest of a few rounds */
  Thread * t = thread_new (0);
  overhead = 1;
  for (int round = 0; round < 5; round++)
  {
    double t0 = now ();
    for (int k = 0; k < 1000; k++)
      if (enter (t, 0))
        leave (t, 0);
    double c = (now () - t0) / 1000;
    overhead = c < overhead ? c : overhead;
  }
  thread_free (t);
  t_start = now ();
  atexit (finish_at_exit);
}

static pthread_once_t once = PTHREAD_ONCE_INIT;

static void start_once (void)
{
  start ();
  started = 1;
}

int profile_enter (int * id, const char * name)
{
  int r = __atomic_load_n (id, __ATOMIC_ACQUIRE);
  if (r < 0)
  {
    /* not under the lock: start () times regions, which may take it */
    pthread_once (&once, start_once);
    pthread_mutex_lock (&lock);
    r = region_of (name);
    pthread_mutex_unlock (&lock);
    __atomic_store_n (id, r, __ATOMIC_RELEASE);
  }
  if (!on || finished)
    return 0;
  if (!self)
    self = thread_new (1);
  return enter (self, r);
}

void profile_leave (ProfileScope * scope)
{
  if (scope->on && self && self->depth)
    leave (self, 1);
}

/* ---------------------------------------------------------------------- */
/* a process's paths: every thread's, merged by path                      */
/* ---------------------------------------------------------------------- */

/*
.. A path is known by a 64-bit hash of its regions' names, root first
.. (FNV-1a, a NUL after each name), the same on every thread and
.. process; its parent by the hash of the path less its last region (0: a
.. root).
*/
typedef struct
{
  uint64_t key, parent;
  double calls, incl, self, longest;
  const char * name;                  /* the last region's (local) */
} Path;

static uint64_t hash_name (uint64_t h, const char * s)
{
  do
    h = (h ^ (unsigned char) *s) * 1099511628211u;
  while (*s++);
  return h;
}

static int path_cmp (const void * a, const void * b)
{
  uint64_t x = ((const Path *) a)->key, y = ((const Path *) b)->key;
  return (x > y) - (x < y);
}

/*
.. Every thread's paths, merged: calls summed, times the largest thread's.
.. Sorted by key; *n of them. The threads go on while this reads them: a
.. call in progress is not counted until it ends, and a number read as it
.. changes is a word written whole on the machines in use.
*/
static Path * local_paths (int * n)
{
  pthread_mutex_lock (&lock);
  int total = 0, k = 0;
  for (Thread * t = threads; t; t = t->next)
    total += t->nnode - 1;
  Path * p = alloc ((size_t) total * sizeof (Path));
  uint64_t * key = NULL;
  int ckey = 0;
  for (Thread * t = threads; t; t = t->next)
  {
    int nn = t->nnode;
    key = grow (key, &ckey, nn, sizeof (uint64_t));
    key [0] = 14695981039346656037u;   /* FNV's offset basis: the root */
    for (int i = 1; i < nn; i++)       /* parents come before children */
    {
      Node * x = &t->node [i];
      key [i] = hash_name (key [x->parent], names [x->region]);
      p [k++] = (Path) { key [i], x->parent ? key [x->parent] : 0,
                         (double) x->calls, x->incl, x->self, x->longest,
                         names [x->region] };
    }
  }
  pthread_mutex_unlock (&lock);
  free (key);
  qsort (p, (size_t) k, sizeof (Path), path_cmp);
  int m = 0;
  for (int i = 0; i < k; i++)
    if (m && p [m - 1].key == p [i].key)
    {
      Path * q = &p [m - 1];
      q->calls += p [i].calls;
      q->incl = fmax (q->incl, p [i].incl);
      q->self = fmax (q->self, p [i].self);
      q->longest = fmax (q->longest, p [i].longest);
    }
    else
      p [m++] = p [i];
  *n = m;
  return p;
}

/* ---------------------------------------------------------------------- */
/* a snapshot: the paths of every process, on rank 0                      */
/* ---------------------------------------------------------------------- */

/*
.. On rank 0: every path any process has, sorted by key, with its parent,
.. its last region's name, and the statistics over the processes.
*/
enum { SUM_CALLS, SUM_INCL, SUM_SELF, NSUM };
enum { MAX_INCL, MAX_SELF, MAX_LONGEST, NMAX };
typedef struct
{
  int n, nprocs, number;
  double wall;                        /* since the first region, rank 0 */
  uint64_t * key, * parent;
  char ** name;
  double * sum, * max, * min;         /* NSUM, NMAX, 1 a path */
} Report;

static void report_free (Report * r)
{
  for (int i = 0; i < r->n; i++)
    free (r->name [i]);
  free (r->key), free (r->parent), free (r->name);
  free (r->sum), free (r->max), free (r->min);
  memset (r, 0, sizeof * r);
}

/* the index of key in a sorted list, or -1 */
static int find_key (const uint64_t * k, int n, uint64_t key)
{
  int lo = 0, hi = n - 1;
  while (lo <= hi)
  {
    int mid = (lo + hi) / 2;
    if (k [mid] == key)
      return mid;
    if (k [mid] < key)
      lo = mid + 1;
    else
      hi = mid - 1;
  }
  return -1;
}

/* the statistics of local paths p [0..n) in the order of keys [0..m) */
static void stats_in_order (const Path * p, int n, const uint64_t * keys,
                            int m, double * sum, double * max, double * min)
{
  memset (sum, 0, (size_t) m * NSUM * sizeof (double));
  memset (max, 0, (size_t) m * NMAX * sizeof (double));
  memset (min, 0, (size_t) m * sizeof (double));
  for (int i = 0; i < n; i++)
  {
    int j = find_key (keys, m, p [i].key);
    if (j < 0)
      continue;
    sum [j * NSUM + SUM_CALLS] = p [i].calls;
    sum [j * NSUM + SUM_INCL] = p [i].incl;
    sum [j * NSUM + SUM_SELF] = p [i].self;
    max [j * NMAX + MAX_INCL] = p [i].incl;
    max [j * NMAX + MAX_SELF] = p [i].self;
    max [j * NMAX + MAX_LONGEST] = p [i].longest;
    min [j] = p [i].incl;
  }
}

/* the report of one process alone, from its paths */
static Report report_local (Path * p, int n, int number)
{
  Report r = { .n = n, .nprocs = 1, .number = number,
               .wall = now () - t_start };
  r.key = alloc ((size_t) n * sizeof (uint64_t));
  r.parent = alloc ((size_t) n * sizeof (uint64_t));
  r.name = alloc ((size_t) n * sizeof (char *));
  for (int i = 0; i < n; i++)
  {
    r.key [i] = p [i].key, r.parent [i] = p [i].parent;
    r.name [i] = strdup (p [i].name);
  }
  r.sum = alloc ((size_t) n * NSUM * sizeof (double));
  r.max = alloc ((size_t) n * NMAX * sizeof (double));
  r.min = alloc ((size_t) n * sizeof (double));
  stats_in_order (p, n, r.key, n, r.sum, r.max, r.min);
  return r;
}

/* ---------------------------------------------------------------------- */
/* the outputs, from a report                                             */
/* ---------------------------------------------------------------------- */

/* the children of path i (-1: the roots), largest mean inclusive first */
static int * children (const Report * r, int i, int * nc)
{
  uint64_t parent = i < 0 ? 0 : r->key [i];
  int * c = alloc ((size_t) r->n * sizeof (int)), k = 0;
  for (int j = 0; j < r->n; j++)
    if (r->parent [j] == parent)
      c [k++] = j;
  for (int a = 1; a < k; a++)         /* insertion sort: few children */
    for (int b = a; b > 0 && r->sum [c [b] * NSUM + SUM_INCL]
                              > r->sum [c [b - 1] * NSUM + SUM_INCL]; b--)
    {
      int x = c [b];
      c [b] = c [b - 1], c [b - 1] = x;
    }
  *nc = k;
  return c;
}

static void table_rows (FILE * f, const Report * r, int i, int depth)
{
  int nc;
  int * c = children (r, i, &nc);
  double P = r->nprocs;
  for (int k = 0; k < nc; k++)
  {
    int j = c [k];
    const double * s = &r->sum [j * NSUM], * m = &r->max [j * NMAX];
    double mean = s [SUM_INCL] / P;
    /* indented by depth; past 8 levels, the depth written instead */
    char label [96];
    if (depth <= 8)
      snprintf (label, sizeof label, "%*s%s", 2 * depth, "", r->name [j]);
    else
      snprintf (label, sizeof label, "%16s[%d] %s", "", depth, r->name [j]);
    double calls = s [SUM_CALLS] / P;     /* a mean: may not be whole */
    fprintf (f, "%-36.36s %10.*f %10.4f %10.4f %10.4f %10.4f %6.2f %6.1f\n",
             label, calls == floor (calls) ? 0 : 2, calls, s [SUM_SELF] / P,
             m [MAX_SELF], mean,
             m [MAX_INCL], mean > 0 ? m [MAX_INCL] / mean : 1,
             r->wall > 0 ? 100 * mean / r->wall : 0);
    table_rows (f, r, j, depth + 1);
  }
  free (c);
}

static void table (FILE * f, const Report * r)
{
  char which [32] = "";
  if (r->number >= 0)
    snprintf (which, sizeof which, ", snapshot %d", r->number);
  fprintf (f, "=== profile%s: %d process%s, wall %.4f s (rank 0), "
           "%.0f ns an enter and leave ===\n", which, r->nprocs,
           r->nprocs > 1 ? "es" : "", r->wall, overhead * 1e9);
  fprintf (f, "%-36s %10s %10s %10s %10s %10s %6s %6s\n", "region",
           "calls", "self mean", "self max", "incl mean", "incl max",
           "imbal", "%wall");
  table_rows (f, r, -1, 0);
}

/* the path of i, names joined by ';', into buf */
static void path_string (const Report * r, int i, char * buf, size_t len)
{
  int chain [DEPTH + 1], n = 0;
  for (int j = i; j >= 0 && n <= DEPTH; )
  {
    chain [n++] = j;
    j = r->parent [j] ? find_key (r->key, r->n, r->parent [j]) : -1;
  }
  buf [0] = '\0';
  for (int k = n - 1; k >= 0; k--)
  {
    size_t used = strlen (buf);
    snprintf (buf + used, len - used, "%s%s", k < n - 1 ? ";" : "",
              r->name [chain [k]]);
  }
}

static void flame (FILE * f, const Report * r)
{
  char path [4096];
  for (int i = 0; i < r->n; i++)
  {
    long us = lround (1e6 * r->sum [i * NSUM + SUM_SELF] / r->nprocs);
    if (us <= 0)
      continue;
    path_string (r, i, path, sizeof path);
    fprintf (f, "%s %ld\n", path, us);
  }
}

static FILE * open_output (const char * suffix, const char * mode)
{
  char fn [512];
  snprintf (fn, sizeof fn, "%s%s", prefix, suffix);
  FILE * f = fopen (fn, mode);
  if (!f)
    fprintf (stderr, "profile: cannot write %s\n", fn);
  return f;
}

/* the final outputs of rank 0 (or of a process alone: suffix) */
static void write_outputs (const Report * r, const char * suffix)
{
  char s [64];
  if (want_table)
  {
    if (!suffix [0])
      table (stdout, r), fflush (stdout);
    snprintf (s, sizeof s, "%s.txt", suffix);
    FILE * f = open_output (s, "w");
    if (f)
      table (f, r), fclose (f);
  }
  if (want_flame)
  {
    snprintf (s, sizeof s, "%s.folded", suffix);
    FILE * f = open_output (s, "w");
    if (f)
      flame (f, r), fclose (f);
  }
}

/* ---------------------------------------------------------------------- */
/* the trace: every process's events, rank 0 writing them                 */
/* ---------------------------------------------------------------------- */

/* JSON's string: quotes and backslashes escaped, control bytes dropped */
static void json_string (FILE * f, const char * s)
{
  fputc ('"', f);
  for (; *s; s++)
    if (*s == '"' || *s == '\\')
      fprintf (f, "\\%c", *s);
    else if ((unsigned char) *s >= 0x20)
      fputc (*s, f);
  fputc ('"', f);
}

/*
.. A process's events as bytes: the names (nnames, then each with its
.. NUL), then a record an event, its thread first. Times from t0, in
.. microseconds.
*/
typedef struct { int region, tid; double ts, dur; } TraceEvent;

static char * trace_pack (double t0, int * bytes, long * dropped)
{
  pthread_mutex_lock (&lock);
  size_t nb = sizeof (int) * 2;
  long ne = 0;
  *dropped = 0;
  for (int i = 0; i < nnames; i++)
    nb += strlen (names [i]) + 1;
  for (Thread * t = threads; t; t = t->next)
    ne += t->nevent, *dropped += t->dropped;
  nb += (size_t) ne * sizeof (TraceEvent);
  char * b = alloc (nb), * at = b;
  int h [2] = { nnames, (int) ne };
  memcpy (at, h, sizeof h), at += sizeof h;
  for (int i = 0; i < nnames; i++)
  {
    size_t l = strlen (names [i]) + 1;
    memcpy (at, names [i], l), at += l;
  }
  for (Thread * t = threads; t; t = t->next)
    for (int i = 0; i < t->nevent; i++)
    {
      TraceEvent e = { t->event [i].region, t->tid,
                       1e6 * (t->event [i].t0 - t0),
                       1e6 * (t->event [i].t1 - t->event [i].t0) };
      memcpy (at, &e, sizeof e), at += sizeof e;
    }
  pthread_mutex_unlock (&lock);
  *bytes = (int) nb;
  return b;
}

/* the events of process pid, from its bytes, into the JSON file */
static void trace_write (FILE * f, int pid, const char * b, int * first)
{
  int h [2];
  memcpy (h, b, sizeof h);
  const char ** nm = alloc ((size_t) (h [0] ? h [0] : 1) * sizeof (char *));
  const char * at = b + sizeof h;
  for (int i = 0; i < h [0]; i++)
    nm [i] = at, at += strlen (at) + 1;
  fprintf (f, "%s{\"name\":\"process_name\",\"ph\":\"M\",\"pid\":%d,"
           "\"args\":{\"name\":\"rank %d\"}}", *first ? "" : ",\n", pid, pid);
  *first = 0;
  for (int i = 0; i < h [1]; i++)
  {
    TraceEvent e;
    memcpy (&e, at + (size_t) i * sizeof e, sizeof e);
    fprintf (f, ",\n{\"name\":");
    json_string (f, nm [e.region]);
    fprintf (f, ",\"ph\":\"X\",\"pid\":%d,\"tid\":%d,\"ts\":%.3f,"
             "\"dur\":%.3f}", pid, e.tid, e.ts, e.dur);
  }
  free (nm);
}

/* ---------------------------------------------------------------------- */
/* MPI: snapshots by non-blocking collectives                             */
/* ---------------------------------------------------------------------- */

static Report last;                   /* rank 0's latest, or a process's */
static int snapshots;                 /* begun */

#ifdef _MPI
/*
.. A snapshot is a chain of non-blocking collectives on a communicator of
.. the profiler's own (MPI_COMM_WORLD duplicated: never matched with the
.. program's messages). Every process issues the same collectives in the
.. same order, each when it gets to a call of profile_snapshot () or
.. profile_finish (): the processes need not be there at once.
..   0. the first time, the communicator: MPI_Comm_idup (MPI_Comm_dup
..      would wait for every process);
..   1. the paths are copied (the snapshot's moment); MPI_Iallgather of
..      each process's count of paths new to the agreed list, and of the
..      bytes of their names;
..   2. MPI_Iallgatherv of the new paths' keys and parents, to every
..      process; MPI_Igatherv of their names, to rank 0;
..   3. the agreed list made larger by the new keys (sorted: the same on
..      every process); the statistics in its order; three MPI_Ireduce,
..      to rank 0: sums (calls, times), maxima, minima;
..   4. on rank 0, the report.
.. After the first snapshots, new paths are rare: step 2 is small.
*/
typedef struct
{
  int phase;                          /* 0 idle, 1-3 waiting for that
                                         step, 4 for the communicator */
  MPI_Request req [3];
  int nreq;
  Path * path;                        /* the copy: this process's */
  int npath;
  int mine [2];                       /* new paths, their name bytes */
  int * counts;                       /* 2 a process */
  uint64_t * mynew, * allnew;         /* key, parent pairs */
  int * kc, * kd;                     /* counts and displacements */
  char * mynames, * allnames;
  int * nc, * nd;
  double * sum, * max, * min, * rsum, * rmax, * rmin;
  int number;
} Snapshot;

static MPI_Comm comm = MPI_COMM_NULL;
static int comm_ready, rank, nprocs;
static Snapshot snap;
static uint64_t * agreed, * agreed_parent;   /* every process: the list */
static int nagreed, cagreed, cparent, cname;
static char ** agreed_name;           /* rank 0: a name a key */

/* MPI running? (the profiler's communicator comes with the first
.. snapshot) */
static int mpi_on (void)
{
  int init, done;
  MPI_Initialized (&init);
  MPI_Finalized (&done);
  if (!init || done)
    return 0;
  MPI_Comm_rank (MPI_COMM_WORLD, &rank);
  MPI_Comm_size (MPI_COMM_WORLD, &nprocs);
  return 1;
}

/* step 1, its collective: the counts of new paths */
static void snap_counts (void)
{
  Snapshot * s = &snap;
  MPI_Iallgather (s->mine, 2, MPI_INT, s->counts, 2, MPI_INT, comm,
                  &s->req [0]);
  s->nreq = 1, s->phase = 1;
}

/* step 1: the copy, and the new paths; the communicator first, the first
.. time (step 0) */
static void snap_begin (void)
{
  Snapshot * s = &snap;
  memset (s, 0, sizeof * s);
  s->number = snapshots++;
  s->path = local_paths (&s->npath);
  s->mynew = alloc ((size_t) (2 * (s->npath ? s->npath : 1))
                    * sizeof (uint64_t));
  size_t nb = 0;
  for (int i = 0; i < s->npath; i++)
    if (find_key (agreed, nagreed, s->path [i].key) < 0)
      nb += strlen (s->path [i].name) + 1;
  s->mynames = alloc (nb);
  for (int i = 0, at = 0; i < s->npath; i++)
    if (find_key (agreed, nagreed, s->path [i].key) < 0)
    {
      s->mynew [2 * s->mine [0]] = s->path [i].key;
      s->mynew [2 * s->mine [0] + 1] = s->path [i].parent;
      s->mine [0]++;
      size_t l = strlen (s->path [i].name) + 1;
      memcpy (s->mynames + at, s->path [i].name, l);
      at += (int) l;
    }
  s->mine [1] = (int) nb;
  s->counts = alloc ((size_t) (2 * nprocs) * sizeof (int));
  if (comm_ready)
  {
    snap_counts ();
    return;
  }
  MPI_Comm_idup (MPI_COMM_WORLD, &comm, &s->req [0]);
  s->nreq = 1, s->phase = 4;
}

/* step 2: the new keys to all, their names to rank 0 */
static void snap_keys (void)
{
  Snapshot * s = &snap;
  s->kc = alloc ((size_t) nprocs * sizeof (int));
  s->kd = alloc ((size_t) nprocs * sizeof (int));
  s->nc = alloc ((size_t) nprocs * sizeof (int));
  s->nd = alloc ((size_t) nprocs * sizeof (int));
  int nk = 0, nn = 0;
  for (int r = 0; r < nprocs; r++)
  {
    s->kc [r] = 2 * s->counts [2 * r], s->kd [r] = nk, nk += s->kc [r];
    s->nc [r] = s->counts [2 * r + 1], s->nd [r] = nn, nn += s->nc [r];
  }
  s->allnew = alloc ((size_t) nk * sizeof (uint64_t));
  s->allnames = alloc ((size_t) (rank == 0 ? nn : 1));
  MPI_Iallgatherv (s->mynew, 2 * s->mine [0], MPI_UINT64_T, s->allnew, s->kc,
                   s->kd, MPI_UINT64_T, comm, &s->req [0]);
  MPI_Igatherv (s->mynames, s->mine [1], MPI_CHAR, s->allnames, s->nc, s->nd,
                MPI_CHAR, 0, comm, &s->req [1]);
  s->nreq = 2, s->phase = 2;
}

/* step 3: the agreed list larger, the statistics in its order, reduced */
static void snap_reduce (void)
{
  Snapshot * s = &snap;
  const char * nm = s->allnames;
  for (int r = 0; r < nprocs; r++)
    for (int i = 0; i < s->counts [2 * r]; i++)
    {
      uint64_t key = s->allnew [s->kd [r] + 2 * i];
      uint64_t parent = s->allnew [s->kd [r] + 2 * i + 1];
      const char * name = nm;
      if (rank == 0)
        nm += strlen (nm) + 1;
      if (find_key (agreed, nagreed, key) >= 0)
        continue;                     /* new to two processes at once */
      /* insert in order: the same list on every process */
      int at = nagreed;
      while (at > 0 && agreed [at - 1] > key)
        at--;
      agreed = grow (agreed, &cagreed, nagreed, sizeof (uint64_t));
      agreed_parent = grow (agreed_parent, &cparent, nagreed,
                            sizeof (uint64_t));
      memmove (agreed + at + 1, agreed + at,
               (size_t) (nagreed - at) * sizeof (uint64_t));
      memmove (agreed_parent + at + 1, agreed_parent + at,
               (size_t) (nagreed - at) * sizeof (uint64_t));
      agreed [at] = key, agreed_parent [at] = parent;
      if (rank == 0)
      {
        agreed_name = grow (agreed_name, &cname, nagreed, sizeof (char *));
        memmove (agreed_name + at + 1, agreed_name + at,
                 (size_t) (nagreed - at) * sizeof (char *));
        agreed_name [at] = strdup (name);
      }
      nagreed++;
    }
  int m = nagreed;
  s->sum = alloc ((size_t) m * NSUM * sizeof (double));
  s->max = alloc ((size_t) m * NMAX * sizeof (double));
  s->min = alloc ((size_t) m * sizeof (double));
  stats_in_order (s->path, s->npath, agreed, m, s->sum, s->max, s->min);
  if (rank == 0)
  {
    s->rsum = alloc ((size_t) m * NSUM * sizeof (double));
    s->rmax = alloc ((size_t) m * NMAX * sizeof (double));
    s->rmin = alloc ((size_t) m * sizeof (double));
  }
  MPI_Ireduce (s->sum, s->rsum, m * NSUM, MPI_DOUBLE, MPI_SUM, 0, comm,
               &s->req [0]);
  MPI_Ireduce (s->max, s->rmax, m * NMAX, MPI_DOUBLE, MPI_MAX, 0, comm,
               &s->req [1]);
  MPI_Ireduce (s->min, s->rmin, m, MPI_DOUBLE, MPI_MIN, 0, comm,
               &s->req [2]);
  s->nreq = 3, s->phase = 3;
}

/* step 4, rank 0: the report; every process frees the snapshot */
static void snap_end (void)
{
  Snapshot * s = &snap;
  if (rank == 0)
  {
    report_free (&last);
    int m = nagreed;
    last = (Report) { .n = m, .nprocs = nprocs, .number = s->number,
                      .wall = now () - t_start, .sum = s->rsum,
                      .max = s->rmax, .min = s->rmin };
    last.key = alloc ((size_t) m * sizeof (uint64_t));
    last.parent = alloc ((size_t) m * sizeof (uint64_t));
    last.name = alloc ((size_t) m * sizeof (char *));
    memcpy (last.key, agreed, (size_t) m * sizeof (uint64_t));
    memcpy (last.parent, agreed_parent, (size_t) m * sizeof (uint64_t));
    for (int i = 0; i < m; i++)
      last.name [i] = strdup (agreed_name [i]);
    if (want_table && !finished)
    {
      FILE * f = open_output ("_snapshots.txt", s->number ? "a" : "w");
      if (f)
        table (f, &last), fputc ('\n', f), fclose (f);
    }
  }
  free (s->path), free (s->counts), free (s->mynew), free (s->allnew);
  free (s->kc), free (s->kd), free (s->nc), free (s->nd);
  free (s->mynames), free (s->allnames);
  free (s->sum), free (s->max), free (s->min);
  s->phase = 0;
}

/* the snapshot moved on as far as its collectives are done (wait: until
.. it is over) */
static void snap_progress (int wait)
{
  while (snap.phase)
  {
    int done = 1;
    if (wait)
      MPI_Waitall (snap.nreq, snap.req, MPI_STATUSES_IGNORE);
    else
      MPI_Testall (snap.nreq, snap.req, &done, MPI_STATUSES_IGNORE);
    if (!done)
      return;
    if (snap.phase == 4)
      comm_ready = 1, snap_counts ();
    else if (snap.phase == 1)
      snap_keys ();
    else if (snap.phase == 2)
      snap_reduce ();
    else
      snap_end ();
  }
}
#endif

void profile_snapshot (void)
{
  if (!started || finished)
    return;
#ifdef _MPI
  if (mpi_on ())
  {
    /* the last one finished first, waited for if need be: every process
     * then begins as many snapshots, so their collectives match */
    snap_progress (1);
    snap_begin ();
    snap_progress (0);
    return;
  }
#endif
  /* one process: the snapshot at once */
  int n;
  Path * p = local_paths (&n);
  report_free (&last);
  last = report_local (p, n, snapshots++);
  free (p);
  if (want_table)
  {
    FILE * f = open_output ("_snapshots.txt", last.number ? "a" : "w");
    if (f)
      table (f, &last), fputc ('\n', f), fclose (f);
  }
}

void profile_finish (void)
{
  if (!started || finished)
    return;
#ifdef _MPI
  int init, done;
  MPI_Initialized (&init);
  MPI_Finalized (&done);
  if (init && done)
  {
    /* too late for collectives: each process its own table */
    int n, r;
    MPI_Comm_rank (MPI_COMM_SELF, &r);
    Path * p = local_paths (&n);
    Report one = report_local (p, n, -1);
    free (p);
    char suffix [32];
    snprintf (suffix, sizeof suffix, "_rank%d", r);
    finished = 1;
    write_outputs (&one, suffix);
    report_free (&one);
    return;
  }
  if (mpi_on ())
  {
    snap_progress (1);                /* the last one, finished */
    finished = 1;                     /* no more timing from here */
    snap_begin ();
    snap_progress (1);
    last.number = -1;
    if (rank == 0)
      write_outputs (&last, "");
    if (want_trace)
    {
      /* the clocks aligned at a barrier: everyone leaves it together */
      MPI_Barrier (comm);
      double t0 = now ();
      long dropped, all;
      int bytes;
      char * b = trace_pack (t0, &bytes, &dropped);
      int * size = rank == 0 ? alloc ((size_t) nprocs * sizeof (int)) : NULL;
      int * disp = rank == 0 ? alloc ((size_t) nprocs * sizeof (int)) : NULL;
      MPI_Gather (&bytes, 1, MPI_INT, size, 1, MPI_INT, 0, comm);
      MPI_Reduce (&dropped, &all, 1, MPI_LONG, MPI_SUM, 0, comm);
      char * every = NULL;
      if (rank == 0)
      {
        long total = 0;
        for (int r = 0; r < nprocs; r++)
          disp [r] = (int) total, total += size [r];
        every = alloc ((size_t) total);
      }
      MPI_Gatherv (b, bytes, MPI_BYTE, every, size, disp, MPI_BYTE, 0, comm);
      if (rank == 0)
      {
        FILE * f = open_output (".trace.json", "w");
        if (f)
        {
          int first = 1;
          fprintf (f, "{\"displayTimeUnit\":\"ms\",\"traceEvents\":[\n");
          for (int r = 0; r < nprocs; r++)
            trace_write (f, r, every + disp [r], &first);
          fprintf (f, "\n]}\n");
          fclose (f);
        }
        if (all)
          fprintf (stderr, "profile: %ld calls not in the trace (more than "
                   "PROFILE_TRACE_EVENTS = %d a thread)\n", all, trace_max);
      }
      free (b), free (size), free (disp), free (every);
    }
    MPI_Comm_free (&comm);
    comm_ready = 0;
    return;
  }
#endif
  finished = 1;
  int n;
  Path * p = local_paths (&n);
  report_free (&last);
  last = report_local (p, n, -1);
  free (p);
  write_outputs (&last, "");
  if (want_trace)
  {
    long dropped;
    int bytes, first = 1;
    char * b = trace_pack (t_start, &bytes, &dropped);
    FILE * f = open_output (".trace.json", "w");
    if (f)
    {
      fprintf (f, "{\"displayTimeUnit\":\"ms\",\"traceEvents\":[\n");
      trace_write (f, 0, b, &first);
      fprintf (f, "\n]}\n");
      fclose (f);
    }
    if (dropped)
      fprintf (stderr, "profile: %ld calls not in the trace (more than "
               "PROFILE_TRACE_EVENTS = %d a thread)\n", dropped, trace_max);
    free (b);
  }
}

static void finish_at_exit (void)
{
  profile_finish ();
}

int profile_query (const char * path, double * calls, double * incl,
                   double * self_time, double * imbalance)
{
  char p [4096];
  for (int i = 0; i < last.n; i++)
  {
    path_string (&last, i, p, sizeof p);
    if (strcmp (p, path))
      continue;
    double P = last.nprocs;
    if (calls)
      *calls = last.sum [i * NSUM + SUM_CALLS] / P;
    if (incl)
      *incl = last.sum [i * NSUM + SUM_INCL] / P;
    if (self_time)
      *self_time = last.sum [i * NSUM + SUM_SELF] / P;
    if (imbalance)
    {
      double mean = last.sum [i * NSUM + SUM_INCL] / P;
      *imbalance = mean > 0 ? last.max [i * NMAX + MAX_INCL] / mean : 1;
    }
    return 1;
  }
  return 0;
}
