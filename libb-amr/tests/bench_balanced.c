/*
.. memindex.[ch] against voxel-memindex.[ch] on a 2:1 balanced tree, built
.. as tree.c builds one: cells that touch a front are refined, a family at
.. a time, and before a cell is refined every index of its level around it
.. (the 3^dimension - 1 box) must hold a node, the coarser leaf covering it
.. being refined first (tree_refine () and covering () in tree.c).
..
..   bench_balanced IMPL PATTERN DEPTH
..     IMPL     mem | voxel
..     PATTERN  shell   a circle (2D) or sphere (3D) of radius 0.35
..              river   a meandering line (2D) or tube (3D), two finest cells wide
..              dense   every cell refined: a uniform tree
..     DEPTH    the finest level
..
.. Every level of the tree goes in its own array, nodes in pre-order (Morton
.. order, parents first). Each node's data is a record {flags, id}, the
.. records in the same order, as pool allocation of families would put
.. them. Phases (ns per node, or per leaf, and per lookup):
..   build     insert every node
..   dfs       depth-first traversal from the root, in Morton order: the
..             node read, then its children looked up (children_of ())
..   list      the nodes in pre-order, each looked up (index_node () over
..             a level list, as ghosts_update () does)
..   nbr       for every leaf, its 3^dimension - 1 neighbours at its level,
..             each falling back to the coarser node covering it
..             (block_ghosts () and covering ())
..   free      every node, children first
.. For voxel, dfs, list and nbr run again with one cursor per level
.. (dfscur, listcur, nbrcur). One line of key=value pairs out:
.. tests/bench_balanced.sh tabulates.
.. Build: gcc -O2 -Ddimension=D -D_GNU_SOURCE -I. tests/bench_balanced.c
..        memindex.c voxel-memindex.c pool.c -lm
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include <time.h>
#include <sys/resource.h>

#include "memindex.h"
#include "voxel-memindex.h"

#define D    dimension
#define NCH  (1 << D)
#define MAXL 22

#if D == 1
# define P(p)  (p).x [0]
# define BOX   3
#elif D == 2
# define P(p)  (p).x [0], (p).x [1]
# define BOX   9
#else
# define P(p)  (p).x [0], (p).x [1], (p).x [2]
# define BOX   27
#endif

#define CALL(f, ...)  f (__VA_ARGS__)

typedef struct { int x [D]; int l; } Node;
typedef struct { int flags; int id; } Rec;      /* flags & 1: a leaf */

/* ---- the measured index, one array per level ------------------------ */

static int use_voxel;
static struct _Memindex mi [MAXL + 1];
static struct _VMemindex vi [MAXL + 1];

static inline Rec * get (Node p)
{
  if (use_voxel)
    return vmem_get (&vi [p.l], P (p));
  return mem_allocated (&mi [p.l], P (p)) ?
    (Rec *) CALL (mem_data, &mi [p.l], P (p)) : NULL;
}

static inline Rec * get_cursor (Node p, VMemCursor * cur)
{
  return vmem_get_cursor (&vi [p.l], &cur [p.l], P (p));
}

static void assign (Node p, void * b)
{
  if (use_voxel)
    vmem_assign (&vi [p.l], P (p), 1 << p.l, b);
  else
    mem_assign (&mi [p.l], P (p), 1 << p.l, b);
}

static void release (Node p)
{
  if (use_voxel)
    vmem_free (&vi [p.l], P (p), 1 << p.l, NULL);
  else
    mem_free (&mi [p.l], P (p), 1 << p.l, NULL);
}

/* ---- nodes ------------------------------------------------------------ */

static Node child (Node p, int c)
{
  Node q;
  for (int d = 0; d < D; d++)
    q.x [d] = 2 * p.x [d] + ((c >> d) & 1);
  q.l = p.l + 1;
  return q;
}

static Node parent (Node p)
{
  for (int d = 0; d < D; d++)
    p.x [d] >>= 1;
  p.l--;
  return p;
}

/* the m-th index of the 3^D box around p, i fastest, as index_in_box () */
static Node in_box (Node p, int m)
{
  for (int d = 0; d < D; d++, m /= 3)
    p.x [d] += m % 3 - 1;
  return p;
}

static int inside (Node p)
{
  for (int d = 0; d < D; d++)
    if (p.x [d] < 0 || p.x [d] >= 1 << p.l)
      return 0;
  return 1;
}

/* ---- the front -------------------------------------------------------- */

static int pattern;     /* 0 shell, 1 river, 2 dense */
static int depth;

static double meander (double x, double f, double a, double b)
{
  return 0.5 + a * sin (6.283185 * f * x) + b * sin (6.283185 * 17 * x);
}

/* does the cell p touch the front? */
static int front (Node p)
{
  double h = 1. / (1 << p.l), c [D], half = 0.5 * h * sqrt (D);
  for (int d = 0; d < D; d++)
    c [d] = (p.x [d] + 0.5) * h;
  if (pattern == 2)
    return 1;
  if (pattern == 0) {
    double r = 0;
    for (int d = 0; d < D; d++)
      r += (c [d] - 0.5) * (c [d] - 0.5);
    return fabs (sqrt (r) - 0.35) <= half;
  }
#if D == 1
  return fabs (c [0] - 0.5) <= half + 0.01;
#else
  /* the curve's range over the cell's x interval, sampled, against the
  .. cell grown by the width: steep parts of the meander are kept */
  double ylo = 2, yhi = -2, zlo = 2, zhi = -2, w = 2. / (1 << depth);
  for (int s = 0; s <= 16; s++) {
    double x = (p.x [0] + s / 16.) * h, y = meander (x, 3, 0.25, 0.05);
    ylo = fmin (ylo, y);
    yhi = fmax (yhi, y);
    double z = 0.5 + 0.1 * cos (6.283185 * 5 * x);
    zlo = fmin (zlo, z);
    zhi = fmax (zhi, z);
  }
  int touch = yhi >= p.x [1] * h - w && ylo <= (p.x [1] + 1) * h + w;
# if D == 3
  touch &= zhi >= p.x [2] * h - w && zlo <= (p.x [2] + 1) * h + w;
# else
  (void) zlo;
  (void) zhi;
# endif
  (void) half;
  return touch;
#endif
}

/* ---- building the tree, with its own (voxel) arrays ------------------- */

static struct _VMemindex gi [MAXL + 1];
static Node * lv [MAXL + 1];
static size_t nlv [MAXL + 1], caplv [MAXL + 1];
static char mark;       /* any non-NULL value for the generator's arrays */

static int g_has (Node p)
{
  return vmem_get (&gi [p.l], P (p)) != NULL;
}

static void g_add (Node p)
{
  vmem_assign (&gi [p.l], P (p), 1 << p.l, &mark);
  if (nlv [p.l] == caplv [p.l]) {
    caplv [p.l] = caplv [p.l] ? 2 * caplv [p.l] : 1024;
    lv [p.l] = realloc (lv [p.l], caplv [p.l] * sizeof (Node));
    if (!lv [p.l])
      abort ();
  }
  lv [p.l] [nlv [p.l]++] = p;
}

static Node covering (Node q)
{
  while (!g_has (q))
    q = parent (q);
  return q;
}

/* tree_refine (): the box around p first, then p's family */
static void refine (Node p)
{
  if (g_has (child (p, 0)))
    return;
  for (int m = 0; m < BOX; m++) {
    Node q = in_box (p, m);
    if (inside (q))
      while (!g_has (q))
        refine (covering (q));
  }
  for (int c = 0; c < NCH; c++)
    g_add (child (p, c));
}

/* pre-order: Morton key at the finest level, then the coarser first */
typedef struct { uint64_t key; Node p; } Keyed;


static uint64_t preorder_key (Node p)
{
  uint64_t k = 0;
  for (int b = 0; b < depth; b++)
    for (int d = 0; d < D; d++)
      k |= (uint64_t) (((p.x [d] << (depth - p.l)) >> b) & 1) << (b * D + d);
  return k << 5 | p.l;
}

static int by_key (const void * a, const void * b)
{
  uint64_t x = ((const Keyed *) a)->key, y = ((const Keyed *) b)->key;
  return x < y ? -1 : x > y;
}

/* ---- measuring -------------------------------------------------------- */

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

/* repeats of a pass of n operations, for about 4M operations in all */
static int reps (double n)
{
  int r = (int) (4e6 / n + 0.5);
  return r < 1 ? 1 : r;
}

/* ---- the phases ------------------------------------------------------- */

static Node * order;    /* every node, pre-order */
static size_t N;

/* depth-first from the root, Morton order: nodes visited, sum of ids */
static size_t dfs (int cursor, uintptr_t * sum)
{
  static Node stack [MAXL * NCH + 1];
  VMemCursor cur [MAXL + 1];
  memset (cur, 0, sizeof cur);
  size_t visited = 0;
  int sp = 0;
  stack [sp++] = order [0];
  while (sp) {
    Node p = stack [--sp];
    Rec * r = cursor ? get_cursor (p, cur) : get (p);
    visited++;
    *sum += r->id;
    if (!(r->flags & 1))
      for (int c = NCH - 1; c >= 0; c--)
        stack [sp++] = child (p, c);
  }
  return visited;
}

static void list (int cursor, uintptr_t * sum)
{
  VMemCursor cur [MAXL + 1];
  memset (cur, 0, sizeof cur);
  for (size_t q = 0; q < N; q++) {
    Rec * r = cursor ? get_cursor (order [q], cur) : get (order [q]);
    *sum += r->id;
  }
}

/* the box around every leaf: lookups made, of which coarser */
static void nbr (int cursor, Node * leaves, size_t nleaves, uintptr_t * sum,
                 size_t * lookups, size_t * coarser)
{
  VMemCursor cur [MAXL + 1];
  memset (cur, 0, sizeof cur);
  for (size_t q = 0; q < nleaves; q++)
    for (int m = 0; m < BOX; m++) {
      if (m == BOX / 2)
        continue;
      Node b = in_box (leaves [q], m);
      if (!inside (b))
        continue;
      Rec * r;
      (*lookups)++;
      while (!(r = cursor ? get_cursor (b, cur) : get (b))) {
        b = parent (b);
        (*lookups)++;
        (*coarser)++;
      }
      *sum += r->id;
    }
}

int main (int argc, char ** argv)
{
  if (argc < 4) {
    fprintf (stderr, "usage: %s mem|voxel shell|river|dense DEPTH\n", argv [0]);
    return 1;
  }
  use_voxel = !strcmp (argv [1], "voxel");
  pattern = !strcmp (argv [2], "shell") ? 0 : !strcmp (argv [2], "river") ? 1 :
    !strcmp (argv [2], "dense") ? 2 : -1;
  depth = atoi (argv [3]);
  if (pattern < 0 || depth < 1 || depth > MAXL || D * depth > 58) {
    fprintf (stderr, "bad pattern or depth\n");
    return 1;
  }

  /* the tree, level by level: the leaves of level l on the front refined */
  Node root = {{0}, 0};
  g_add (root);
  for (int l = 0; l < depth; l++)
    for (size_t q = 0; q < nlv [l]; q++)
      if (front (lv [l] [q]))
        refine (lv [l] [q]);

  for (int l = 0; l <= depth; l++)
    N += nlv [l];
  Keyed * k = malloc (N * sizeof (Keyed));
  size_t n = 0;
  for (int l = 0; l <= depth; l++)
    for (size_t q = 0; q < nlv [l]; q++) {
      k [n].p = lv [l] [q];
      k [n++].key = preorder_key (lv [l] [q]);
    }
  qsort (k, N, sizeof (Keyed), by_key);
  order = malloc (N * sizeof (Node));
  Rec * rec = malloc (N * sizeof (Rec));
  size_t nleaves = 0;
  for (size_t q = 0; q < N; q++) {
    order [q] = k [q].p;
    rec [q].id = q + 1;
    rec [q].flags = order [q].l == depth || !g_has (child (order [q], 0));
    nleaves += rec [q].flags;
  }
  free (k);
  Node * leaves = malloc (nleaves * sizeof (Node));
  nleaves = 0;
  for (size_t q = 0; q < N; q++)
    if (rec [q].flags)
      leaves [nleaves++] = order [q];
  size_t perlevel [MAXL + 1];
  for (int l = 0; l <= depth; l++) {
    perlevel [l] = nlv [l];
    free (lv [l]);
    vmem_destroy (&gi [l], 1 << l);
  }

  printf ("impl=%s dim=%d pattern=%s depth=%d nodes=%zu leaves=%zu "
          "finest=%zu", argv [1], D, argv [2], depth, N, nleaves,
          perlevel [depth]);
  long rss0 = status_kb ("VmRSS"), pte0 = status_kb ("VmPTE");

  /* build */
  double t = now ();
  long f = minflt ();
  for (size_t q = 0; q < N; q++)
    assign (order [q], &rec [q]);
  printf (" build_ns=%.2f build_flt=%ld", 1e9 * (now () - t) / N,
          minflt () - f);
  long rss = status_kb ("VmRSS") - rss0, pte = status_kb ("VmPTE") - pte0;
  printf (" rss_kb=%ld pte_kb=%ld bytes_per_node=%.1f", rss, pte,
          1024. * (rss + pte) / N);

  uintptr_t check = 0;
  for (int cursor = 0; cursor <= use_voxel; cursor++) {
    const char * suffix = cursor ? "cur" : "";
    uintptr_t sum = 0;
    int r = reps (N);

    t = now ();
    f = minflt ();
    size_t visited = 0;
    for (int e = 0; e < r; e++)
      visited += dfs (cursor, &sum);
    if (visited != (size_t) r * N) {
      fprintf (stderr, "dfs visited %zu of %zu\n", visited, (size_t) r * N);
      return 1;
    }
    printf (" dfs%s_ns=%.2f dfs%s_flt=%ld", suffix,
            1e9 * (now () - t) / visited, suffix, minflt () - f);

    t = now ();
    for (int e = 0; e < r; e++)
      list (cursor, &sum);
    printf (" list%s_ns=%.2f", suffix, 1e9 * (now () - t) / ((double) r * N));

    size_t lookups = 0, coarser = 0;
    r = reps ((BOX - 1) * (double) nleaves);
    t = now ();
    f = minflt ();
    for (int e = 0; e < r; e++)
      nbr (cursor, leaves, nleaves, &sum, &lookups, &coarser);
    double el = now () - t;
    printf (" nbr%s_leaf_ns=%.2f nbr%s_lookup_ns=%.2f nbr%s_flt=%ld", suffix,
            1e9 * el / ((double) r * nleaves), suffix, 1e9 * el / lookups,
            suffix, minflt () - f);
    if (!cursor)
      printf (" lookups_per_leaf=%.2f coarser_share=%.3f",
              (double) lookups / ((double) r * nleaves),
              (double) coarser / lookups);
    if (cursor && sum != check) {
      fprintf (stderr, "cursor lookups differ\n");
      return 1;
    }
    check = sum;
  }

  /* free, children first */
  t = now ();
  for (size_t q = N; q-- > 0;)
    release (order [q]);
  printf (" free_ns=%.2f", 1e9 * (now () - t) / N);
  for (size_t q = 0; q < N; q += 101)
    if (get (order [q])) {
      fprintf (stderr, "node %zu still there after free\n", q);
      return 1;
    }
  printf (" checksum=%lx\n", (unsigned long) check);
  return 0;
}
