/* tree_restart () (tree-output.h): a tree refined over half its domain
 * (moved and resized), global scalars f, g, h and a local one, written
 * with tree_output (), then:
 *   refused, the tree and its values left as they were: no file, and
 *   copies of the file broken one way each (the magic, the dimension, m,
 *   cut short, a level past depth max, two leaves out of order, a codec,
 *   a leaf missing), and a true file of blocks of 32 cells for a tree of 16
 *   restarted onto another tree (uniform, another domain): f and g filled
 *   (2 returned), the local one and h (not in the file) not; the domain
 *   and the leaves the file's (as many, each where the file has it); f
 *   the function and g the cell's size at every cell, bit for bit; their
 *   ghosts marked stale
 *   a tree of one leaf, the root (serial)
 *   restarted from a copy in the other precision: within a float's
 *   rounding
 *   refused with no tree
 * With -D_MPI on processes (tests/mpi-fork) as without; the file is
 * BAMR_OUT (default test_restart.bamr), the broken copies beside it. */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>

#include "tree.h"
#include "block-ghosts.h"
#include "tree-output.h"
#ifdef _MPI
  #include <mpi.h>
  #include "tree-mpi.h"
#endif

scalar f, g, h;                /* global: named "f", "g", "h" by cwrap */

static int fails = 0, checks = 0;

#define CHECK(c, ...) do {                                        \
    checks++;                                                     \
    if (!(c)) {                                                   \
      fprintf (stderr, "rank %d, %s:%d: ", tree_rank (), __FILE__, \
               __LINE__);                                         \
      fprintf (stderr, __VA_ARGS__);                              \
      fputc ('\n', stderr);                                       \
      fails++;                                                    \
    }                                                             \
  } while (0)

static const double X0 = -1.25, Y0 = 2., Z0 = 0.5, H0 = 3.;

static double fun (double x, double y, double z)
{
  return sin (2 * x) + x * y - 0.5 * z + 1.;
}

static void refine_half (int to)
{
  for (int pass = 0; pass < 4 * to; pass++) {
    const TreeCaches * tc = tree_caches ();
    int n = 0;
    Index * todo = malloc ((size_t) tc->nleaves * sizeof (Index));
    for (int q = 0; q < tc->nleaves; q++) {
      Index p = tc->leaves [q].p;
      if (p.level < to && (p.i + 0.5) / (1 << p.level) < 0.5)
        todo [n++] = p;
    }
    for (int a = 0; a < n; a++)
      tree_refine (todo [a]);
    free (todo);
    if (!n)
      break;
  }
}

static void barrier (void)
{
#ifdef _MPI
  int a = 0, b = 0;
  MPI_Allreduce (&a, &b, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
#endif
}

static long total_leaves (void)
{
  long n = tree_caches ()->nowned;
#ifdef _MPI
  long t = 0;
  MPI_Allreduce (&n, &t, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
  n = t;
#endif
  return n;
}

/* a leaf's key at depth 21 */
static uint64_t key_of (Index p)
{
  uint64_t k = 0, q [3] = { (uint64_t) p.i, 0, 0 };
#if dimension > 1
  q [1] = (uint64_t) p.j;
#endif
#if dimension > 2
  q [2] = (uint64_t) p.k;
#endif
  for (int d = 0; d < 3; d++)
    q [d] <<= 21 - p.level;
  for (int b = 0; b < 21; b++)
    for (int d = 0; d < dimension; d++)
      k |= (q [d] >> b & 1) << (dimension * b + d);
  return k;              /* a tree's leaves have distinct first keys */
}

/* the cells where f is not the function or g not the cell's size (tol 0:
   bit for bit), and where h or a are not zero */
static long bad_cells (double tol, scalar a)
{
  long bad = 0;
  foreach (reduction (+:bad)) {
    double want = fun (x, y, z);
    if (tol == 0. ? (f [] != (Real) want || g [] != (Real) Delta)
        : (fabs (f [] - want) > tol * (1 + fabs (want)) ||
           fabs (g [] - Delta) > tol * Delta))
      bad++;
    if (a >= 0 && (h [] != 0 || a [] != 0))
      bad++;
  }
  return bad;
}

static unsigned char * read_all (const char * path, size_t * n)
{
  FILE * fp = fopen (path, "rb");
  if (!fp)
    return NULL;
  fseek (fp, 0, SEEK_END);
  long len = ftell (fp);
  fseek (fp, 0, SEEK_SET);
  unsigned char * p = malloc ((size_t) len + 1);
  if (p && fread (p, 1, (size_t) len, fp) != (size_t) len)
    free (p), p = NULL;
  fclose (fp);
  *n = (size_t) len;
  return p;
}

static void write_all (const char * path, const unsigned char * p, size_t n)
{
  FILE * fp = fopen (path, "wb");
  if (fp) {
    fwrite (p, 1, n, fp);
    fclose (fp);
  }
}

static uint32_t get32 (const unsigned char * p)
{
  return (uint32_t) p [0] | (uint32_t) p [1] << 8 | (uint32_t) p [2] << 16 |
    (uint32_t) p [3] << 24;
}

static void put32 (unsigned char * p, uint32_t v)
{
  for (int b = 0; b < 4; b++)
    p [b] = (unsigned char) (v >> (8 * b));
}

int main (void)
{
  const char * out = getenv ("BAMR_OUT");
  char path [512], broken [540], other [540];
  snprintf (path, sizeof (path), "%s", out ? out : "test_restart.bamr");
  snprintf (other, sizeof (other), "%s.other.bamr", path);

  int to = dimension == 1 ? 4 : dimension == 2 ? 3 : 2;
  tree_replicated (to + 1, 4);
  tree_shift_resize (X0, Y0, Z0, H0);
  refine_half (to);
  /* the whole tree's leaves, to find the restarted ones among */
  const TreeCaches * tc = tree_caches ();
  long N0 = tc->nleaves;
  uint64_t * keys = malloc ((size_t) N0 * sizeof (uint64_t));
  for (long q = 0; q < N0; q++)
    keys [q] = key_of (tc->leaves [q].p);
#ifdef _MPI
  CHECK (tree_mpi_partition (NULL) == 0, "partition");
#endif
  scalar loc;                  /* local: named "loc", not restored */
  foreach () {
    f [] = fun (x, y, z);
    g [] = Delta;
    h [] = 7;
    loc [] = 3;
  }
  scalar list [3] = { f, g, loc };
  CHECK (tree_output (path, list, 3) == 0, "tree_output (%s)", path);

  /* refused: the tree and its values as they were */
  long n_before = tree_get_nleaves ();
  CHECK (tree_restart ("/nonexistent/x.bamr") == -1, "no file");
  size_t len = 0;
  unsigned char * file = tree_rank () == 0 ? read_all (path, &len) : NULL;
  const int H = (int) ((80 + 3 * 4 + strlen ("f") + strlen ("g") +
                        strlen ("loc") + 7) / 8 * 8);
  const char * what [8] = { "the magic", "the dimension", "m", "cut short",
                            "a level past depth max", "two leaves swapped",
                            "a codec", "a leaf missing" };
  uint64_t mm = 1;
  for (int d = 0; d < dimension; d++)
    mm *= 16;
  for (int b = 0; b < 8; b++) {
    snprintf (broken, sizeof (broken), "%s.broken%d", path, b);
    if (tree_rank () == 0 && file) {
      unsigned char * c = malloc (len);
      memcpy (c, file, len);
      size_t n = len;
      switch (b) {
      case 0: c [0] = 'X'; break;
      case 1: put32 (c + 16, dimension % 3 + 1); break;
      case 2: put32 (c + 20, 8); break;
      case 3: n = len - 8; break;
      case 4: put32 (c + H, 99); break;
      case 5: {
        unsigned char t [16];
        memcpy (t, c + H, 16);
        memcpy (c + H, c + H + 16, 16);
        memcpy (c + H + 16, t, 16);
        break;
      }
      case 6: put32 (c + 72, 1); break;
      case 7: {
        /* the last leaf out, the file consistent otherwise: only the
           leaves' cover tells */
        uint64_t N = (uint64_t) N0, rs = get32 (file + 24), at = 0;
        memcpy (c, file, (size_t) H);
        put32 (c + 32, (uint32_t) (N - 1));
        at = (uint64_t) H;
        memcpy (c + at, file + H, 16 * (N - 1));
        at += 16 * (N - 1);
        for (int k = 0; k < 3; k++) {
          memcpy (c + at, file + H + 16 * N + (uint64_t) k * N * mm * rs,
                  (N - 1) * mm * rs);
          at += (N - 1) * mm * rs;
        }
        n = at;
        break;
      }
      }
      write_all (broken, c, n);
      free (c);
    }
    barrier ();
    CHECK (tree_restart (broken) == -1, "refused: %s", what [b]);
    CHECK (tree_get_nleaves () == n_before && bad_cells (0., -1) == 0,
           "%s: the tree left as it was", what [b]);
    if (tree_rank () == 0)
      remove (broken);
  }

  /* the other precision, for later */
  if (tree_rank () == 0 && file) {
    uint32_t rs = get32 (file + 24), rs2 = rs == 8 ? 4 : 8;
    uint64_t N = (uint64_t) N0;
    uint64_t head = (uint64_t) H + 16 * N, cells = 3 * N * mm;
    size_t n2 = head + cells * rs2;
    unsigned char * c = malloc (n2);
    memcpy (c, file, head);
    put32 (c + 24, rs2);
    for (uint64_t i = 0; i < cells; i++) {
      double v;
      if (rs == 8)
        memcpy (&v, file + head + 8 * i, 8);
      else {
        float w;
        memcpy (&w, file + head + 4 * i, 4);
        v = w;
      }
      if (rs2 == 8)
        memcpy (c + head + 8 * i, &v, 8);
      else {
        float w = (float) v;
        memcpy (c + head + 4 * i, &w, 4);
      }
    }
    write_all (other, c, n2);
    free (c);
  }
  free (file);
  barrier ();

  /* a file of blocks of 32 cells, refused by a tree of 16 */
  char path32 [540];
  snprintf (path32, sizeof (path32), "%s.m32.bamr", path);
  tree_destroy ();
  tree_replicated (to + 1, 5);
  foreach ()
    f [] = 1;
  CHECK (tree_output (path32, &f, 1) == 0, "the file of 32 cells written");
  tree_destroy ();
  tree_replicated (to + 1, 4);
  long n16 = tree_get_nleaves ();
  CHECK (tree_restart (path32) == -1 && tree_get_nleaves () == n16 &&
         tree_get_block_m () == 16, "refused: blocks of 32 cells");
  if (tree_rank () == 0)
    remove (path32);

  /* restarted onto another tree; f and g's ghosts updated first, so that
     their marks are the restart's */
  tree_destroy ();
  tree_replicated (to + 1, 4);
  tree_shift_resize (0., 0., 0., 1.);
  ghosts_update (f), ghosts_update (g);
  CHECK (!ghosts_is_dirty (f) && !ghosts_is_dirty (g), "marks cleared");
  int r = tree_restart (path);
  CHECK (r == 2, "%d scalars filled, not 2 (f, g)", r);
  CHECK (tree_get_origin (0) == X0 && tree_get_size () == H0
         && (dimension < 2 || tree_get_origin (1) == Y0)
         && (dimension < 3 || tree_get_origin (2) == Z0), "the domain");
  CHECK (total_leaves () == N0, "%ld leaves, %ld written", total_leaves (),
         N0);
  tc = tree_caches ();
  int lost = 0;
  for (int q = 0; q < tc->nowned; q++) {
    uint64_t k = key_of (tc->owned [q].p);
    long lo = 0, hi = N0;
    while (lo < hi) {
      long c = lo + (hi - lo) / 2;
      if (keys [c] < k)
        lo = c + 1;
      else
        hi = c;
    }
    lost += lo == N0 || keys [lo] != k;
  }
  CHECK (!lost, "%d leaves not the file's", lost);
  long bad = bad_cells (0., loc);
  CHECK (bad == 0, "%ld cells wrong (f, g bit for bit; h, loc zero)", bad);
  CHECK (ghosts_is_dirty (f) && ghosts_is_dirty (g), "f, g's ghosts not "
         "marked stale");
  ghosts_update (f);

  /* from the other precision */
  r = tree_restart (other);
  CHECK (r == 2, "the other precision: %d filled", r);
  bad = bad_cells (2e-7, -1);
  CHECK (bad == 0, "the other precision: %ld cells off", bad);
  if (tree_rank () == 0)
    remove (other);

#ifndef _MPI
  /* one leaf, the root: written, then restarted onto a tree of level 1 */
  tree_destroy ();
  tree_replicated (to + 1, 4);
  Index root;
  memset (&root, 0, sizeof (root));
  CHECK (tree_coarsen (root) == 1, "the root coarsened");
  foreach () {
    f [] = fun (x, y, z);
    g [] = Delta;
  }
  scalar two [2] = { f, g };
  CHECK (tree_output (path, two, 2) == 0, "the root written");
  tree_destroy ();
  tree_replicated (to + 1, 4);
  CHECK (tree_restart (path) == 2 && tree_get_nleaves () == 1 &&
         tree_caches ()->leaves [0].p.level == 0 && bad_cells (0., -1) == 0,
         "the root restarted");
#endif

  /* no tree */
  tree_destroy ();
  CHECK (tree_restart (path) == -1, "no tree");
  free (keys);
  if (tree_rank () == 0) {
    if (fails)
      printf ("%d of %d checks failed\n", fails, checks);
    else
      printf ("%d checks passed, %ld leaves\n", checks, N0);
  }
  return fails != 0;
}
