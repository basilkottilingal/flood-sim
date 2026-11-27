/* tree_output () (tree-output.h): a tree refined over half its domain
 * (moved and resized), two scalars, written, then read back with a reader
 * of the test's own:
 *   the header field by field; the leaves along the curve (keys strictly
 *   increasing), as many as the tree has, covering the domain (their
 *   volumes sum to h0^D)
 *   the first scalar the function at each cell's centre, the second the
 *   cell's size, both from the leaf's index alone
 *   each process's own part of the file, bit for bit, its blocks' Reals
 *   ".bamr" added to a name with no extension, not to one with
 *   the names the scalars' own, "f" and "g"
 *   with -D_MPI, a tree not partitioned (tree_replicated ()): its leaves
 *   once, not once a process
 *   refused: n < 0, no list, a scalar not in use, a scalar with no name,
 *   one name twice, a file that cannot be
 *   opened; an older, longer file at the name cut to the new one
 * With -D_MPI on processes (tests/mpi-fork) as without: tests/run.sh then
 * compares the files of 1 and 3 processes and of the serial build, byte
 * for byte. The file is BAMR_OUT (default test_output.bamr). */

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

/* the leaves coarser than `to` in the lower half along x refined */
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

static uint32_t get32 (const unsigned char * p)
{
  return (uint32_t) p [0] | (uint32_t) p [1] << 8 | (uint32_t) p [2] << 16 |
    (uint32_t) p [3] << 24;
}

static uint64_t get64 (const unsigned char * p)
{
  return (uint64_t) get32 (p) | (uint64_t) get32 (p + 4) << 32;
}

static double get_double (const unsigned char * p)
{
  uint64_t v = get64 (p);
  double x;
  memcpy (&x, &v, sizeof (x));
  return x;
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

/* a leaf's place on the curve at depth 21: its first finest cell's key */
static uint64_t key (const int32_t * v)
{
  uint64_t k = 0, q [3];
  for (int d = 0; d < 3; d++)
    q [d] = (uint64_t) v [1 + d] << (21 - v [0]);
  for (int b = 0; b < 21; b++)
    for (int d = 0; d < dimension; d++)
      k |= (q [d] >> b & 1) << (dimension * b + d);
  return k;
}

int main (void)
{
  const char * base = getenv ("BAMR_OUT");
  char path [512], noext [512], other [520];
  snprintf (path, sizeof (path), "%s", base ? base : "test_output.bamr");
  /* path without its extension, and another name with one */
  snprintf (noext, sizeof (noext), "%s", path);
  char * dot = strrchr (noext, '.');
  if (dot)
    *dot = '\0';
  snprintf (other, sizeof (other), "%s.out", noext);

  int to = dimension == 1 ? 4 : dimension == 2 ? 3 : 2;
  tree_replicated (to + 1, 4);
  tree_shift_resize (X0, Y0, Z0, H0);
  refine_half (to);
  long nleaves = tree_caches ()->nleaves;
#ifdef _MPI
  CHECK (tree_mpi_partition (NULL) == 0, "partition");
#endif
  scalar f, g;                 /* named "f" and "g" by cwrap */
  foreach () {
    f [] = fun (x, y, z);
    g [] = Delta;
  }
  scalar list [2] = { f, g };

  /* refused, on every process alike */
  scalar gone = scalar_new ();
  scalar_free (gone);
  scalar anon = scalar_new ();          /* no name */
  scalar twice [2] = { f, f };
  CHECK (tree_output (noext, list, -1) == -1, "n < 0");
  CHECK (tree_output (noext, &anon, 1) == -1, "a scalar with no name");
  CHECK (tree_output (noext, twice, 2) == -1, "a name twice");
  scalar_free (anon);
  CHECK (tree_output (noext, NULL, 2) == -1, "no list");
  CHECK (tree_output (noext, &gone, 1) == -1, "a scalar not in use");
  CHECK (tree_output ("/nonexistent/dir/x", list, 2) == -1,
         "a file that cannot be opened");

  /* written: ".bamr" added to noext, nothing to other; an older, longer
     file there first, which must not be left at the end */
  char added [520];
  snprintf (added, sizeof (added), "%s.bamr", noext);
  if (tree_rank () == 0) {
    FILE * old = fopen (added, "wb");
    for (int i = 0; old && i < (1 << 20); i++)
      fputc ('x', old);
    if (old)
      fclose (old);
  }
  CHECK (tree_output (noext, list, 2) == 0, "tree_output (%s)", noext);
  CHECK (tree_output (other, list, 2) == 0, "tree_output (%s)", other);

  size_t len = 0, len2 = 0;
  uint64_t N = 0;
  unsigned char * file = read_all (added, &len);
  unsigned char * file2 = read_all (other, &len2);
  CHECK (file && file2 && len == len2 && !memcmp (file, file2, len),
         "%s and %s the same", added, other);
  free (file2);
  if (!file)
    goto done;

  /* the header */
  const int m = tree_get_block_m ();
  uint64_t mm = 1;
  for (int d = 0; d < dimension; d++)
    mm *= (uint64_t) m;
  CHECK (!memcmp (file, "BAMRGRID", 8) && get32 (file + 8) == 1 &&
         get32 (file + 12) == 0x01020304u, "magic, version, byte order");
  CHECK (get32 (file + 16) == dimension && get32 (file + 20) == (uint32_t) m
         && get32 (file + 24) == sizeof (Real) && get32 (file + 28) == 2,
         "dimension, m, Real, scalars");
  N = get64 (file + 32);
  CHECK (N == (uint64_t) nleaves, "%llu leaves, the tree has %ld",
         (unsigned long long) N, nleaves);
  const double o [3] = { X0, dimension > 1 ? Y0 : 0., dimension > 2 ? Z0
                         : 0. };
  for (int d = 0; d < 3; d++)
    CHECK (get_double (file + 40 + 8 * d) == o [d], "origin %d", d);
  CHECK (get_double (file + 64) == H0 && get32 (file + 72) == 0,
         "size, codec");
  uint64_t at = 80;
  const char * want [2] = { "f", "g" };
  for (int k = 0; k < 2; k++) {
    uint32_t l = get32 (file + at);
    CHECK (l == strlen (want [k]) && !memcmp (file + at + 4, want [k], l),
           "name %d", k);
    at += 4 + l;
  }
  uint64_t H = (at + 7) / 8 * 8;
  CHECK (len == H + 16 * N + 2 * N * mm * sizeof (Real), "the file's size");

  /* the leaves, then every cell from the leaf's index alone */
  double vol = 0, worst_f = 0, worst_g = 0;
  uint64_t last = 0;
  int sorted = 1;
  for (uint64_t q = 0; q < N; q++) {
    int32_t v [4];
    for (int c = 0; c < 4; c++)
      v [c] = (int32_t) get32 (file + H + 16 * q + 4 * (uint64_t) c);
    uint64_t kq = key (v);
    if (q && kq <= last)
      sorted = 0;
    last = kq;
    double h = H0 / (double) (1L << v [0]), cell = h / m;
    vol += pow (h, dimension);
    for (uint64_t c = 0; c < mm; c++) {
      double xc [3] = { 0, 0, 0 };
      uint64_t r = c;
      for (int d = 0; d < dimension; d++) {
        xc [d] = o [d] + v [1 + d] * h + ((double) (r % (uint64_t) m) + 0.5)
          * cell;
        r /= (uint64_t) m;
      }
      Real rf, rg;
      memcpy (&rf, file + H + 16 * N + (q * mm + c) * sizeof (Real),
              sizeof (Real));
      memcpy (&rg, file + H + 16 * N + (N * mm + q * mm + c) * sizeof (Real),
              sizeof (Real));
      worst_f = fmax (worst_f, fabs (rf - fun (xc [0], xc [1], xc [2])));
      worst_g = fmax (worst_g, fabs (rg - cell) / cell);
    }
  }
  double tol = sizeof (Real) == sizeof (float) ? 2e-5 : 1e-12;
  CHECK (sorted, "the leaves not along the curve");
  CHECK (fabs (vol - pow (H0, dimension)) < 1e-12 * pow (H0, dimension),
         "the leaves cover %.17g of %.17g", vol, pow (H0, dimension));
  CHECK (worst_f < tol, "a cell off the function by %g", worst_f);
  CHECK (worst_g < tol, "a cell's size off by %g", worst_g);

  /* this process's own part, bit for bit */
  const TreeCaches * tc = tree_caches ();
  uint64_t before = 0;
#ifdef _MPI
  unsigned long long own = (unsigned long long) tc->nowned, b = 0;
  MPI_Exscan (&own, &b, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
  before = tree_rank () ? b : 0;
#endif
  BlockShape sh = tree_block_shape ();
  int bad = 0;
  for (int q = 0; q < tc->nowned; q++) {
    Index p = tc->owned [q].p;
    const unsigned char * e = file + H + 16 * (before + (uint64_t) q);
    int32_t v [4] = { p.level, p.i, 0, 0 };
#if dimension > 1
    v [2] = p.j;
#endif
#if dimension > 2
    v [3] = p.k;
#endif
    for (int c = 0; c < 4; c++)
      bad += get32 (e + 4 * c) != (uint32_t) v [c];
    for (int k = 0; k < 2; k++) {
      const Real * blk = block_data (tc->owned [q].node->s, list [k]);
      for (uint64_t c = 0; c < mm; c++) {
        int x [3] = { 0, 0, 0 };
        uint64_t r = c;
        for (int d = 0; d < dimension; d++)
          x [d] = (int) (r % (uint64_t) m), r /= (uint64_t) m;
        bad += memcmp (file + H + 16 * N + (((uint64_t) k * N + before +
                       (uint64_t) q) * mm + c) * sizeof (Real),
                       &blk [block_at (sh, x)], sizeof (Real)) != 0;
      }
    }
  }
  CHECK (!bad, "%d values of this process's leaves differ", bad);
  free (file);

#ifdef _MPI
  /* a tree not partitioned (whole on every process): written once */
  tree_destroy ();
  tree_replicated (to + 1, 4);
  long nrep = tree_caches ()->nleaves;
  scalar one [1] = { f };
  foreach ()
    f [] = 1;
  CHECK (tree_output (other, one, 1) == 0, "the replicated tree written");
  file = read_all (other, &len);
  CHECK (file && get64 (file + 32) == (uint64_t) nrep &&
         len == 88 + 16 * (uint64_t) nrep + (uint64_t) nrep * mm *
         sizeof (Real), "the replicated tree's %ld leaves once, not %llu",
         nrep, file ? (unsigned long long) get64 (file + 32) : 0ull);
  free (file);
#endif

done:
  tree_destroy ();
  if (tree_rank () == 0) {
    if (fails)
      printf ("%d of %d checks failed\n", fails, checks);
    else
      printf ("%d checks passed, %llu leaves\n", checks,
              (unsigned long long) N);
  }
  return fails != 0;
}
