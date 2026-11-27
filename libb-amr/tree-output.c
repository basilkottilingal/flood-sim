/*
.. The tree's leaves and scalars to a .bamr file, and back: see
.. tree-output.h.
..
.. The header is written byte by byte, little-endian, not as a struct (no
.. padding to depend on). Without MPI, one process writes in order. With
.. MPI, each process writes its own leaves (tree_caches ()->owned, a run of
.. the curve): its first leaf's place from MPI_Exscan of the counts, the
.. total from MPI_Allreduce, then MPI_File_write_at, in pieces under 1 GiB
.. (counts are ints). Rank 0 writes the header. The file's size is set
.. first (MPI_File_set_size), which also cuts what an older, longer file
.. left. The errors are agreed by one MPI_Allreduce at the end, which also
.. makes every process's writes come before any returns.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#ifdef _MPI
  #include <mpi.h>
  #include "tree-mpi.h"
#endif

#include "tree.h"
#include "block-ghosts.h"
#include "tree-output.h"

enum { HEAD = 80, VERSION = 1 };

static void put32 (unsigned char * p, uint32_t v)
{
  for (int b = 0; b < 4; b++)
    p [b] = (unsigned char) (v >> (8 * b));
}

static void put64 (unsigned char * p, uint64_t v)
{
  for (int b = 0; b < 8; b++)
    p [b] = (unsigned char) (v >> (8 * b));
}

static void put_double (unsigned char * p, double x)
{
  uint64_t v;
  memcpy (&v, &x, sizeof (v));
  put64 (p, v);
}

/* name, with ".bamr" added if its last component has no '.' */
static char * file_name (const char * name)
{
  const char * base = strrchr (name, '/');
  base = base ? base + 1 : name;
  int ext = strchr (base, '.') != NULL;
  size_t n = strlen (name);
  char * s = malloc (n + 6);
  if (!s)
    return NULL;
  memcpy (s, name, n);
  strcpy (s + n, ext ? "" : ".bamr");
  return s;
}

/* 1 if every scalar of the list has a name, none the same as another's:
   a restart finds the fields by their names */
static int names_ok (const scalar * list, int n)
{
  for (int k = 0; k < n; k++) {
    if (!*scalar_name (list [k]))
      return 0;
    for (int l = 0; l < k; l++)
      if (!strcmp (scalar_name (list [k]), scalar_name (list [l])))
        return 0;
  }
  return 1;
}

/* the header: H bytes, H a multiple of 8 */
static unsigned char * header (const scalar * list, int n, uint64_t N,
                               uint64_t * H)
{
  uint64_t h = HEAD;
  for (int k = 0; k < n; k++)
    h += 4 + strlen (scalar_name (list [k]));
  h = (h + 7) / 8 * 8;
  unsigned char * p = calloc (1, h);
  if (!p)
    return NULL;
  memcpy (p, "BAMRGRID", 8);
  put32 (p + 8, VERSION);
  put32 (p + 12, 0x01020304u);
  put32 (p + 16, (uint32_t) dimension);
  put32 (p + 20, (uint32_t) tree_get_block_m ());
  put32 (p + 24, (uint32_t) sizeof (Real));
  put32 (p + 28, (uint32_t) n);
  put64 (p + 32, N);
  for (int d = 0; d < 3; d++)
    put_double (p + 40 + 8 * d, d < dimension ? tree_get_origin (d) : 0.);
  put_double (p + 64, tree_get_size ());
  put32 (p + 72, 0), put32 (p + 76, 0);
  uint64_t at = HEAD;
  for (int k = 0; k < n; k++) {
    const char * s = scalar_name (list [k]);
    uint32_t len = (uint32_t) strlen (s);
    put32 (p + at, len);
    memcpy (p + at + 4, s, len);
    at += 4 + len;
  }
  *H = h;
  return p;
}

typedef struct {
#ifdef _MPI
  MPI_File fh;
#else
  FILE * fp;
#endif
  int failed;
} Out;

/* bytes at offset */
static void write_at (Out * o, uint64_t offset, const void * buf,
                      uint64_t bytes)
{
  if (o->failed || !bytes)
    return;
#ifdef _MPI
  const char * c = buf;
  while (bytes) {
    uint64_t chunk = bytes < (1u << 30) ? bytes : (1u << 30);
    MPI_Status st;
    if (MPI_File_write_at (o->fh, (MPI_Offset) offset, c, (int) chunk,
                           MPI_BYTE, &st) != MPI_SUCCESS) {
      o->failed = 1;
      return;
    }
    c += chunk, offset += chunk, bytes -= chunk;
  }
#else
  if (fseeko (o->fp, (off_t) offset, SEEK_SET) ||
      fwrite (buf, 1, bytes, o->fp) != bytes)
    o->failed = 1;
#endif
}

/* 1 if every process says 0 */
static int all_ok (int failed)
{
#ifdef _MPI
  int any = 0;
  MPI_Allreduce (&failed, &any, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
  return !any;
#else
  return !failed;
#endif
}

int tree_output (const char * name, const scalar * list, int n)
{
  /* arguments: the same on every process, so every one returns -1 */
  int bad = !name || n < 0 || (n > 0 && !list);
  for (int k = 0; !bad && k < n; k++)
    bad = !scalar_in_use (list [k]);
  bad = bad || !names_ok (list, n);
  const uint32_t probe = 1;
  bad = bad || *(const unsigned char *) &probe != 1;   /* not little-endian */
  if (!all_ok (bad))
    return -1;

  const TreeCaches * tc = tree_caches ();
  uint64_t own = (uint64_t) tc->nowned, before = 0, N = own;
#ifdef _MPI
  /* a tree not partitioned yet (tree_replicated ()) is whole on every
     process, each owning every leaf: the first process writes it, once */
  if (!tree_get_partition () && tree_rank () > 0)
    own = 0;
  unsigned long long o = own, b = 0, t = 0;
  N = own;
  MPI_Exscan (&o, &b, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce (&o, &t, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
  before = tree_rank () ? b : 0, N = t;
#endif
  const int m = tree_get_block_m ();
  uint64_t mm = 1;
  for (int d = 0; d < dimension; d++)
    mm *= (uint64_t) m;
  const uint64_t rs = sizeof (Real);

  uint64_t H = 0;
  unsigned char * head = header (list, n, N, &H);
  char * path = file_name (name);
  int failed = !head || !path;
  Out out = { .failed = failed };
#ifdef _MPI
  if (!all_ok (failed)) {
    free (head), free (path);
    return -1;
  }
  out.failed = MPI_File_open (MPI_COMM_WORLD, path,
                              MPI_MODE_CREATE | MPI_MODE_WRONLY,
                              MPI_INFO_NULL, &out.fh) != MPI_SUCCESS;
  if (!all_ok (out.failed)) {
    free (head), free (path);
    return -1;
  }
  out.failed = MPI_File_set_size (out.fh, (MPI_Offset) (H + 16 * N +
                                  (uint64_t) n * N * mm * rs))
    != MPI_SUCCESS;
#else
  if (!failed && !(out.fp = fopen (path, "wb")))
    out.failed = 1;
  if (out.failed) {
    free (head), free (path);
    return -1;
  }
#endif

  /* the header, by the first process */
  if (tree_rank () == 0)
    write_at (&out, 0, head, H);

  /* the own leaves' indices */
  unsigned char * idx = malloc ((size_t) (own ? own : 1) * 16);
  if (!idx)
    out.failed = 1;
  for (uint64_t q = 0; idx && q < own; q++) {
    Index p = tc->owned [q].p;
    int32_t v [4] = { p.level, p.i, 0, 0 };
#if dimension > 1
    v [2] = p.j;
#endif
#if dimension > 2
    v [3] = p.k;
#endif
    for (int c = 0; c < 4; c++)
      put32 (idx + 16 * q + 4 * (uint64_t) c, (uint32_t) v [c]);
  }
  if (idx)
    write_at (&out, H + 16 * before, idx, 16 * own);
  free (idx);

  /* each scalar's cells, the own leaves' */
  Real * v = malloc ((size_t) ((own ? own : 1) * mm) * sizeof (Real));
  if (!v)
    out.failed = 1;
  BlockShape sh = tree_block_shape ();
  for (int k = 0; v && k < n; k++) {
    for (uint64_t q = 0; q < own; q++) {
      const Real * b = block_data (tc->owned [q].node->s, list [k]);
      for (uint64_t c = 0; c < mm; c++) {
        int x [3] = { 0, 0, 0 };
        uint64_t r = c;
        for (int d = 0; d < dimension; d++)
          x [d] = (int) (r % (uint64_t) m), r /= (uint64_t) m;
        v [q * mm + c] = b [block_at (sh, x)];
      }
    }
    write_at (&out, H + 16 * N + (uint64_t) k * N * mm * rs +
              before * mm * rs, v, own * mm * rs);
  }
  free (v), free (head), free (path);

#ifdef _MPI
  if (MPI_File_close (&out.fh) != MPI_SUCCESS)
    out.failed = 1;
#else
  if (fclose (out.fp))
    out.failed = 1;
#endif
  return all_ok (out.failed) ? 0 : -1;
}

/* -------- reading back: tree_restart () -------- */

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

typedef struct {
#ifdef _MPI
  MPI_File fh;
#else
  FILE * fp;
#endif
  uint64_t size;
  int failed;
} In;

static int in_open (In * in, const char * path)
{
  memset (in, 0, sizeof (In));
#ifdef _MPI
  MPI_Offset size = 0;
  in->failed = MPI_File_open (MPI_COMM_WORLD, path, MPI_MODE_RDONLY,
                              MPI_INFO_NULL, &in->fh) != MPI_SUCCESS;
  if (!all_ok (in->failed))
    return -1;
  in->failed = MPI_File_get_size (in->fh, &size) != MPI_SUCCESS;
  in->size = (uint64_t) size;
#else
  if (!(in->fp = fopen (path, "rb")) || fseeko (in->fp, 0, SEEK_END)) {
    if (in->fp)
      fclose (in->fp);
    return -1;
  }
  in->size = (uint64_t) ftello (in->fp);
#endif
  return 0;
}

static void in_close (In * in)
{
#ifdef _MPI
  MPI_File_close (&in->fh);
#else
  fclose (in->fp);
#endif
}

/* bytes at offset, within the file */
static void read_at (In * in, uint64_t offset, void * buf, uint64_t bytes)
{
  if (in->failed || !bytes)
    return;
  if (offset > in->size || bytes > in->size - offset) {
    in->failed = 1;
    return;
  }
#ifdef _MPI
  char * c = buf;
  while (bytes) {
    uint64_t chunk = bytes < (1u << 30) ? bytes : (1u << 30);
    MPI_Status st;
    if (MPI_File_read_at (in->fh, (MPI_Offset) offset, c, (int) chunk,
                          MPI_BYTE, &st) != MPI_SUCCESS) {
      in->failed = 1;
      return;
    }
    c += chunk, offset += chunk, bytes -= chunk;
  }
#else
  if (fseeko (in->fp, (off_t) offset, SEEK_SET) ||
      fread (buf, 1, bytes, in->fp) != bytes)
    in->failed = 1;
#endif
}

/* a leaf's place on the curve: its first cell's key at depth L */
static uint64_t leaf_key (const int32_t * v, int L)
{
  uint64_t k = 0;
  for (int b = 0; b < L; b++)
    for (int d = 0; d < dimension; d++)
      k |= ((uint64_t) v [1 + d] << (L - v [0]) >> b & 1) <<
        (dimension * b + d);
  return k;
}

static Index index_of (const int32_t * v)
{
  Index p;
  memset (&p, 0, sizeof (p));
  p.level = v [0], p.i = v [1];
#if dimension > 1
  p.j = v [2];
#endif
#if dimension > 2
  p.k = v [3];
#endif
  return p;
}

static int same_index (Index a, const int32_t * v)
{
  int same = a.level == v [0] && a.i == v [1];
#if dimension > 1
  same = same && a.j == v [2];
#endif
#if dimension > 2
  same = same && a.k == v [3];
#endif
  return same;
}

/* the file's header and leaves, checked: 0, or -1 (nothing kept) */
typedef struct {
  int m, rs, n;
  uint64_t N, H;
  double origin [3], size;
  char ** names;
  int32_t * leaf;              /* 4 N */
  uint64_t * key;              /* N */
  int depth;                   /* the finest level */
} File;

static void file_free (File * F)
{
  for (int k = 0; F->names && k < F->n; k++)
    free (F->names [k]);
  free (F->names), free (F->leaf), free (F->key);
}

static int read_file (In * in, File * F)
{
  memset (F, 0, sizeof (File));
  unsigned char h [HEAD];
  read_at (in, 0, h, HEAD);
  if (in->failed || memcmp (h, "BAMRGRID", 8) || get32 (h + 8) != VERSION ||
      get32 (h + 12) != 0x01020304u || get32 (h + 16) != dimension ||
      get32 (h + 20) != (uint32_t) tree_get_block_m () ||
      (get32 (h + 24) != 4 && get32 (h + 24) != 8) || get32 (h + 72) != 0 ||
      get32 (h + 28) > 4096)
    return -1;
  F->m = (int) get32 (h + 20), F->rs = (int) get32 (h + 24);
  F->n = (int) get32 (h + 28), F->N = get64 (h + 32);
  for (int d = 0; d < 3; d++)
    F->origin [d] = get_double (h + 40 + 8 * d);
  F->size = get_double (h + 64);
  if (F->N < 1 || F->N > (1ull << 40) || !(F->size > 0.))
    return -1;
  /* the names */
  F->names = calloc ((size_t) (F->n ? F->n : 1), sizeof (char *));
  uint64_t at = HEAD;
  for (int k = 0; F->names && k < F->n; k++) {
    unsigned char b [4];
    read_at (in, at, b, 4);
    uint32_t len = get32 (b);
    if (in->failed || len < 1 || len > 4096 || !(F->names [k] =
                                                  malloc (len + 1)))
      return -1;
    read_at (in, at + 4, F->names [k], len);
    F->names [k] [len] = '\0';
    at += 4 + len;
  }
  if (!F->names || in->failed)
    return -1;
  F->H = (at + 7) / 8 * 8;
  uint64_t mm = 1;
  for (int d = 0; d < dimension; d++)
    mm *= (uint64_t) F->m;
  if (in->size != F->H + 16 * F->N + (uint64_t) F->n * F->N * mm *
      (uint64_t) F->rs)
    return -1;
  /* the leaves: in range, along the curve, covering the domain */
  unsigned char * raw = malloc (16 * F->N);
  F->leaf = malloc (16 * F->N);
  F->key = malloc (8 * F->N);
  if (!raw || !F->leaf || !F->key) {
    free (raw);
    return -1;
  }
  read_at (in, F->H, raw, 16 * F->N);
  int bad = in->failed, depthmax = tree_get_depth_max ();
  for (uint64_t q = 0; !bad && q < F->N; q++)
    for (int c = 0; c < 4; c++) {
      int32_t v = (int32_t) get32 (raw + 16 * q + 4 * (uint64_t) c);
      F->leaf [4 * q + (uint64_t) c] = v;
      if (c == 0)
        bad |= v < 0 || v > depthmax;
      else                     /* its level checked first: in 0 ... 31 */
        bad |= v < 0 || (c > dimension ? v != 0 :
                         bad || v >= (1 << F->leaf [4 * q]));
    }
  free (raw);
  if (bad)
    return -1;
  int L = 0;
  for (uint64_t q = 0; q < F->N; q++)
    L = F->leaf [4 * q] > L ? F->leaf [4 * q] : L;
  F->depth = L;
  uint64_t cover = 0;
  for (uint64_t q = 0; q < F->N; q++) {
    F->key [q] = leaf_key (F->leaf + 4 * q, L);
    if (q && F->key [q] <= F->key [q - 1])
      return -1;
    cover += 1ull << (dimension * (L - F->leaf [4 * q]));
  }
  return cover == 1ull << (dimension * L) ? 0 : -1;
}

/* the file's leaves made: from level 1, each leaf's ancestor at each
   level refined, coarser levels first (the file's tree is 2:1 balanced, so
   the balancing refinements stay within it) */
static void rebuild (const File * F)
{
  int blockdepth = 0;
  while ((1 << blockdepth) < F->m)
    blockdepth++;
  int depthmax = tree_get_depth_max ();
  tree_destroy ();
  tree_replicated (depthmax, blockdepth);
  tree_shift_resize (F->origin [0], F->origin [1], F->origin [2], F->size);
  if (F->depth == 0) {               /* one leaf: the root */
    Index root;
    memset (&root, 0, sizeof (root));
    tree_coarsen (root);
    return;
  }
  for (int l = 1; l < F->depth; l++) {
    int32_t prev [4] = { -1, -1, -1, -1 };
    for (uint64_t q = 0; q < F->N; q++) {
      const int32_t * v = F->leaf + 4 * q;
      if (v [0] <= l)
        continue;
      int32_t a [4] = { l, v [1] >> (v [0] - l), v [2] >> (v [0] - l),
                        v [3] >> (v [0] - l) };
      if (memcmp (a, prev, sizeof (a))) {
        tree_refine (index_of (a));
        memcpy (prev, a, sizeof (a));
      }
    }
  }
}

/* 1 if a tree is there: leaves on some process (a process may hold none
   of a distributed tree) */
static int tree_there (void)
{
  int here = tree_get_nleaves () > 0, any = here;
#ifdef _MPI
  MPI_Allreduce (&here, &any, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
#endif
  return any;
}

int tree_restart (const char * path)
{
  if (!tree_there () || !all_ok (!path))
    return -1;
  In in;
  if (in_open (&in, path))
    return -1;
  File F;
  int failed = read_file (&in, &F) || in.failed;
  if (!all_ok (failed)) {
    file_free (&F);
    in_close (&in);
    return -1;
  }

  /* past this point the program's tree is replaced */
  rebuild (&F);
  const TreeCaches * tc = tree_caches ();
  failed = (uint64_t) tc->nleaves != F.N;
  for (uint64_t q = 0; !failed && q < F.N; q++)
    failed = !same_index (tc->leaves [q].p, F.leaf + 4 * q);
#ifdef _MPI
  if (all_ok (failed))
    failed = tree_mpi_partition (NULL) != 0;
  tc = tree_caches ();
#endif
  if (!all_ok (failed)) {
    file_free (&F);
    in_close (&in);
    return -1;
  }

  /* this process's leaves: a run of the file's, found by its first key */
  uint64_t own = (uint64_t) tc->nowned, first = 0;
  if (own) {
    int32_t v [4] = { tc->owned [0].p.level, tc->owned [0].p.i, 0, 0 };
#if dimension > 1
    v [2] = tc->owned [0].p.j;
#endif
#if dimension > 2
    v [3] = tc->owned [0].p.k;
#endif
    uint64_t k = leaf_key (v, F.depth), lo = 0, hi = F.N;
    while (lo < hi) {
      uint64_t c = lo + (hi - lo) / 2;
      if (F.key [c] < k)
        lo = c + 1;
      else
        hi = c;
    }
    first = lo;
    for (uint64_t q = 0; !failed && q < own; q++)
      failed = first + q >= F.N ||
        !same_index (tc->owned [q].p, F.leaf + 4 * (first + q));
  }

  /* each global scalar named in the file, its cells */
  const int m = F.m;
  uint64_t mm = 1;
  for (int d = 0; d < dimension; d++)
    mm *= (uint64_t) m;
  const uint64_t rs = (uint64_t) F.rs;
  unsigned char * raw = malloc ((size_t) ((own ? own : 1) * mm * rs));
  failed |= !raw;
  BlockShape sh = tree_block_shape ();
  int restored = 0;
  for (int k = 0; !failed && k < F.n; k++) {
    scalar s = scalar_find (F.names [k]);
    if (s < 0 || s >= n_global_scalars || !scalar_in_use (s))
      continue;                        /* not a global scalar: left */
    read_at (&in, F.H + 16 * F.N + ((uint64_t) k * F.N + first) * mm * rs,
             raw, own * mm * rs);
    for (uint64_t q = 0; !in.failed && q < own; q++) {
      Real * b = block_data (tc->owned [q].node->s, s);
      for (uint64_t c = 0; c < mm; c++) {
        int x [3] = { 0, 0, 0 };
        uint64_t r = c;
        for (int d = 0; d < dimension; d++)
          x [d] = (int) (r % (uint64_t) m), r /= (uint64_t) m;
        const unsigned char * e = raw + (q * mm + c) * rs;
        if (rs == 8) {
          double y;
          memcpy (&y, e, 8);
          b [block_at (sh, x)] = (Real) y;
        }
        else {
          float y;
          memcpy (&y, e, 4);
          b [block_at (sh, x)] = (Real) y;
        }
      }
    }
    ghosts_mark_dirty (s);
    restored++;
  }
  failed |= in.failed;
  free (raw);
  file_free (&F);
  in_close (&in);
  return all_ok (failed) ? restored : -1;
}
