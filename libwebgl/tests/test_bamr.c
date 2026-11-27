/*
.. webgl_bamr () against libb-amr's writer, 2D: a tree refined over half its
.. domain (moved and resized), global scalars u and w, written by
.. tree_output (), read by webgl_bamr () (libwebgl's reader, not libb-amr's),
.. written as HTML and read back:
..   the boxes the leaves', the levels, the names as in the file, the layer
..   named after the file
..   every cell of u and w, bit for bit, the blocks' Reals in single
..   precision
..   refused, nothing added: no file, and copies broken one way each (the
..   magic, the dimension, the version, a codec, cut short, a level past
..   30, an index out of its level)
.. argv[1] the HTML (default test_bamr_api.html); the .bamr beside it,
.. named test_bamr.bamr, for tests/run.sh's bamr2html check. Prints the
.. HTML's path, then the summary.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "tree.h"
#include "block-ghosts.h"
#include "tree-output.h"
#include "webgl-bamr.h"
#include "decode.h"

scalar u, w;                   /* global: named "u", "w" by cwrap */

static int fails = 0, checks = 0;

#define CHECK(c, ...) do {                                        \
    checks++;                                                     \
    if (!(c)) {                                                   \
      fprintf (stderr, "%s:%d: ", __FILE__, __LINE__);            \
      fprintf (stderr, __VA_ARGS__);                              \
      fputc ('\n', stderr);                                       \
      fails++;                                                    \
    }                                                             \
  } while (0)

static const double X0 = -1.25, Y0 = 2., H0 = 3.;

static void refine_half (int to)
{
  for (int pass = 0; pass < 4 * to; pass++) {
    const TreeCaches * tc = tree_caches ();
    int n = 0;
    Index * todo = malloc ((size_t) tc->nleaves * sizeof (Index));
    for (int q = 0; q < tc->nleaves; q++) {
      Index p = tc->leaves[q].p;
      if (p.level < to && (p.i + 0.5) / (1 << p.level) < 0.5)
        todo[n++] = p;
    }
    for (int a = 0; a < n; a++)
      tree_refine (todo[a]);
    free (todo);
    if (!n)
      break;
  }
}

static void * blob (const char * html, int id, size_t * n, size_t esize)
{
  char name[16];
  snprintf (name, sizeof (name), "b%d", id);
  size_t len = 0;
  void * p = blob_bytes (html, name, &len);
  *n = len / esize;
  return p;
}

static unsigned get32 (const unsigned char * p)
{
  return (unsigned) p[0] | (unsigned) p[1] << 8 | (unsigned) p[2] << 16 |
    (unsigned) p[3] << 24;
}

static void put32 (unsigned char * p, unsigned v)
{
  for (int b = 0; b < 4; b++)
    p[b] = (unsigned char) (v >> (8 * b));
}

int main (int argc, char ** argv)
{
  const char * out = argc > 1 ? argv[1] : "test_bamr_api.html";
  char path[600];
  snprintf (path, sizeof (path), "%s", out);
  char * slash = strrchr (path, '/');
  snprintf (slash ? slash + 1 : path, sizeof (path) - (size_t) (slash ?
            slash + 1 - path : 0), "test_bamr.bamr");

  tree_replicated (4, 4);
  tree_shift_resize (X0, Y0, 0., H0);
  refine_half (3);
  foreach () {
    u[] = sin (2 * x) + x * y;
    w[] = Delta;
  }
  scalar list[2] = { u, w };
  CHECK (tree_output (path, list, 2) == 0, "tree_output (%s)", path);

  /* refused: nothing added */
  WebglScene * S = webgl_scene_new ("test_bamr.bamr");
  CHECK (webgl_bamr (S, "/nonexistent/x.bamr") == -1, "no file");
  size_t len = 0;
  unsigned char * file = (unsigned char *) read_file (path, &len);
  CHECK (file, "read %s", path);
  const char * what[7] = { "the magic", "the dimension", "the version",
                           "a codec", "cut short", "a level past 30",
                           "an index out of its level" };
  const size_t H = (80 + 4 + 1 + 4 + 1 + 7) / 8 * 8;    /* names u, w */
  char broken[620];
  snprintf (broken, sizeof (broken), "%s.broken", path);
  for (int b = 0; file && b < 7; b++) {
    unsigned char * c = malloc (len);
    memcpy (c, file, len);
    size_t n = len;
    switch (b) {
    case 0: c[0] = 'X'; break;
    case 1: put32 (c + 16, 3); break;
    case 2: put32 (c + 8, 2); break;
    case 3: put32 (c + 72, 1); break;
    case 4: n = len - 4; break;
    case 5: put32 (c + H, 40); break;   /* 1 << 40: past an int */
    case 6: put32 (c + H + 4, 1u << get32 (c + H)); break;
    }
    FILE * fp = fopen (broken, "wb");
    if (fp) {
      fwrite (c, 1, n, fp);
      fclose (fp);
    }
    free (c);
    CHECK (webgl_bamr (S, broken) == -1, "refused: %s", what[b]);
  }
  remove (broken);
  free (file);

  /* read, written, read back */
  int layer = webgl_bamr (S, path);
  CHECK (layer == 0, "the layer %d, 0 (nothing added by the refusals)",
         layer);
  CHECK (webgl_scene_write (S, out) == 0, "write %s", out);
  webgl_scene_free (S);
  char * html = read_file (out, NULL);
  CHECK (html, "read %s", out);
  if (!html)
    return 1;
  CHECK (strstr (html, "\"name\":\"test_bamr.bamr\"") &&
         strstr (html, "\"name\":\"u\",\"where\":\"cell\"") &&
         strstr (html, "\"name\":\"w\",\"where\":\"cell\"") &&
         strstr (html, "\"name\":\"level\",\"where\":\"block\","
                 "\"categorical\":1"), "the names");

  const TreeCaches * tc = tree_caches ();
  int nl = tc->nleaves, m = tree_get_block_m ();
  size_t n, mm = (size_t) m * (size_t) m;
  float * box = blob (html, 0, &n, sizeof (float));
  CHECK (n == 4 * (size_t) nl, "%zu box values, %d leaves", n, nl);
  float * fu = blob (html, 1, &n, sizeof (float));
  CHECK (n == (size_t) nl * mm, "%zu values of u", n);
  float * fw = blob (html, 2, &n, sizeof (float));
  float * lev = blob (html, 3, &n, sizeof (float));
  int bad_box = 0, bad_cell = 0;
  BlockShape sh = tree_block_shape ();
  for (int q = 0; box && fu && fw && lev && q < nl; q++) {
    Index p = tc->leaves[q].p;
    double h = H0 / (double) (1L << p.level);
    if (box[4 * q] != (float) (X0 + p.i * h - X0) ||
        box[4 * q + 1] != (float) (Y0 + p.j * h - Y0) ||
        box[4 * q + 2] != (float) h || lev[q] != (float) p.level)
      bad_box++;
    const Real * bu = block_data (tc->leaves[q].node->s, u);
    const Real * bw = block_data (tc->leaves[q].node->s, w);
    for (int j = 0; j < m; j++)
      for (int i = 0; i < m; i++) {
        int x[3] = { i, j, 0 };
        size_t c = (size_t) q * mm + (size_t) (j * m + i);
        if (fu[c] != (float) bu[block_at (sh, x)] ||
            fw[c] != (float) bw[block_at (sh, x)])
          bad_cell++;
      }
  }
  CHECK (!bad_box, "%d boxes or levels off", bad_box);
  CHECK (!bad_cell, "%d cells off", bad_cell);
  free (box), free (fu), free (fw), free (lev), free (html);
  tree_destroy ();
  printf ("%s\n", out);
  if (fails)
    printf ("%d of %d checks failed\n", fails, checks);
  else
    printf ("%d checks passed, %d leaves\n", checks, nl);
  return fails != 0;
}
