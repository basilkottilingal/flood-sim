/*
.. libwebgl's adapters, 2D: a tree refined over half its domain (moved and
.. resized), a scalar on it, and a circle from libhmesh with fields on its
.. points and edges, written and read back:
..   the blocks: each leaf's box as its Index gives it; each cell's value
..   the function at the cell's centre (so the cells' order, i fastest
..   then j, is right); levels and ranks
..   the mesh, a few edges removed (holes in the ids): its points where
..   hm_x () has them, each segment between the
..   points hm_vertex () names (by the ids written with them), lengths as
..   hm_measure () gives, a scalar on the points, a vector on the edges as
..   two fields
.. Prints the HTML's path, then the summary.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "tree.h"
#include "block-ghosts.h"
#include "hmesh-tree.h"
#include "webgl-tree.h"
#include "webgl-hmesh.h"
#include "decode.h"

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

static double fu (double x, double y)
{
  return sin (3 * x) * cos (2 * y) + x;
}

static double circle (double x, double y)
{
  double dx = x - (X0 + 0.5 * H0), dy = y - (Y0 + 0.5 * H0);
  return sqrt (dx * dx + dy * dy) - 0.3 * H0;
}

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

/* blob id decoded as floats or ints, its count in *n */
static void * blob (const char * html, int id, size_t * n, size_t esize)
{
  char name[16];
  snprintf (name, sizeof (name), "b%d", id);
  size_t len = 0;
  void * p = blob_bytes (html, name, &len);
  *n = len / esize;
  return p;
}

int main (int argc, char ** argv)
{
  const char * out = argc > 1 ? argv[1] : "test_scene.html";
  tree (4, 4);
  tree_shift_resize (X0, Y0, 0., H0);
  refine_half (3);
  scalar u = scalar_new ();
  foreach ()
    u[] = fu (x, y);

  Hmesh * M = hmesh_new (1);
  CHECK (hm_tree_from_function (M, circle) > 0, "the circle");
  /* a few edges out, with their points: holes in the ids, so the points
     must be renumbered */
  for (uint32_t v = 5; v < 60; v += 23) {
    uint32_t at[8], na = 0;
    hm_foreach (M, 1, e)
      if (na < 8 && (hm_vertex (M, 1, e, 0) == v ||
                     hm_vertex (M, 1, e, 1) == v))
        at[na++] = e;
    for (uint32_t i = 0; i < na; i++)
      hm_remove_closure (M, 1, at[i]);
  }
  CHECK (hm_count (M, 0) < hm_end (M, 0), "no hole in the points' ids");
  int fh = hm_scalar_new (M, 0, "height");
  int fq = hm_vector_new (M, 1, "flux");
  hm_foreach (M, 0, v)
    hm_val (M, fh, 0)[v] = (HmReal) hm_x (M, 1)[v];
  hm_foreach (M, 1, e) {
    hm_val (M, fq, 0)[e] = (HmReal) e;
    hm_val (M, fq, 1)[e] = (HmReal) (-0.5 * e);
  }

  WebglScene * S = webgl_scene_new ("tree and interface");
  int lg = webgl_tree (S, "grid", 1, &u);
  int lm = webgl_hmesh (S, "interface", M);
  CHECK (lg == 0 && lm == 1, "layers %d %d", lg, lm);
  CHECK (webgl_scene_write (S, out) == 0, "write %s", out);
  webgl_scene_free (S);
  char * html = read_file (out, NULL);
  CHECK (html, "read back");
  if (!html)
    return 1;

  /* the corner is the tree's origin: everything lies in the domain */
  const double lo[2] = { X0, Y0 };
  const TreeCaches * tc = tree_caches ();
  int nl = tc->nleaves, m = tree_get_block_m ();
  size_t n, mm = (size_t) m * (size_t) m;

  /* b0 boxes, b1 u, b2 level, b3 pid */
  float * box = blob (html, 0, &n, sizeof (float));
  CHECK (n == 4 * (size_t) nl, "%zu box values for %d leaves", n, nl);
  float * fv = blob (html, 1, &n, sizeof (float));
  CHECK (n == (size_t) nl * mm, "%zu cell values", n);
  float * lev = blob (html, 2, &n, sizeof (float));
  float * pid = blob (html, 3, &n, sizeof (float));
  double worst = 0;
  int bad_box = 0, bad_level = 0;
  for (int q = 0; box && fv && lev && pid && q < nl; q++) {
    Index p = tc->leaves[q].p;
    double h = H0 / (double) (1L << p.level);
    if (box[4 * q] != (float) (X0 + p.i * h - lo[0]) ||
        box[4 * q + 1] != (float) (Y0 + p.j * h - lo[1]) ||
        box[4 * q + 2] != (float) h || box[4 * q + 3] != (float) h)
      bad_box++;
    if (lev[q] != (float) p.level || pid[q] != (float) tc->leaves[q].node->pid)
      bad_level++;
    for (int j = 0; j < m; j++)
      for (int i = 0; i < m; i++) {
        double xc = X0 + p.i * h + (i + 0.5) * h / m;
        double yc = Y0 + p.j * h + (j + 0.5) * h / m;
        worst = fmax (worst, fabs (fv[(size_t) q * mm + (size_t) (j * m + i)]
                                   - fu (xc, yc)));
      }
  }
  CHECK (!bad_box, "%d boxes off", bad_box);
  CHECK (!bad_level, "%d levels or ranks off", bad_level);
  CHECK (worst < 1e-5, "a cell off its centre's value by %g", worst);
  free (box), free (fv), free (lev), free (pid);

  /* b4 xy, b5 segments, b6 point ids, b7 height, b8 edge ids, b9 length,
     b10 flux.x, b11 flux.y */
  uint32_t np = hm_count (M, 0), ne = hm_count (M, 1);
  float * xy = blob (html, 4, &n, sizeof (float));
  CHECK (n == 2 * (size_t) np, "%zu coordinates for %u points", n, np);
  uint32_t * seg = blob (html, 5, &n, sizeof (uint32_t));
  CHECK (n == 2 * (size_t) ne, "%zu indices for %u edges", n, ne);
  float * vid = blob (html, 6, &n, sizeof (float));
  float * ht = blob (html, 7, &n, sizeof (float));
  float * eid = blob (html, 8, &n, sizeof (float));
  float * len = blob (html, 9, &n, sizeof (float));
  float * qx = blob (html, 10, &n, sizeof (float));
  float * qy = blob (html, 11, &n, sizeof (float));
  int bad_pt = 0, bad_seg = 0;
  for (uint32_t a = 0; xy && vid && ht && a < np; a++) {
    uint32_t v = (uint32_t) vid[a];
    if (!hm_alive (M, 0, v) || xy[2 * a] != (float) (hm_x (M, 0)[v] - lo[0])
        || xy[2 * a + 1] != (float) (hm_x (M, 1)[v] - lo[1])
        || ht[a] != (float) hm_val (M, fh, 0)[v])
      bad_pt++;
  }
  for (uint32_t b = 0; seg && eid && len && qx && qy && vid && b < ne; b++) {
    uint32_t e = (uint32_t) eid[b];
    if (!hm_alive (M, 1, e) ||
        (uint32_t) vid[seg[2 * b]] != hm_vertex (M, 1, e, 0) ||
        (uint32_t) vid[seg[2 * b + 1]] != hm_vertex (M, 1, e, 1) ||
        len[b] != (float) hm_measure (M, 1, e) ||
        qx[b] != (float) hm_val (M, fq, 0)[e] ||
        qy[b] != (float) hm_val (M, fq, 1)[e])
      bad_seg++;
  }
  CHECK (!bad_pt, "%d points off", bad_pt);
  CHECK (!bad_seg, "%d segments off", bad_seg);
  CHECK (strstr (html, "\"name\":\"flux.y\""),
         "the vector's components named");
  CHECK (strstr (html, "\"name\":\"level\",\"where\":\"block\","
                 "\"categorical\":1"), "levels categorical");
  free (xy), free (seg), free (vid), free (ht), free (eid), free (len);
  free (qx), free (qy), free (html);

  printf ("%s: %d leaves of %d^2 cells, %u points, %u edges\n", out, nl, m,
          np, ne);
  hmesh_free (M);
  scalar_free (u);
  tree_destroy ();
  printf ("%s\n", out);
  if (fails)
    printf ("%d of %d checks failed\n", fails, checks);
  else
    printf ("%d checks passed\n", checks);
  return fails != 0;
}
