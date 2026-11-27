/*
.. libwebgl's scene, written and read back:
..   every kind of layer (blocks, lines, triangles, points), far from the
..   origin, with fields of every kind (cells with NaN and inf, blocks
..   categorical, points, elements)
..   the calls refused: a field of the wrong kind for its layer, an index
..   out of range, NULLs, a layer out of range; nothing added by them
..   the HTML: the marks all filled, the title and names escaped (a name
..   holding </script> and quotes), as many <script> opened as closed
..   every blob decoded (a decoder of the test's own) and compared bit for
..   bit with what it must hold: offsets from the scene's corner in single
..   precision, indices, field values
..   three.js and OrbitControls in the file byte for byte (a build with
..   them), or the import map pointing at jsDelivr's (a build without);
..   the viewer in the file, byte for byte
.. Prints the HTML's path, then the summary, for tests/run.sh.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "webgl.h"
#include "webgl-assets.h"
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

/* the blob's bytes against the expected, bit for bit */
static void check_blob (const char * html, int id, const void * want,
                        size_t bytes, const char * what)
{
  char name[16];
  snprintf (name, sizeof (name), "b%d", id);
  size_t len = 0;
  unsigned char * got = blob_bytes (html, name, &len);
  CHECK (got && len == bytes && !memcmp (got, want, bytes),
         "blob %s (%s): %zu bytes, %zu wanted, %s", name, what, len, bytes,
         got ? (len == bytes ? "different" : "wrong size") : "not found");
  free (got);
}

/* the text of the import map's entry, decoded, against the asset */
static void check_import (const char * html, const char * key,
                          const unsigned char * want, size_t n)
{
  const char * a = strstr (html, key);
  CHECK (a, "import map: %s missing", key);
  if (!a)
    return;
  a = strstr (a, "base64,");
  const char * b = a ? strchr (a, '"') : NULL;
  CHECK (a && b, "import map: %s not a data: URL", key);
  if (!a || !b)
    return;
  a += 7;
  size_t len = 0;
  unsigned char * got = b64_decode (a, (size_t) (b - a), &len);
  CHECK (got && len == n && !memcmp (got, want, n),
         "import map: %s not the vendored file", key);
  free (got);
}

int main (int argc, char ** argv)
{
  const char * out = argc > 1 ? argv[1] : "test_webgl.html";
  WebglScene * S = webgl_scene_new ("a <b>\"scene\"</b> & more");

  /* blocks: three of 4 x 4, far from the origin */
  const int n = 3, m = 4;
  double box[12] = { 1000.25, -500.5, 2., 2.,
                     1002.25, -500.5, 2., 2.,
                     1000.25, -498.5, 1., 1. };
  int lb = webgl_blocks (S, "grid </script> \"q\"", n, m, box);
  CHECK (lb == 0, "blocks layer %d", lb);
  float cell[48], blk[3] = { 3, 4, 5 };
  for (int i = 0; i < 48; i++)
    cell[i] = (float) (i * 0.37 - 2.);
  cell[5] = NAN, cell[17] = INFINITY;
  CHECK (webgl_field (S, lb, "u", WEBGL_CELL, cell, 0) == 0, "cell field");
  CHECK (webgl_field (S, lb, "level", WEBGL_BLOCK, blk, 1) == 0,
         "block field");

  /* lines: a curve of 5 points */
  double lxy[10] = { 1000.5, -499., 1001., -498.75, 1001.5, -498.6,
                     1002., -498.75, 1003.75, -497.5 };
  uint32_t le[8] = { 0, 1, 1, 2, 2, 3, 3, 4 };
  int ll = webgl_lines (S, "curve", 5, lxy, 4, le);
  CHECK (ll == 1, "lines layer %d", ll);
  float lp[5] = { 0, 1, 2, 3, 4 }, ls[4] = { 0.5f, -1.5f, 2.5f, NAN };
  CHECK (webgl_field (S, ll, "h", WEBGL_POINT, lp, 0) == 0, "point field");
  CHECK (webgl_field (S, ll, "q", WEBGL_ELEMENT, ls, 0) == 0,
         "element field");

  /* triangles: two over a square */
  double txy[8] = { 1000.25, -500.5, 1001.25, -500.5, 1001.25, -499.5,
                    1000.25, -499.5 };
  uint32_t tt[6] = { 0, 1, 2, 0, 2, 3 };
  int lt = webgl_triangles (S, "patch", 4, txy, 2, tt);
  CHECK (lt == 2, "triangles layer %d", lt);
  float tf[2] = { 1e-7f, 7e6f };
  CHECK (webgl_field (S, lt, "a", WEBGL_ELEMENT, tf, 0) == 0,
         "triangle field");

  /* points */
  double pxy[6] = { 1004.25, -497.25, 999.75, -501., 1001., -499. };
  int lpnt = webgl_points (S, "dots", 3, pxy);
  CHECK (lpnt == 3, "points layer %d", lpnt);
  float pf[3] = { -1, 0, 1 };
  CHECK (webgl_field (S, lpnt, "r", WEBGL_POINT, pf, 1) == 0, "points field");

  /* refused, nothing added */
  CHECK (webgl_field (S, ll, "x", WEBGL_CELL, cell, 0) == -1,
         "a cell field on lines");
  CHECK (webgl_field (S, lb, "x", WEBGL_POINT, lp, 0) == -1,
         "a point field on blocks");
  CHECK (webgl_field (S, lpnt, "x", WEBGL_ELEMENT, pf, 0) == -1,
         "an element field on points");
  CHECK (webgl_field (S, 9, "x", WEBGL_POINT, lp, 0) == -1, "layer 9");
  CHECK (webgl_field (S, -1, "x", WEBGL_POINT, lp, 0) == -1, "layer -1");
  CHECK (webgl_field (S, ll, NULL, WEBGL_POINT, lp, 0) == -1, "no name");
  CHECK (webgl_field (S, ll, "x", WEBGL_POINT, NULL, 0) == -1, "no values");
  CHECK (webgl_field (S, ll, "x", 7, lp, 0) == -1, "where 7");
  uint32_t bad[2] = { 0, 5 };
  CHECK (webgl_lines (S, "bad", 5, lxy, 1, bad) == -1, "index 5 of 5");
  uint32_t badt[3] = { 0, 1, 4 };
  CHECK (webgl_triangles (S, "bad", 4, txy, 1, badt) == -1, "index 4 of 4");
  CHECK (webgl_blocks (S, "bad", -1, 4, box) == -1, "n < 0");
  CHECK (webgl_blocks (S, "bad", 1, 0, box) == -1, "m = 0");
  CHECK (webgl_points (S, NULL, 1, pxy) == -1, "no name");
  CHECK (webgl_points (S, "e", 0, NULL) == 4, "an empty layer is fine");
  CHECK (webgl_scene_write (S, NULL) == -1, "no path");
  CHECK (webgl_scene_write (S, "/nonexistent/dir/x.html") == -1,
         "a path that cannot be written");

  CHECK (webgl_scene_write (S, out) == 0, "write %s", out);
  webgl_scene_free (S);

  size_t len;
  char * html = read_file (out, &len);
  CHECK (html, "read %s back", out);
  if (!html)
    goto done;

  /* the marks filled, the text escaped, scripts balanced */
  for (const char * mk = "@TITLE@\0@IMPORTMAP@\0@DATA@\0@VIEWER@\0"; *mk;
       mk += strlen (mk) + 1)
    CHECK (!strstr (html, mk), "mark %s left in the file", mk);
  CHECK (strstr (html, "<title>a &lt;b&gt;&quot;scene&quot;&lt;/b&gt; "
                 "&amp; more</title>"), "the title escaped");
  CHECK (strstr (html, "\"grid \\u003c/script\\u003e \\\"q\\\"\""),
         "the name escaped in the manifest");
  int open = 0, close = 0;
  for (const char * p = html; (p = strstr (p, "<script")); p++)
    open++;
  for (const char * p = html; (p = strstr (p, "</script>")); p++)
    close++;
  CHECK (open == close, "%d <script> opened, %d closed", open, close);

  /* the corner: the smallest coordinate of all */
  const double lo[2] = { 999.75, -501. };

  /* blocks: b0 the boxes, b1 the cells, b2 the blocks' field */
  float fbox[12];
  for (int i = 0; i < 12; i++)
    fbox[i] = (float) (box[i] - ((i % 4) < 2 ? lo[i % 4] : 0.));
  check_blob (html, 0, fbox, sizeof (fbox), "boxes");
  check_blob (html, 1, cell, sizeof (cell), "cells");
  check_blob (html, 2, blk, sizeof (blk), "blocks");
  /* lines: b3 positions, b4 segments, b5 h, b6 q */
  float fl[10];
  for (int i = 0; i < 10; i++)
    fl[i] = (float) (lxy[i] - lo[i % 2]);
  check_blob (html, 3, fl, sizeof (fl), "line points");
  check_blob (html, 4, le, sizeof (le), "segments");
  check_blob (html, 5, lp, sizeof (lp), "h");
  check_blob (html, 6, ls, sizeof (ls), "q");
  /* triangles: b7, b8, b9 */
  float ft[8];
  for (int i = 0; i < 8; i++)
    ft[i] = (float) (txy[i] - lo[i % 2]);
  check_blob (html, 7, ft, sizeof (ft), "triangle points");
  check_blob (html, 8, tt, sizeof (tt), "triangles");
  check_blob (html, 9, tf, sizeof (tf), "a");
  /* points: b10, b11; the empty layer: b12 */
  float fp[6];
  for (int i = 0; i < 6; i++)
    fp[i] = (float) (pxy[i] - lo[i % 2]);
  check_blob (html, 10, fp, sizeof (fp), "points");
  check_blob (html, 11, pf, sizeof (pf), "r");
  size_t elen = 1;
  unsigned char * empty = blob_bytes (html, "b12", &elen);
  CHECK (empty && elen == 0, "the empty layer's blob");
  free (empty);
  CHECK (!strstr (html, "id=\"b13\""), "no blob more");

  /* the code in the file: three.js embedded, or its URLs */
  if (webgl_offline ()) {
    check_import (html, "\"three\":", vendor_three_module_min_js,
                  vendor_three_module_min_js_len);
    check_import (html, "\"three/addons/controls/OrbitControls.js\":",
                  vendor_OrbitControls_js, vendor_OrbitControls_js_len);
  }
  else {
    CHECK (strstr (html, "\"three\":\"" WEBGL_THREE_CDN
                   "build/three.module.min.js\"") &&
           strstr (html, "\"three/addons/controls/OrbitControls.js\":\""
                   WEBGL_THREE_CDN "examples/jsm/controls/OrbitControls.js\""),
           "online: the import map's URLs");
    CHECK (!strstr (html, "data:text/javascript"), "online: nothing embedded");
  }
  const char * js = strstr (html, "<script type=\"module\">\n");
  CHECK (js && !memcmp (js + 23, viewer_viewer_js, viewer_viewer_js_len),
         "the viewer inline");
  free (html);

done:
  printf ("%s\n", out);
  if (fails)
    printf ("%d of %d checks failed\n", fails, checks);
  else
    printf ("%d checks passed (three.js %s)\n", checks,
            webgl_offline () ? "embedded" : "from jsDelivr");
  return fails != 0;
}
