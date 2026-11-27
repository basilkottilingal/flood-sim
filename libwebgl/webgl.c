/*
.. A scene written as one HTML file: see webgl.h.
..
.. The file is the viewer's template (viewer/viewer.html) with four marks
.. filled in:
..
..   @TITLE@      the title, escaped for HTML
..   @IMPORTMAP@  three.js and OrbitControls as data: URLs (base64), so that
..                the viewer's imports resolve with no network and no file
..                beside the page; or, when the build has no three.js
..                (vendor/ empty), jsDelivr's URLs of the same version
..   @DATA@       the manifest (JSON: the layers, their fields, the range of
..                each) and one <script type="application/octet-stream"> a
..                typed array, base64, little-endian: f32 for boxes,
..                positions and fields, u32 for indices
..   @VIEWER@     viewer/viewer.js, inline as a module
..
.. The template, the viewer and three.js are compiled in (webgl-assets.c,
.. made by xxd -i), so a program needs no file at run time.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "webgl.h"
#include "webgl-assets.h"

enum { BLOCKS, LINES, TRIANGLES, POINTS };

typedef struct {
  char  * name;
  int     where, categorical;
  float * v;
  size_t  n;
} Field;

typedef struct {
  int        type;
  char     * name;
  int        n, m;            /* blocks: their number, cells a side */
  double   * box;             /* blocks: 4 n */
  int        np, ne;          /* points, segments or triangles */
  double   * xy;              /* 2 np */
  uint32_t * e;               /* 2 ne (lines) or 3 ne (triangles) */
  Field    * field;
  int        nfield;
} Layer;

struct WebglScene {
  char  * title;
  Layer * layer;
  int     nlayer;
};

static void * xmalloc (size_t n)
{
  void * p = malloc (n ? n : 1);
  if (!p) {
    fprintf (stderr, "webgl: out of memory\n");
    abort ();
  }
  return p;
}

static void * copy (const void * p, size_t n)
{
  void * q = xmalloc (n);
  if (n)
    memcpy (q, p, n);
  return q;
}

static char * copy_string (const char * s)
{
  return copy (s, strlen (s) + 1);
}

WebglScene * webgl_scene_new (const char * title)
{
  WebglScene * S = xmalloc (sizeof (WebglScene));
  memset (S, 0, sizeof (WebglScene));
  S->title = copy_string (title ? title : "");
  return S;
}

void webgl_scene_free (WebglScene * S)
{
  if (!S)
    return;
  for (int l = 0; l < S->nlayer; l++) {
    Layer * L = &S->layer[l];
    for (int f = 0; f < L->nfield; f++)
      free (L->field[f].name), free (L->field[f].v);
    free (L->field), free (L->name), free (L->box), free (L->xy);
    free (L->e);
  }
  free (S->layer), free (S->title), free (S);
}

static Layer * new_layer (WebglScene * S, int type, const char * name)
{
  Layer * l = realloc (S->layer, (size_t) (S->nlayer + 1) * sizeof (Layer));
  if (!l) {
    fprintf (stderr, "webgl: out of memory\n");
    abort ();
  }
  S->layer = l;
  Layer * L = &S->layer[S->nlayer++];
  memset (L, 0, sizeof (Layer));
  L->type = type;
  L->name = copy_string (name);
  return L;
}

int webgl_blocks (WebglScene * S, const char * name, int n, int m,
                  const double * box)
{
  if (!S || !name || n < 0 || m < 1 || (n && !box))
    return -1;
  Layer * L = new_layer (S, BLOCKS, name);
  L->n = n, L->m = m;
  L->box = copy (box, 4 * (size_t) n * sizeof (double));
  return S->nlayer - 1;
}

/* a layer of points and elements of k points each (0: points alone) */
static int geometry (WebglScene * S, int type, const char * name, int np,
                     const double * xy, int ne, const uint32_t * e, int k)
{
  if (!S || !name || np < 0 || ne < 0 || (np && !xy) || (ne && !e))
    return -1;
  for (size_t i = 0; i < (size_t) ne * (size_t) k; i++)
    if (e[i] >= (uint32_t) np)
      return -1;
  Layer * L = new_layer (S, type, name);
  L->np = np, L->ne = ne;
  L->xy = copy (xy, 2 * (size_t) np * sizeof (double));
  L->e = copy (e, (size_t) ne * (size_t) k * sizeof (uint32_t));
  return S->nlayer - 1;
}

int webgl_lines (WebglScene * S, const char * name, int np, const double * xy,
                 int ne, const uint32_t * e)
{
  return geometry (S, LINES, name, np, xy, ne, e, 2);
}

int webgl_triangles (WebglScene * S, const char * name, int np,
                     const double * xy, int ne, const uint32_t * t)
{
  return geometry (S, TRIANGLES, name, np, xy, ne, t, 3);
}

int webgl_points (WebglScene * S, const char * name, int np, const double * xy)
{
  return geometry (S, POINTS, name, np, xy, 0, NULL, 0);
}

/* may a field live there on this layer: cells and blocks on blocks, points
   on the others, elements on lines and triangles */
static int field_allowed (const Layer * L, int where)
{
  switch (where) {
  case WEBGL_CELL: case WEBGL_BLOCK: return L->type == BLOCKS;
  case WEBGL_POINT: return L->type != BLOCKS;
  case WEBGL_ELEMENT: return L->type == LINES || L->type == TRIANGLES;
  }
  return 0;
}

/* the number of values of a field allowed there */
static size_t field_size (const Layer * L, int where)
{
  switch (where) {
  case WEBGL_CELL: return (size_t) L->n * (size_t) L->m * (size_t) L->m;
  case WEBGL_BLOCK: return (size_t) L->n;
  case WEBGL_POINT: return (size_t) L->np;
  default: return (size_t) L->ne;
  }
}

int webgl_field (WebglScene * S, int layer, const char * name, int where,
                 const float * v, int categorical)
{
  if (!S || !name || layer < 0 || layer >= S->nlayer)
    return -1;
  Layer * L = &S->layer[layer];
  if (!field_allowed (L, where))
    return -1;
  size_t n = field_size (L, where);
  if (n && !v)
    return -1;
  Field * f = realloc (L->field, (size_t) (L->nfield + 1) * sizeof (Field));
  if (!f) {
    fprintf (stderr, "webgl: out of memory\n");
    abort ();
  }
  L->field = f;
  Field * F = &L->field[L->nfield++];
  F->name = copy_string (name);
  F->where = where, F->categorical = categorical != 0;
  F->n = n;
  F->v = copy (v, n * sizeof (float));
  return 0;
}

/* -------- writing -------- */

static const char b64 [] =
  "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";

static void base64 (FILE * fp, const unsigned char * p, size_t n)
{
  size_t i = 0;
  for (; i + 2 < n; i += 3) {
    uint32_t w = (uint32_t) p[i] << 16 | (uint32_t) p[i + 1] << 8 | p[i + 2];
    fputc (b64[w >> 18], fp), fputc (b64[w >> 12 & 63], fp);
    fputc (b64[w >> 6 & 63], fp), fputc (b64[w & 63], fp);
  }
  if (i + 1 == n) {
    uint32_t w = (uint32_t) p[i] << 16;
    fputc (b64[w >> 18], fp), fputc (b64[w >> 12 & 63], fp);
    fputs ("==", fp);
  }
  else if (i + 2 == n) {
    uint32_t w = (uint32_t) p[i] << 16 | (uint32_t) p[i + 1] << 8;
    fputc (b64[w >> 18], fp), fputc (b64[w >> 12 & 63], fp);
    fputc (b64[w >> 6 & 63], fp), fputc ('=', fp);
  }
}

/* a JSON string: quotes, backslashes and controls escaped, and < too, so
   that no name ends the <script> it is in */
static void json_string (FILE * fp, const char * s)
{
  fputc ('"', fp);
  for (; *s; s++) {
    unsigned char c = (unsigned char) *s;
    if (c == '"' || c == '\\')
      fprintf (fp, "\\%c", c);
    else if (c < 0x20 || c == '<' || c == '>' || c == '&')
      fprintf (fp, "\\u%04x", c);
    else
      fputc (c, fp);
  }
  fputc ('"', fp);
}

/* the title for HTML text */
static void html_text (FILE * fp, const char * s)
{
  for (; *s; s++)
    switch (*s) {
    case '<': fputs ("&lt;", fp); break;
    case '>': fputs ("&gt;", fp); break;
    case '&': fputs ("&amp;", fp); break;
    case '"': fputs ("&quot;", fp); break;
    default: fputc (*s, fp);
    }
}

/* a JSON number, null if not finite */
static void json_number (FILE * fp, double x)
{
  if (isfinite (x))
    fprintf (fp, "%.17g", x);
  else
    fputs ("null", fp);
}

/* the box of everything, over the finite coordinates */
static void scene_box (const WebglScene * S, double * lo, double * hi)
{
  lo[0] = lo[1] = HUGE_VAL, hi[0] = hi[1] = -HUGE_VAL;
  for (int l = 0; l < S->nlayer; l++) {
    const Layer * L = &S->layer[l];
    for (int b = 0; b < L->n; b++)
      for (int d = 0; d < 2; d++) {
        double a = L->box[4 * b + d], w = L->box[4 * b + 2 + d];
        if (isfinite (a) && isfinite (w))
          lo[d] = fmin (lo[d], a), hi[d] = fmax (hi[d], a + w);
      }
    for (int i = 0; i < L->np; i++)
      for (int d = 0; d < 2; d++)
        if (isfinite (L->xy[2 * i + d]))
          lo[d] = fmin (lo[d], L->xy[2 * i + d]),
          hi[d] = fmax (hi[d], L->xy[2 * i + d]);
  }
  for (int d = 0; d < 2; d++)
    if (lo[d] > hi[d])
      lo[d] = hi[d] = 0.;
}

typedef struct {
  FILE * fp;
  int    nblob;
} Out;

/* one typed array as a <script> of base64, its id returned */
static int blob (Out * o, const void * p, size_t bytes)
{
  int id = o->nblob++;
  fprintf (o->fp, "<script id=\"b%d\" type=\"application/octet-stream\">", id);
  base64 (o->fp, p, bytes);
  fputs ("</script>\n", o->fp);
  return id;
}

/* coordinates as offsets from the scene's corner, in single precision */
static float * offsets (const double * x, size_t n, int stride,
                        const double * lo)
{
  float * f = xmalloc (n * sizeof (float));
  for (size_t i = 0; i < n; i++)
    f[i] = (float) (x[i] - lo[(i % (size_t) stride) % 2]);
  return f;
}

static const char * type_name [] = { "blocks", "lines", "triangles",
                                     "points" };
static const char * where_name [] = { "cell", "block", "point", "element" };

/* the blobs of every layer, then the manifest naming them */
static void write_data (const WebglScene * S, FILE * fp)
{
  Out o = { fp, 0 };
  double lo[2], hi[2];
  scene_box (S, lo, hi);
  int (* id)[3] = xmalloc ((size_t) (S->nlayer ? S->nlayer : 1) *
                           sizeof (*id));
  int nf = 0;
  for (int l = 0; l < S->nlayer; l++)
    nf += S->layer[l].nfield;
  int * fid = xmalloc ((size_t) (nf ? nf : 1) * sizeof (int)), k = 0;
  for (int l = 0; l < S->nlayer; l++) {
    const Layer * L = &S->layer[l];
    id[l][0] = id[l][1] = id[l][2] = -1;
    if (L->type == BLOCKS) {
      /* x0, y0 offsets; widths as they are */
      float * f = xmalloc (4 * (size_t) L->n * sizeof (float));
      for (int b = 0; b < L->n; b++)
        for (int c = 0; c < 4; c++)
          f[4 * b + c] = (float) (L->box[4 * b + c] - (c < 2 ? lo[c] : 0.));
      id[l][0] = blob (&o, f, 4 * (size_t) L->n * sizeof (float));
      free (f);
    }
    else {
      float * f = offsets (L->xy, 2 * (size_t) L->np, 2, lo);
      id[l][1] = blob (&o, f, 2 * (size_t) L->np * sizeof (float));
      free (f);
      if (L->type != POINTS)
        id[l][2] = blob (&o, L->e, (size_t) L->ne *
                         (L->type == LINES ? 2u : 3u) * sizeof (uint32_t));
    }
    for (int f = 0; f < L->nfield; f++)
      fid[k++] = blob (&o, L->field[f].v, L->field[f].n * sizeof (float));
  }

  fputs ("<script id=\"manifest\" type=\"application/json\">{\"title\":", fp);
  json_string (fp, S->title);
  fputs (",\"origin\":[", fp);
  json_number (fp, lo[0]), fputc (',', fp), json_number (fp, lo[1]);
  fputs ("],\"size\":[", fp);
  json_number (fp, hi[0] - lo[0]), fputc (',', fp);
  json_number (fp, hi[1] - lo[1]);
  fputs ("],\"layers\":[", fp);
  k = 0;
  for (int l = 0; l < S->nlayer; l++) {
    const Layer * L = &S->layer[l];
    fprintf (fp, "%s{\"type\":\"%s\",\"name\":", l ? "," : "",
             type_name[L->type]);
    json_string (fp, L->name);
    if (L->type == BLOCKS)
      fprintf (fp, ",\"n\":%d,\"m\":%d,\"box\":\"b%d\"", L->n, L->m,
               id[l][0]);
    else {
      fprintf (fp, ",\"np\":%d,\"ne\":%d,\"xy\":\"b%d\"", L->np, L->ne,
               id[l][1]);
      if (L->type != POINTS)
        fprintf (fp, ",\"e\":\"b%d\"", id[l][2]);
    }
    fputs (",\"fields\":[", fp);
    for (int f = 0; f < L->nfield; f++) {
      const Field * F = &L->field[f];
      double mn = HUGE_VAL, mx = -HUGE_VAL;
      for (size_t i = 0; i < F->n; i++)
        if (isfinite (F->v[i]))
          mn = fmin (mn, F->v[i]), mx = fmax (mx, F->v[i]);
      fprintf (fp, "%s{\"name\":", f ? "," : "");
      json_string (fp, F->name);
      fprintf (fp, ",\"where\":\"%s\",\"categorical\":%d,\"data\":\"b%d\","
               "\"min\":", where_name[F->where], F->categorical, fid[k++]);
      json_number (fp, mn), fputs (",\"max\":", fp), json_number (fp, mx);
      fputc ('}', fp);
    }
    fputs ("]}", fp);
  }
  fputs ("]}</script>\n", fp);
  free (id), free (fid);
}

/* the template up to the next mark, the mark's number returned (-1 at the
   end) */
static int template_to_mark (FILE * fp, const unsigned char * t, size_t n,
                             size_t * at)
{
  static const char * marks [] = { "@TITLE@", "@IMPORTMAP@", "@DATA@",
                                    "@VIEWER@" };
  size_t i = *at;
  while (i < n) {
    if (t[i] == '@')
      for (int k = 0; k < 4; k++) {
        size_t len = strlen (marks[k]);
        if (i + len <= n && !memcmp (t + i, marks[k], len)) {
          fwrite (t + *at, 1, i - *at, fp);
          *at = i + len;
          return k;
        }
      }
    i++;
  }
  fwrite (t + *at, 1, n - *at, fp);
  *at = n;
  return -1;
}

int webgl_offline (void)
{
  return vendor_three_module_min_js_len > 0 && vendor_OrbitControls_js_len > 0;
}

int webgl_scene_write (const WebglScene * S, const char * path)
{
  if (!S || !path)
    return -1;
  FILE * fp = fopen (path, "wb");
  if (!fp)
    return -1;
  size_t at = 0;
  int mark;
  while ((mark = template_to_mark (fp, viewer_viewer_html,
                                   viewer_viewer_html_len, &at)) >= 0)
    switch (mark) {
    case 0:
      html_text (fp, S->title);
      break;
    case 1:
      fputs ("{\"imports\":{\"three\":\"", fp);
      if (webgl_offline ()) {
        fputs ("data:text/javascript;base64,", fp);
        base64 (fp, vendor_three_module_min_js,
                vendor_three_module_min_js_len);
      }
      else
        fputs (WEBGL_THREE_CDN "build/three.module.min.js", fp);
      fputs ("\",\"three/addons/controls/OrbitControls.js\":\"", fp);
      if (webgl_offline ()) {
        fputs ("data:text/javascript;base64,", fp);
        base64 (fp, vendor_OrbitControls_js, vendor_OrbitControls_js_len);
      }
      else
        fputs (WEBGL_THREE_CDN "examples/jsm/controls/OrbitControls.js", fp);
      fputs ("\"}}", fp);
      break;
    case 2:
      write_data (S, fp);
      break;
    case 3:
      fwrite (viewer_viewer_js, 1, viewer_viewer_js_len, fp);
      break;
    }
  return fclose (fp) ? -1 : 0;
}
