/*
.. A manifold as a layer: see webgl-hmesh.h.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "webgl-hmesh.h"

static void * xmalloc (size_t n)
{
  void * p = malloc (n ? n : 1);
  if (!p) {
    fprintf (stderr, "webgl: out of memory\n");
    abort ();
  }
  return p;
}

/* each HmReal field on level k as a field of the layer, where */
static void fields_of (WebglScene * S, int layer, Hmesh * M, int k,
                       int where, const uint32_t * ids, uint32_t n)
{
  float * v = xmalloc ((size_t) n * sizeof (float));
  for (int f = 0; f < M->nfield; f++) {
    const HmField * F = &M->field[f];
    if (F->level != k || !F->real)
      continue;
    for (int c = 0; c < F->ncomp; c++) {
      const HmReal * a = hm_val (M, f, c);
      for (uint32_t i = 0; i < n; i++)
        v[i] = (float) a[ids[i]];
      char name[128];
      if (F->ncomp > 1)
        snprintf (name, sizeof (name), "%s.%c", F->name ? F->name : "f",
                  "xyz"[c]);
      else
        snprintf (name, sizeof (name), "%s", F->name ? F->name : "f");
      webgl_field (S, layer, name, where, v, 0);
    }
  }
  free (v);
}

int webgl_hmesh (WebglScene * S, const char * name, Hmesh * M)
{
#if HM_DIMENSION != 2
  (void) S, (void) name, (void) M;
  return -1;
#else
  if (!S || !name || !M || M->dim > 2)
    return -1;
  int m = M->dim;
  uint32_t np = hm_count (M, 0), end = hm_end (M, 0), a = 0;
  uint32_t * local = xmalloc ((size_t) end * sizeof (uint32_t));
  uint32_t * pid = xmalloc ((size_t) np * sizeof (uint32_t));
  double * xy = xmalloc (2 * (size_t) np * sizeof (double));
  hm_foreach (M, 0, v) {
    local[v] = a;
    pid[a] = v;
    xy[2 * a] = hm_x (M, 0)[v], xy[2 * a + 1] = hm_x (M, 1)[v];
    a++;
  }
  uint32_t ne = m ? hm_count (M, m) : 0, b = 0;
  uint32_t * e = xmalloc ((size_t) ne * (size_t) (m + 1) * sizeof (uint32_t));
  uint32_t * sid = xmalloc ((size_t) ne * sizeof (uint32_t));
  float * measure = xmalloc ((size_t) ne * sizeof (float));
  if (m)
    hm_foreach (M, m, s) {
      for (int i = 0; i <= m; i++)
        e[(size_t) b * (size_t) (m + 1) + (size_t) i] =
          local[hm_vertex (M, m, s, i)];
      sid[b] = s;
      measure[b] = (float) hm_measure (M, m, s);
      b++;
    }
  int layer = m == 0 ? webgl_points (S, name, (int) np, xy)
    : m == 1 ? webgl_lines (S, name, (int) np, xy, (int) ne, e)
    : webgl_triangles (S, name, (int) np, xy, (int) ne, e);
  if (layer >= 0) {
    float * idf = xmalloc ((size_t) (np > ne ? np : ne) * sizeof (float));
    for (uint32_t i = 0; i < np; i++)
      idf[i] = (float) pid[i];
    webgl_field (S, layer, "id", WEBGL_POINT, idf, 0);
    fields_of (S, layer, M, 0, WEBGL_POINT, pid, np);
    if (m) {
      for (uint32_t i = 0; i < ne; i++)
        idf[i] = (float) sid[i];
      webgl_field (S, layer, "id", WEBGL_ELEMENT, idf, 0);
      webgl_field (S, layer, m == 1 ? "length" : "area", WEBGL_ELEMENT,
                   measure, 0);
      fields_of (S, layer, M, m, WEBGL_ELEMENT, sid, ne);
    }
    free (idf);
  }
  free (local), free (pid), free (xy), free (e), free (sid), free (measure);
  return layer;
#endif
}
