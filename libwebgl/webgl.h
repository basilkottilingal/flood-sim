#ifndef _WEBGL_H_
#define _WEBGL_H_

  /*
  .. A scene for the browser: layers of geometry and fields on them,
  .. written as one HTML file that draws them with three.js (WebGL), with
  .. nothing to install: the viewer and the data are inside the file, and
  .. three.js too when the library was built with it (sh vendor/fetch.sh,
  .. then the pages need no network); otherwise the pages load three.js
  .. r160 from jsDelivr. Firefox opens them from disk.
  ..
  .. 2D for now. A scene holds layers:
  ..
  ..   blocks      n blocks of m x m cells, each a box (x0, y0, width,
  ..               height): a block-structured grid, as libb-amr's leaves
  ..   lines       points and segments between them: a curve, a mesh's
  ..               edges
  ..   triangles   points and triangles
  ..   points      points alone
  ..
  .. and fields on a layer, by where they live:
  ..
  ..   WEBGL_CELL     a value a cell of a blocks layer, n m^2, block by
  ..                  block, i fastest within a block
  ..   WEBGL_BLOCK    a value a block (a level, an owner)
  ..   WEBGL_POINT    a value a point
  ..   WEBGL_ELEMENT  a value a segment or triangle
  ..
  .. A categorical field (levels, ranks) is drawn with distinct colours,
  .. the others with a colour map over their range. NaN is drawn grey.
  ..
  .. The scene copies what it is given: the caller's arrays may go after
  .. each call. Positions are stored as offsets from the scene's first
  .. point, in single precision, and drawn so; picking reports them back
  .. with the offset added in double. Field values are single precision.
  ..
  ..   webgl_scene_new (title)            an empty scene
  ..   webgl_scene_free (S)
  ..   webgl_blocks (S, name, n, m, box)  a blocks layer: box[4 n] the
  ..                                      blocks' x0, y0, width, height
  ..   webgl_lines (S, name, np, xy, ne, e)       xy[2 np], e[2 ne]
  ..   webgl_triangles (S, name, np, xy, ne, t)   t[3 ne]
  ..   webgl_points (S, name, np, xy)
  ..   webgl_field (S, layer, name, where, v, categorical)
  ..                                      a field on a layer, v as long as
  ..                                      where says
  ..   webgl_scene_write (S, path)        the HTML file
  ..   webgl_offline ()                   1 if three.js is in the pages, 0
  ..                                      if they load it from jsDelivr
  ..
  .. Each layer call returns the layer's number, each other call 0, and
  .. every call -1 for arguments it cannot use (nothing changed): an index
  .. out of range, a field of the wrong kind for its layer, a NULL.
  */

  #include <stdint.h>

  typedef struct WebglScene WebglScene;

  enum { WEBGL_CELL, WEBGL_BLOCK, WEBGL_POINT, WEBGL_ELEMENT };

  WebglScene * webgl_scene_new (const char * title);
  void         webgl_scene_free (WebglScene * S);
  int          webgl_blocks (WebglScene * S, const char * name, int n, int m,
                             const double * box);
  int          webgl_lines (WebglScene * S, const char * name, int np,
                            const double * xy, int ne, const uint32_t * e);
  int          webgl_triangles (WebglScene * S, const char * name, int np,
                                const double * xy, int ne,
                                const uint32_t * t);
  int          webgl_points (WebglScene * S, const char * name, int np,
                             const double * xy);
  int          webgl_field (WebglScene * S, int layer, const char * name,
                            int where, const float * v, int categorical);
  int          webgl_scene_write (const WebglScene * S, const char * path);
  int          webgl_offline (void);

  /* where the pages load three.js from when it is not in them */
  #define WEBGL_THREE_CDN "https://cdn.jsdelivr.net/npm/three@0.160.0/"

#endif
