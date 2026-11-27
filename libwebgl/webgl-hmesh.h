#ifndef _WEBGL_HMESH_H_
#define _WEBGL_HMESH_H_

  /*
  .. A libhmesh manifold in a scene (webgl.h); 2D for now.
  ..
  .. webgl_hmesh (S, name, M): its points and top simplices as a layer:
  .. points (m = 0), lines (m = 1: the edges) or triangles (m = 2). Points
  .. are renumbered 0 ... n - 1 in the order of their ids, holes left out.
  .. Fields:
  ..
  ..   every HmReal field on the points        WEBGL_POINT
  ..   every HmReal field on the top level     WEBGL_ELEMENT
  ..     (a vector's components as name.x, name.y)
  ..   id: each point's id, each simplex's     the mesh's own numbers, to
  ..                                           find them in the code
  ..   length (m = 1) or area (m = 2)          hm_measure ()
  ..
  .. The layer's number, or -1 (not 2D, m > 2, or an argument NULL).
  */

  #include "hmesh.h"
  #include "webgl.h"

  int webgl_hmesh (WebglScene * S, const char * name, Hmesh * M);

#endif
