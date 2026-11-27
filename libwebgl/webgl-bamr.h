#ifndef _WEBGL_BAMR_H_
#define _WEBGL_BAMR_H_

  /*
  .. A .bamr file (libb-amr's tree_output (), the format in
  .. libb-amr/tree-output.h) in a scene (webgl.h): plain C, no libb-amr, so
  .. a saved grid is drawn after the run, from any number of processes'
  .. writing (the file is the same whatever their number). 2D for now.
  ..
  .. webgl_bamr (S, path): the file's leaves as a blocks layer named after
  .. the file (its last path component), a field a scalar (WEBGL_CELL,
  .. named as in the file, single precision whatever the file's), and each
  .. leaf's level (WEBGL_BLOCK, categorical). The file is checked first:
  .. the magic, version 1, the byte order, codec 0, its size against its
  .. sections, the levels (0 to 30) and the indices in range. The layer's
  .. number, or -1 with nothing added (no such file, not 2D, a check
  .. failed). The file holds no owner ranks: it does not depend on their
  .. number.
  ..
  .. tools/bamr2html.c makes the command: bamr2html FILE.bamr [OUT.html].
  */

  #include "webgl.h"

  int webgl_bamr (WebglScene * S, const char * path);

#endif
