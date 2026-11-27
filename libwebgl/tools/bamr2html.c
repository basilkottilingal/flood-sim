/*
.. bamr2html FILE.bamr [OUT.html]: a saved grid (libb-amr's tree_output ())
.. as one HTML file drawn by three.js (libwebgl): its leaves, a field a
.. scalar, the leaves' levels. OUT defaults to FILE with .html for .bamr.
.. three.js is in the page if libwebgl was built with it (vendor/fetch.sh),
.. else the page loads it from jsDelivr. 2D for now.
..
..   gcc -O2 -I.. bamr2html.c ../webgl.c ../webgl-bamr.c ../webgl-assets.c \
..       -o bamr2html -lm                   (or: make bamr2html, libwebgl)
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "webgl-bamr.h"

int main (int argc, char ** argv)
{
  if (argc < 2 || argc > 3 || !strcmp (argv[1], "-h") ||
      !strcmp (argv[1], "--help")) {
    fprintf (stderr, "usage: bamr2html FILE.bamr [OUT.html]\n");
    return 2;
  }
  const char * in = argv[1];
  char * out;
  if (argc == 3)
    out = strdup (argv[2]);
  else {
    size_t n = strlen (in);
    out = malloc (n + 6);
    if (out) {
      strcpy (out, in);
      if (n > 5 && !strcmp (in + n - 5, ".bamr"))
        out[n - 5] = '\0';
      strcat (out, ".html");
    }
  }
  if (!out)
    return 1;
  const char * base = strrchr (in, '/');
  WebglScene * S = webgl_scene_new (base ? base + 1 : in);
  int layer = webgl_bamr (S, in);
  if (layer < 0) {
    fprintf (stderr, "bamr2html: %s: not a 2D .bamr file, or not "
             "readable\n", in);
    webgl_scene_free (S), free (out);
    return 1;
  }
  if (webgl_scene_write (S, out)) {
    fprintf (stderr, "bamr2html: %s: not written\n", out);
    webgl_scene_free (S), free (out);
    return 1;
  }
  printf ("%s (three.js %s)\n", out,
          webgl_offline () ? "inside" : "from jsDelivr");
  webgl_scene_free (S), free (out);
  return 0;
}
