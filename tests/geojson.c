#include "libgeojson/geojson.h"

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

/*
.. $ make geojson.s && ./geojson.s ../france-regions.geojson nom
*/

int main (int argc, const char ** argv)
{
  assert (argc >= 3);
  GeoJSON g;
  if (geojson (argv [1], &g, argv [2]))
  {
    printf ("cannot find parse geojson");
    exit (-1);
  }
  ast_free (g.tree);
}
