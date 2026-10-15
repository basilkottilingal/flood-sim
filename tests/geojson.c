#include "geojson.h"
#include "ast.h"

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
  geojson (argv [1], &g, argv [2]); 
  ast_free_all ();
}
