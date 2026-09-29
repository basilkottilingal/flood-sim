#include "graph-analysis.h"
#include "geotiff.h"
#include "filemap.h"
#include "flow.h"
#include "min-heap.h"
#include "pgm.h"
#include "grid.h"

#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <stdlib.h>
#include <errno.h>
#include <assert.h>
#include <math.h>
#include <float.h>

#define error(e) filemap_close_all(e)

int flow_graph (DEM dem, FlowGraph * graph)
{
  FlowNetwork network;
  if (flow_network (dem, FLOW_D8, &network))
  {
    error ("flow_graph () : flow_network ();");
    return -1;
  }

  MinPQ sink;
  if (flow_accumulation (&network, dem, &sink))
  {
    error ("flow_graph () : flow_accumulation ()");
    return -1;
  }

  int w = network.w, h = network.h, hw = w * h;
  double ** elevation = dem.raster;
  uint8_t ** dir         = network.dir;
  uint8_t ** inverse     = grid (uint8_t, w, h);
  uint8_t ** nans        = grid (uint8_t, w, h);
  uint8_t ** visited     = grid (uint8_t, w, h);
  uint8_t ** pit_centers = grid (uint8_t, w, h);
  uint8_t ** tree_tag    = grid (uint8_t, w, h);

  /* inverse is basically inversing the "dir" evaluated by D8 */
  for (int y=0; y<h; ++y)
    for (int x=0; x<w; ++x)
    {
if(isnan (elevation [y][x]))
nans [y][x] = 255u;
      uint8_t code = dir [y][x];
      if (code == 0)
        continue;
      int nbr = pop_bit (&code), xnbr = flow_neighbor [nbr].x, ynbr = flow_neighbor [nbr].y;
      inverse [y+ynbr][x+xnbr] |= (uint8_t) 1 << ( (nbr+4) % 8 );
    }

  unsigned char tag = 2u;
  double elev; int x, y;
  while (pq_pop (&sink, &elev, &x, &y))
  {
    tag = ((int) tag % 245u) + 5u;
    if (elevation [y][x] <= 0.f)
      tag = 0u ;               /* sink node (maybe) sea */
    else if (x >= w || x < 0 || y >= h || y < 0)
      tag = 255u ;             /* sink node outside box */

    pit_centers [y][x] = 128;

    /* tree traverseal */
    int depth = 0;
    do
    {
      while ( inverse [y][x] )
      {
        depth++;
        int nbr = pop_bit ( &inverse [y][x] );
        x += flow_neighbor [nbr].x;
        y += flow_neighbor [nbr].y;
        assert (depth < hw);
      }

      tree_tag [y][x] = tag;
      visited  [y][x] = 255u;

      if (!depth--)
        break; 

      uint8_t code = dir [y][x]; assert (code);
      int nbr = pop_bit ( &code );
      x += flow_neighbor [nbr].x;
      y += flow_neighbor [nbr].y;
    } while (1);

  }
  
#if 0
  for (int y=0; y<h; ++y)
    for (int x=0; x<w; ++x)
      if (!visited [y][x])
        printf ("{%u}", dir [y][x]);
#endif

//test
flow_network_grayscale ((FlowNetwork) {.dir = pit_centers, .w = w, .h = h}, "pits.pgm", w, h);
flow_network_grayscale ((FlowNetwork) {.dir = inverse, .w = w, .h = h}, "dag-inverse.pgm", w, h);
flow_network_grayscale ((FlowNetwork) {.dir = tree_tag, .w = w, .h = h}, "pit-tag.pgm", w, h);
flow_network_grayscale ((FlowNetwork) {.dir = visited, .w = w, .h = h}, "visited.pgm", w, h);
flow_network_grayscale ((FlowNetwork) {.dir = nans, .w = w, .h = h}, "NAN.pgm", w, h);

  /* travering through each DAG rooted about a pit center */
  
//fixme : graph not set. not all grids freed yet.
  flow_network_free (network);
grid_free (nans);
  grid_free (inverse);
  grid_free (pit_centers);
  grid_free (tree_tag);
  pq_free (&sink);
  return 0;
}

#undef error
