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
    error ("flow_graph () : flow_network (dem, FLOW_D8, &network);");
    return -1;
  }

  return 0;
}

#undef error
