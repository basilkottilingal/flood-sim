#ifndef _FLOW_NETWORK_GRAPH_H
#define _FLOW_NETWORK_GRAPH_H

  #include "min-heap.h"
  #include "geotiff.h" 
  #include <stdint.h>
  /*
  .. Flow Graph is a forest (collection of trees / DAGs i.e directed acyclic
  .. graphs). The direction is the opposite of the flow accumuation direction
  .. determined by D8 algorithm. So you can consider pits determined by D8 as
  .. the root nodes of each DAGs in Flow Graph. There are other cases where
  .. a D8 accumulation algorithm will drain to a boundary node or nodes
  .. classified as ocean/sea.
  */ 
  typedef struct
  {
    uint8_t ** dir;
    MinPQ pits;
    MinPQ boundary;
    MinPQ ocean;
    int w, h;
  } FlowGraph;

  /* APIs */
  int  flow_graph ( DEM dem, FlowGraph * g );
  //free
  
#endif
