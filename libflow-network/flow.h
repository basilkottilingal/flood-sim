#ifndef _FLOW_NETWORK_FLOW_H
#define _FLOW_NETWORK_FLOW_H

  #include "geotiff.h"
  #include "min-heap.h"

  /*
  .. Neighbors in this order
  ..  3 2 1
  ..  4 . 0
  ..  5 6 7               
  */
  static const struct { int x, y; }
    flow_neighbor [] = 
      { 
        {1, 0},  {1, 1},   {0, 1},  {-1, 1},
        {-1, 0}, {-1, -1}, {0, -1}, {1, -1}
      };

  typedef enum
  {
    FLOW_D8 = 1,
    FLOW_DINFTY,
    FLOW_MFD,
  } FLOW_DRAIN;

  typedef struct
  {
    uint8_t ** dir;
    double  ** accumulation;
    FLOW_DRAIN type;
    int w, h;
  } FlowNetwork;

  int  pop_bit (uint8_t * v);
  int  flow_network           ( DEM dem, FLOW_DRAIN type, FlowNetwork * network );
  void flow_network_free      ( FlowNetwork );
  int  flow_network_grayscale ( FlowNetwork, const char * out, int w, int h );
  int  flow_accumulation_grayscale ( FlowNetwork, const char * out, int w, int h );
  int  flow_accumulation      ( FlowNetwork *, DEM, MinPQ * pq );
  int  flow_remove_pits       ( DEM * );
  void flow_network_error     ( int type );
  
#endif
