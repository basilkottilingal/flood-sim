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
        {1, 0}, {1, 1}, {0, 1}, {-1, 1},
        {-1, 0}, {-1, -1}, {0, -1}, {1, -1}
      };

  typedef enum
  {
    FN_GENERAL               = 0,
    FN_BOUNDARY              = 32,
    FN_OCEAN                 = 128,
    FN_FLAT                  = 64,
    //FN_RESERVOIR           = 1,
    //FN_RESERVOIR_CATCHMENT = 2,
    //FN_RESERVOIR_BREACH    = 4,
  } FLOW_NODE;

  typedef enum
  {
    FLOW_D8 = 1,
    FLOW_DINFTY,
    FLOW_MFD,
  } FLOW_DRAIN;

  typedef struct
  {
    uint8_t ** dir;
    uint8_t ** invDir;
    double  ** accumulation;
    FLOW_DRAIN type;
    int w, h;
  } FlowNetwork;

  int  pop_bit (uint8_t * v);
  int  flow_network           ( DEM dem, FLOW_DRAIN type, FlowNetwork * network );
  void flow_network_free      ( FlowNetwork );
  int  flow_accumulation      ( FlowNetwork *, DEM, MinPQ * pq );
  int  flow_remove_pits       ( DEM * );
  void flow_network_error     ( int type );
  int  flow_routine           ( DEM );
  
#endif
