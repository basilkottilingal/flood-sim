#ifndef _FLOW_NETWORK_FLOW_H_
#define _FLOW_NETWORK_FLOW_H_

  #include "geotiff.h"
  #include "min-heap.h"

  typedef enum
  {
    FLOW_D8 = 1,
    FLOW_DINFTY,
    FLOW_MFD,
  } FLOW_DRAIN;

  typedef struct
  {
    const double * const * raster; /* original elevation raster          */
    double  ** elevation;          /* pit filled elevation               */
    uint8_t ** dir;                /* flow direction. only D8 iplemented */
    uint8_t ** invDir;
    double  ** accumulation;
    FLOW_DRAIN type;
    int w, h;
  } FlowNetwork;

  /* API */
  int  flow_network       ( DEM dem, FLOW_DRAIN type, FlowNetwork * network );
  void flow_network_free  ( FlowNetwork * );
  void flow_network_error ( int err );
  
#endif
