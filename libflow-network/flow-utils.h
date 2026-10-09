#ifndef _FLOW_NETWORK_FLOW_UTILS_H_
#define _FLOW_NETWORK_FLOW_UTILS_H_

  #include <stdint.h>

  /*
  .. Neighbors in this order
  ..  3 2 1
  ..  4 . 0
  ..  5 6 7               
  */
  static const struct { int x, y; }
    flow_neighbor [] = 
      { 
        { 1, 0}, { 1,  1}, {0,  1}, {-1,  1},
        {-1, 0}, {-1, -1}, {0, -1}, { 1, -1}
      };

  /*
  .. APIs
  .. 1. pop the last bit of a uin8_t number and return the bit position [0:7].
  ..    Return -1 if v contains 0
  .. 2. validity of a forest of direct acyclic graph (or directed trees)
  .. 3. see if all the water drains for a pit filled elevation
  .. 4. DAG forest corresponding to D8 flow network deduced from elevation.  
  .. 5. forest of directed graphs (not necessarily DAG)
  ..    corresponding to DInfty flow network deduced from elevation.  
  */
  int  pop_bit        (uint8_t * v);
  void DAG_validity   (uint8_t ** dir, uint8_t ** invDir);
  int  drain_validity (double ** elevation, int ** tag, int ntags);
  void D8             (double ** elevation, uint8_t ** dir, uint8_t ** invDir);
  void DInfty         (double ** elevation, uint8_t ** dir, uint8_t ** invDir);
  
#endif
