#ifndef _FLOW_NETWORK_MIN_PRIORITY_QUEUE_H
#define _FLOW_NETWORK_MIN_PRIORITY_QUEUE_H

  #include <stdint.h>

  typedef struct
  {
    float key;
    int x, y;
  } PQNode;

  typedef struct
  {
    PQNode *data;
    int size;
    int capacity;
  } MinPQ;

  /* APIs */
  int  pq_create (MinPQ *pq);
  void pq_free   (MinPQ *pq);
  int  pq_push   (MinPQ *pq, float key, int x, int y);
  int  pq_pop    (MinPQ *pq, float * key, int * x, int * y);
#endif
