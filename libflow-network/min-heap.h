#ifndef _FLOW_NETWORK_MIN_PRIORITY_QUEUE_H
#define _FLOW_NETWORK_MIN_PRIORITY_QUEUE_H
  typedef struct
  {
    float key;
    uint16_t i, j;
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
  int  pq_push   (MinPQ *pq, float key, uint16_t i, uint16_t j);
  int  pq_pop    (MinPQ *pq, float * key, uint16_t * i, uint16_t *j);
#endif
