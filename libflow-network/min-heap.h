#ifndef _FLOW_NETWORK_MIN_PRIORITY_QUEUE_H
#define _FLOW_NETWORK_MIN_PRIORITY_QUEUE_H

  #include <stdint.h>

  typedef struct
  {
    double key;
    int x, y;
  } PQNode;

  typedef struct
  {
    PQNode *data;
    int size;
    int capacity;
  } MinPQ;

  /* APIs */
  int  pq_create (MinPQ * pq);
  void pq_free   (MinPQ * pq);
  void pq_truncate (MinPQ * pq);
  int  pq_push   (MinPQ * pq, double key, int x, int y);
  int  pq_pop    (MinPQ * pq, double * key, int * x, int * y);
#endif
