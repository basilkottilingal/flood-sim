#ifndef _FLOW_NETWORK_FIFO_QUEUE_H_
#define _FLOW_NETWORK_FIFO_QUEUE_H_

  #include <stdlib.h>

  typedef struct
  {
    int x;
    int y;
  } QPair;

  typedef struct
  {
    QPair  *buf;
    size_t head;
    size_t count;
    size_t cap;
  } Queue;
 
  /* API */
  int queue_init (Queue *q);
  void queue_destroy (Queue *q);
  int queue_push (Queue *q, int x, int y);
  int queue_pop (Queue *q, int * x, int * y);
#endif
