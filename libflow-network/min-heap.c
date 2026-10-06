#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <stdint.h>
#include <float.h>

#include "min-heap.h"

int pq_create (MinPQ * pq)
{
  size_t page = (size_t) 1 << 14;
  pq->data = malloc (page);
  if (pq->data == NULL)
    return -1;
  pq->size = 0;
  pq->capacity = page / sizeof (PQNode);
  return 0;
}

void pq_truncate (MinPQ * pq)
{
  int n = pq->size ? pq->size : 1;
  PQNode * data = realloc (pq->data, n * sizeof (PQNode));
  if (data != NULL)
  {
    pq->data = data;
    pq->capacity = n;
  }
}

static inline void swap (PQNode *a, PQNode *b)
{
  PQNode tmp = *a;
  *a = *b;
  *b = tmp;
}

static void sift_up (MinPQ *pq, int i)
{
  while (i > 0)
  {
    int parent = (i - 1) / 2;
    if (pq->data[parent].key <= pq->data[i].key) break;
    swap (&pq->data[parent], &pq->data[i]);
    i = parent;
  }
}

static void sift_down(MinPQ *pq, int i)
{
  for (;;)
  {
    int left  = 2 * i + 1;
    int right = 2 * i + 2;
    int smallest = i;

    if (left  < pq->size && pq->data[left].key  < pq->data[smallest].key)
      smallest = left;
    if (right < pq->size && pq->data[right].key < pq->data[smallest].key)
      smallest = right;

    if (smallest == i)
      return;
    swap (&pq->data[i], &pq->data[smallest]);
    i = smallest;
  }
}

int pq_push (MinPQ *pq, double key, int x, int y)
{
  if (pq->size == pq->capacity)
  {
    pq->capacity *= 2;
    pq->data = realloc (pq->data, sizeof(PQNode) * pq->capacity);
    if ( pq->data == NULL )
      return -1;
  }
  pq->data[pq->size].key = key;
  pq->data[pq->size].x = x;
  pq->data[pq->size].y = y;
  sift_up (pq, pq->size);
  pq->size++;
  return 0;
}

double pq_min_key (MinPQ *pq)
{
  if (pq->size == 0)
    return FLT_MAX; /* handle error */
  return pq->data[0].key;
}

int pq_pop (MinPQ *pq, double * key, int * x, int * y)
{
  if (pq->size == 0)
    return 0;

  PQNode top = pq->data[0];
  *key = top.key, *x = top.x, *y = top.y;
  pq->size--;
  pq->data[0] = pq->data[pq->size];
  sift_down (pq, 0);
  
  return 1;
}

void pq_free(MinPQ *pq)
{
  free (pq->data);
}
