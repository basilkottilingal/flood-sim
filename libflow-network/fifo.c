#include "fifo.h"

#include <stdlib.h>

#include <stdint.h>
#include <stdlib.h>
#include <stdbool.h>
#include <stddef.h>

int queue_init (Queue *q)
{
  size_t page = 1<<12;
  q->buf = malloc (page * sizeof(QPair));
  if (!q->buf) return -1;
  q->head  = 0;
  q->count = 0;
  q->cap   = page;
  return 0;
}

void queue_destroy (Queue *q)
{
  free (q->buf);
  q->buf = NULL;
  q->head = q->count = q->cap = 0;
}

static inline size_t queue_size  (const Queue *q)  { return q->count; }
static inline bool   queue_empty (const Queue *q)  { return q->count == 0; }

static int queue_grow (Queue *q)
{
  size_t new_cap = q->cap << 1;
  QPair *nb = malloc (new_cap * sizeof(QPair));
  if (!nb) return -1;

  /* Unwrap into the new buffer so head becomes 0. */
  size_t mask = q->cap - 1;
  for (size_t i = 0; i < q->count; i++)
    nb[i] = q->buf [(q->head + i) & mask];

  free(q->buf);
  q->buf  = nb;
  q->head = 0;
  q->cap  = new_cap;
  return 0;
}

int queue_push (Queue *q, int x, int y)
{
  if (q->count == q->cap && queue_grow(q) != 0) return -1;
  size_t tail = (q->head + q->count) & (q->cap - 1);
  q->buf[tail].x = x;
  q->buf[tail].y = y;
  q->count++;
  return 0;
}

int queue_pop (Queue *q, int * x, int * y)
{
  if (q->count == 0) return 0;
  *x = q->buf[q->head].x;
  *y = q->buf[q->head].y;
  q->head = (q->head + 1) & (q->cap - 1);
  q->count--;
  return 1;
}

#if 0
static inline void queue_clear(Queue *q)
{
  q->head = 0;
  q->count = 0;
}
static const QPair *queue_peek (const Queue *q)
{
  return q->count ? &q->buf[q->head] : NULL;
}

/* Iterate front to back without popping. */
#define QUEUE_FOREACH(q, idx, p) \
  for (size_t idx = 0; \
     idx < (q)->count && ((p) = &(q)->buf[((q)->head + idx) & ((q)->cap - 1)], 1); \
     idx++)
#endif
