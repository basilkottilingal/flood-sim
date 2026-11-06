/* Block macros for the @loop tests */
#include <cwrap/cache.h>

@loop repeat (int n) : k
{
  for (int k = 0; k < n; k++)
    {...}
}

@loop range (int lo, int hi, int step) : i
{
  int count = 0;
  for (int i = lo; i < hi; i += step) {
    count++;
    {...}
  }
}

typedef struct { int a [8]; int n; } Bag;

/* b is a parameter the block sees, and the argument is often named b too */
@loop each (const Bag * b) : x, b
{
  for (int idx = 0; idx < b->n; idx++) {
    int x = b->a [idx];
    {...}
  }
}

typedef struct { int depth, level; } P;

/* a typedef of its own, and locals spelled like members and designators */
@loop walk (P p) : level, q
{
  typedef struct { int depth; } Hidden;
  Hidden h = { p.depth };
  int depth = h.depth;
  for (int level = 0; level < depth; level++) {
    P q = { .depth = depth, .level = level };
    {...}
  }
}

/* uses another loop, and passes on one of the names that loop gives it */
@loop walk2 (int d) : level, sum
{
  int sum = 0;
  walk ((P) { d, 0 }) {
    sum += level;
    (void) q;
    {...}
  }
}

/* a cache and a foreach of its own */
@loop squares (int n) : v
{
  cache (int) c;
  for (int i = 0; i < n; i++)
    cache_push (c, i * i);
  foreach (v in c)
    {...}
}

/* the same definition twice is fine, as from a header included twice */
@loop repeat (int n) : k
{
  for (int k = 0; k < n; k++)
    {...}
}
