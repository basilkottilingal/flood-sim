#include <cwrap/cache.h>
typedef int *int_p;
cache (int) g;
typedef cache (int) ints;
cache (int) f (void)
{
  cache (int) v, w;
  w = v;
  cache (int) u = v;
  struct { cache (int) c; } s;
  s.c = v;
  return v;
}
cache (int *) p;
cache (int_p) q;
