#include <cwrap/cache.h>
static int count (cache (int) c) { return c.n; }
int main (void)
{
  cache (int) v;
  cache_push (v, 1);
  return count (v) - 1 + count ((v)) - 1;
}
