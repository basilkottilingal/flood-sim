#include <stdio.h>
#include <cwrap/cache.h>
typedef struct { double x, y; } point;
cache (int) squares = { .update = NULL };
static int rebuilt;
static void rebuild (cache (int) * c)
{
  rebuilt++;
  cache_clear (*c);
  for (int i = 0; i < 20; i++)
    cache_push (*c, i * i);
}
static int sum (cache (int) * c)
{
  int s = 0;
  for (int i = 0; i < c->n; i++) s += c->a[i];
  return s;
}
int main (void)
{
  cache (point) pts;
  cache (char *) words = { .update = NULL };
  int cache = 3;                        /* still a name */
  squares.update = rebuild;
  squares.is_dirty = 1;
  int total = 0;
  foreach (s in squares)
    total += s;
  foreach (s in squares)
    total += s;
  cache_push (pts, (point) { 1, 2 });
  cache_push (words, "hi");
  total += sum (&squares);
  printf ("%d %d %d %g %s %d\n", total, rebuilt, squares.n, pts.a[0].y, words.a[0], cache);
  cache_free (squares);
  return 0;
}
