#include <stdio.h>
#include <cwrap/cache.h>
typedef struct { double x, y; } point;
int main (void)
{
  cache (point) p1;
  cache (point) p3;
  cache (struct { double x, y; }) p2;
  do {
    point p = (point) {1.0, 1.0};
    /* this works */
    cache_push (p1, p);
    /* this won't, since cache_push () is a macro*/
    cache_push (p3, (point) {1.0, 1.0});
  } while (0);
  foreach (p  in p1)
  {
    printf ("\n{%g %g}", p.x, p.y);
  }
  foreach (p  in p3)
  {
    printf ("\n{%g %g}", p.x, p.y);
  }
}
