#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#ifndef VALUE
#define VALUE 0
#endif
struct point { double x, y; };
int main(void)
{
  struct point p = { .x = 3.0, .y = 4.0 };
  printf("%.1f %d %s\n", sqrt(p.x * p.x + p.y * p.y), VALUE, "a\"b");
  return 0;
}
