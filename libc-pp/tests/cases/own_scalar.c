#include <stdio.h>
typedef double scalar;
typedef struct { double x, y; } vector;
scalar a = 1.5;
vector v;
int main (void)
{
  scalar b;
  b = a;
  v.x = b;
  printf ("%g\n", v.x);
  return 0;
}
