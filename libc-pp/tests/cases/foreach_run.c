#include <stdio.h>
struct ints { int n; int *a; };
int main (void)
{
  int data[] = { 1, 2, 3, 4 };
  struct ints V = { 4, data };
  int sum = 0;
  foreach (v in V)
    sum += v;
  foreach (x
           in V) {
    foreach (y in V)
      sum += x * y;
    foreach (x in V)
      sum += x;
  }
  int unused;
  printf ("%d\n", sum);
  return 0;
}
