#include <stdio.h>
static int cache = 2;
static int foreach (int x) { return x + cache; }
int main (void)
{
  int in = 1, V = 3;
  printf ("%d\n", foreach (in) * foreach (V));
  return 0;
}
