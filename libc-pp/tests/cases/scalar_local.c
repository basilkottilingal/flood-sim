#include <stdio.h>
scalar g;
static int sum (int k)
{
  scalar a;
  int total = 0;
  for (int i = 0; i < 3; i++) {
    scalar b, c;
    total += b + c;
    if (i == k)
      return total + a;
  }
  return total;
}
static scalar counted (void)
{
  static scalar s;
  return s;
}
int main (void)
{
  int unused;
  printf ("%d\n", sum (1));
  scalar kept = scalar_new ();
  printf ("in use %d kept %d gs %d globals %d\n", scalar_count (), kept,
          g + counted (), n_global_scalars);
  return 0;
}
