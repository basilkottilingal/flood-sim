#include <stdio.h>
static int calls;
static void foreach (void) { calls++; }
int main (void)
{
  foreach ();
  foreach ( );
  printf ("%d\n", calls);
  return 0;
}
