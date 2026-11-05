#include <stdio.h>
int main (int argc, char **argv)
{
  scalar a = scalar_new (), b = scalar_new ();
  scalar_free (a);
  scalar c = scalar_new ();
  scalar_free (b);
  scalar_free (b);
  printf ("%d %d %d %d %d\n", a, b, c, scalar_count (), scalar_in_use (b));
  (void)scalar_new ();
  scalar_free (scalar_new ());
  (void)scalar_new ();
  for (int i = scalar_count (); i < NSCALAR_MAX; i++)
    (void)scalar_new ();
  printf ("%d %d\n", scalar_count (), scalar_in_use (NSCALAR_MAX - 1));
  if (argc > 1)
    (void)scalar_new ();
  (void)argv;
  return 0;
}
