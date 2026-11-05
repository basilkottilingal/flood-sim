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
  foreach_scalar (s)
    printf ("%d,", s);
  printf ("\n");
  foreach_scalar (s)
    if (s == 1)
      scalar_free (s);
  foreach_scalar (s)
    printf ("%d,", s);
  printf ("\n");
  for (int i = scalar_count (); i < NSCALAR_MAX; i++)
    (void)scalar_new ();
  int n = 0, last = -1;
  foreach_scalar (s) {
    n++;
    last = s;
  }
  printf ("%d %d %d %d\n", scalar_count (), scalar_in_use (NSCALAR_MAX - 1), n, last);
  if (argc > 1)
    (void)scalar_new ();
  (void)argv;
  return 0;
}
