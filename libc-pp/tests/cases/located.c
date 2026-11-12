#include <stdio.h>

face_vector F;
vertex_scalar v;
scalar c;

static void show (const char * what, scalar s)
{
  printf ("%s %s %d; ", what, scalar_name (s), scalar_location (s));
}

static int counted (void)
{
  static face_vector kept;
  return scalar_location (kept.x);
}

int main (void)
{
  show ("F.x", F.x);
#if CWRAP_DIMENSION > 1
  show ("F.y", F.y);
#endif
#if CWRAP_DIMENSION > 2
  show ("F.z", F.z);
#endif
  show ("v", v);
  show ("c", c);
  int before = scalar_count (), inner;
  {
    face_vector G;
    vertex_scalar w;
    show ("G.x", G.x);
    show ("w", w);
    inner = scalar_count () - before;
  }
  printf ("inner %d after %d static %d\n", inner, scalar_count () - before, counted ());
  return 0;
}
