#include "stencil.h"
scalar s, t;
static void fill (scalar a) { foreach () a[] = 1; }
static void twice (scalar a) { fill (a); fill (a); }
static double one (void) { return 1; }
void f (void)
{
  foreach () { s[] = one (); fill (t); }
  foreach () { twice (t); }
}
