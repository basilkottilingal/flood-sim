#include "loop.h"
struct S { int v; };
int main (void)
{
  struct S s = { 1 };
  range (0, 4,
         1) {
    s.nope = i;
  }
  int after = s.nope2;
  each (&s) {
  }
  return after;
}
