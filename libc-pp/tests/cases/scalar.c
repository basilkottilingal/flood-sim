#include "scalar.h"
scalar a;
scalar b, *p, c = 7, d __attribute__ ((unused)), f (void), arr[2];
extern scalar e;
__extension__ static scalar s;
typedef scalar t;
int g (void) { scalar local = 0; return local; }
scalar last;
int main (void) { return a + b + d + s + last + g (); }
