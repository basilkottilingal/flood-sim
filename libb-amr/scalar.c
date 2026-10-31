#include "tree.h"
#include "scalar.h"
#include "morton-curve.h"

#include <assert.h>
#include <stdlib.h>

static scalar available [NSCALAR_MAX] = {-1};
static int navailable = 0;
static int n = 0;

scalar scalar_new ()
{
  scalar s = navailable ? available [--navailable] : 
    n < NSCALAR_MAX ? n++ : -1;
  assert (s >= 0 && "hit nscalar limit");
  /* fixme : allocate */
  morton_start (0,0,0)
     
  morton_end ();
  return s;
}

void scalar_free (scalar s)
{
  assert (s>=0); 
  available [navailable++] = s;
}

vector vector_new ()
{
  vector v = (vector) {-1,-1,-1};
  scalar * V = &v.x;
  for (int d=0; d < dimension; ++d)
    V [d] = scalar_new ();
  return v;
}

void vector_free (vector v)
{
  scalar * V = &v.x;
  for (int d=0; d < dimension; ++d)
    scalar_free (V[d]);
}

