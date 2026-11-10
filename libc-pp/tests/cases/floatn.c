/* gcc's floating keywords, as <complex.h> declares with _GNU_SOURCE */
#include <complex.h>
#include <tgmath.h>

_Float32 a = 1;
_Complex _Float64 z;
_Float32x w;

int main (void)
{
  float x = 4;
  _Complex _Float32 y = 0;
  return (int) (sqrt (x) + creal (z) + creal (y) + a + w) == 3 ? 0 : 1;
}
