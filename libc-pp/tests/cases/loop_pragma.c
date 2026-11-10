/* A pragma in a loop's body, as _Pragma leaves it: it goes with each use,
 * not where the loop is defined, and lines stay right around it */
@loop scaled (double * a, int n) : x
{
  _Pragma ("GCC ivdep")
  for (int q = 0; q < n; q++) {
    double * x = &a [q];
    {...}
  }
}

int main (void)
{
  double a [4] = { 1, 2, 3, 4 };
  scaled (a, 4)
    *x *= 2;
  scaled (a, 4) *x += 1;
  return (int) (a [0] + a [3]) + nope;
}
