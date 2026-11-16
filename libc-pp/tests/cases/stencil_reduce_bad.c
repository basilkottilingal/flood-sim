/* reduction (): its errors, and a variable written but not reduced: a
 * warning, an error with -fopenmp or -D_MPI. Build with -dimension=2. */
#include "stencil.h"

scalar s;
double g;

int main (void)
{
  double a = 0, b = 0;
  int k = 0;
  foreach (reduction (-:a)) a -= s[];
  foreach (reduction (+:k, k)) k += s[] > 0;
  foreach (reduction (+:nosuch)) s[] = 1;
  foreach (reduction (avg:a)) a += s[];
  foreach () { double w = s[]; w += 1; s[] = w; }
  foreach () b += s[];
  foreach () g = s[] > g ? s[] : g;
  foreach () { k++; k++; }
  return (int) (a + b + k);
}
