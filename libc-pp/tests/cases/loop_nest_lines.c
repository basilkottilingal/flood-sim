#include "loop_nest.h"
#include "loop_nest_lines.h"
int main (void)
{
  int total = 0;
  repeat (2)
    range (0, 3)
      grid (2, 2) {
        total += r + c + i + k + nope1;
      }
  int after = nope2;
  bad_body (2)
    total++;
  return total + after + nope3;
}
