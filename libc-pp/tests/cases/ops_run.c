#include <stdio.h>

typedef struct { int s; double * out; } Pr;

@prolongation p_twice (const Pr * pr)
{
  pr->out [pr->s] = 2 * pr->s;
}

@prolongation p_neg (const Pr * pr) { pr->out [pr->s] = - pr->s; }

@restriction r_one (int s, double * out) { out [s] = 1; }

@prolongation p_neg;                   /* a declaration: the same int */

static void prolong (int id, const Pr * pr)
{
  @prolongation (id) (pr);
}

int main (void)
{
  double out [3] = { 0, 0, 0 };
  Pr a = { 1, out }, b = { 2, out };
  int op = p_neg;
  prolong (p_twice, &a);
  prolong (op, &b);
  @restriction (r_one) (0, out);
  printf ("%d %d %d | %g %g %g\n", p_twice, p_neg, r_one, out [0], out [1], out [2]);
  return 0;
}
