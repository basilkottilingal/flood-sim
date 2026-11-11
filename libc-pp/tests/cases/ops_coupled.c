#include <stdio.h>

typedef struct { int s; double * out; } Pr;
enum { GHOST_PROLONGATION, GHOST_RESTRICTION, GHOST_BOUNDARY };

/* what libb-amr's ghost_ops_set () is given, kept */
static int last_kind, last_op, last_n;
static scalar last_uses [4];
static int ghost_ops_set (scalar s, int kind, int op, const scalar * uses, int n)
{
  (void) s;
  last_kind = kind, last_op = op, last_n = n;
  for (int k = 0; k < n; k++)
    last_uses [k] = uses [k];
  return 0;
}

scalar zb;
vector u;

@prolongation plain (const Pr * pr) { pr->out [0] = 1; }
@prolongation lake (const Pr * pr) (scalar a, scalar b) { pr->out [0] = 10 * a + b; }
@boundary b_one (int s) { (void) s; }

int main (void)
{
  double out [1];
  Pr pr = { 0, out };
  scalar uses [2] = { 3, 4 };
  @prolongation (lake, uses) (&pr);           /* 10 * 3 + 4 */
  printf ("%g ", out [0]);
  @prolongation (plain, uses) (&pr);
  printf ("%g ", out [0]);
  printf ("%d %d | ", @prolongation uses (lake), @prolongation uses (plain));
  scalar h = scalar_new ();
  ghost_ops (h, prolongation (lake, zb, u.y));
  printf ("%d %d %d %d ", last_kind, last_op, last_n,
          last_uses [0] == zb && last_uses [1] == u.y);
  ghost_ops (h, boundary (b_one));
  printf ("%d %d %d\n", last_kind, last_op, last_n);
  return 0;
}
