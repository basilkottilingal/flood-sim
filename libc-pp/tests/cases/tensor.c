#include <stdio.h>
scalar s; int kept_by_hand (void);
vector v;
tensor t;
static int sum (void)
{
  vector u;
  tensor w;
  static vector kept;
  return u.x + w.y.y + kept.y - kept.x;
}
int main (void)
{
  int in_use = scalar_count ();
  int inner = sum ();
  printf ("%d %d %d %d %d %d %d %d\n", s, v.x, v.y, t.x.x, t.x.y, t.y.x, t.y.y,
          n_global_scalars);
  printf ("%d %d %d %d\n", in_use, inner, scalar_count (), kept_by_hand ());
  t.x.y = 1;
  v = v;
  t.y = v;
  return 0;
}
/* by hand: ids that last until freed */
int kept_by_hand (void)
{
  vector a = vector_new ();
  tensor b = tensor_new ();
  int last = b.y.y, n = scalar_count ();
  tensor_free (b);
  vector_free (a);
  return last * 100 + n - scalar_count ();
}
