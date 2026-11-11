scalar zb;
vector u;
extern scalar far;
tensor t;
void f (void)
{
  static scalar kept;
  scalar mine = scalar_new ();
  scalar_free (zb);
  scalar_free (u.x);
  scalar_free (kept);
  scalar_free (far);
  tensor_free (t);
  scalar_free (mine);
}
