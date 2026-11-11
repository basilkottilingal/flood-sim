@boundary b_one (int s) { (void) s; }
@boundary b_one (int s) { (void) s; }
void f (int id)
{
  @restriction (id) (0);
  int x = nope;
}
