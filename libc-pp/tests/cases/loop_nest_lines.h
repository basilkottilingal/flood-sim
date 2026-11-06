/* a loop with an error in its body, after a loop in it */
@loop bad_body (int n) : k
{
  repeat (n) {
    (void) k;
  }
  int z = missing_name;
  repeat (z)
    {...}
}
