@loop foreach_leaf () : sblocks { void *** sblocks = 0; {...} }
@loop foreach_cell (int shape, int halo) : c, shape { int c = halo; {...} } @foreach foreach () = foreach_leaf () foreach_cell (0, 0);
scalar s;
void f (void)
{
  foreach () s[] = x;
}
