@loop foreach_leaf () : sblocks { void *** sblocks = 0; {...} }
@loop foreach_cell (int shape, int halo) : c, shape { int c = halo; {...} }
scalar s;
void f (void)
{
  foreach () s[] = x;
}
