#include "stencil.h"
@foreach again () = foreach_leaf () foreach_cell (tree_block_shape (), 0);
scalar s;
void f (void)
{
  foreach () { again () s[] = 0; }
  again (1) s[] = 0;
}
@foreach foreach () = foreach_leaf () foreach_cell (tree_block_shape (), 1);
@foreach foreach () = foreach_leaf () foreach_cell (tree_block_shape (), 0);
