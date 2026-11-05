BITREE        l
1,            0    00000000001  {1}            (i,l) -> 2^l + i
1+2 = 3,      1    00000000011  {2,3}
3+4 = 7,      2    00000000111  {4,5,6,7}
7+8 = 15      3    00000001111  {8, ... 15}



QUADTREE                        [starting @  
1,            0    00000000001  {1}             (i,j,l)->
1+4 = 5,      1    00000000101  {2,3,4,5}
5 + 16 = 21,  2    00000010101  {6, .., 21}
21 + 64 = 85, 3    00001010101  {22, .., 85}

OCTREE
1,            0    00000000001  (i,j,k,l) -> (1 << 3*l) | (k << 2*l) | (j<<l) | i
1+8 = 9,      1    00000001001         verify (l=1) ??  {(1<<3)|(k<<2)|(j<<1)|i} for [0,1]^3 = {8
9+64 = 73,    2    00001001001


0,   00000000001 
1,   00000000101
2,   00000010101
3,   00001010101

/*---------------------- Starting index -------------- */

---BITREE--
2^n                                                  or
(1ULL << n)


---QUADTREE--
(4^n-1)/3 + 1                                        or
(4^n-4)/3                                            or
((2ULL << (2 * n)) - 1) & 0x5555555555555555ULL,     or
((1ULL << (2 * n + 2)) - 1) / 3,                    

---OCTREE--
(8^n-1)/7 + 1                                        or
(8^n-8)/7                                            or
((2ULL << (3 * n)) - 1) & 0x9249249249249249ULL      or
((1ULL << (3 * n + 3)) - 1) / 7





Mapping (quadtree traversal . Morton)
{f(0) = 2, f(2) = 1, f(1)  = 3, f(3) = 0} //cyclic version
f(x) = (x ^ (x >> 1)) ^ 2

Mapping (octree traversal) [note that z > y > x in the axis preference in this morton]
{f(0) = 4, f(4) = 2, f(2) = 6, f(6) = 1, f(1) = 5, f(5) = 3, f(3) = 7, f(7) = 0}
can be mapped by
f(x) = (0x5AFAC >> (3*x)) & 7 


  #if dimension == 1
    # define IJK           i
    # define IJKargs       i
    # define IJKarray      [i]
  #elif dimension == 2
    # define IJK           i, j
    # define IJKargs       i, int j
    # define IJKarray      [i][j]
  #else
    # define IJK           i, j, k
    # define IJKargs       i, int j, int k
    # define IJKarray      [i][j][k]
  #endif
  #define IJKL             IJK, level
  #define IJKLargs         IJKargs, int level
  #define IJKEXPAND(...)   __VA_ARGS__
  #define tree_is_allocated(IJKL)                                             \
    (level <= BAMR_MAX_DEPTH && mem_allocated (&grid_info_layers [level], IJK))
  #define tree_node(IJKL)                                                     \
    ((BAMRNode *) mem_data (&grid_info_layers [level], IJK))
  #define mem_data(m,IJK) ((m)->b IJKarray)
  #define tree_has_children(n)                                                 \
    ((n->flags & BAMR_NODE_CHILDREN) &&                                        \
     !(n->flags & BAMR_NODE_PRUNED))
  #define BAMR_NODE_ENV()                                                     \
      BAMRNode * node = tree_is_allocated (IJKL) ?  tree_node (IJKL) : NULL;   \
      (void) node
