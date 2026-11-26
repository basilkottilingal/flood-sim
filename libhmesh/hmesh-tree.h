#ifndef _HMESH_TREE_H_
#define _HMESH_TREE_H_

  /*
  .. The manifold in libb-amr's tree: its points sorted along the tree's
  .. Morton curve, and grouped by the leaf they fall in.
  ..
  .. hm_tree_bind (M): the tree's domain, [x0, x0 + h0]^D from
  .. tree_get_origin () and tree_get_size (), becomes the mesh's box
  .. (hm_set_domain ()): from then on a point outside is refused. -1 if a
  .. point is outside already, or if the tree is deeper than the keys can
  .. tell its cells apart: tree_get_depth_max () + log2 (tree_get_block_m ())
  .. must be at most HM_MORTON_BITS (52 in 1D, 32 in 2D, 21 in 3D: in 3D a
  .. depth of 17 with blocks of 16 cells).
  ..
  .. hm_tree_clip (M): a mesh crossing the tree's domain cut to it
  .. (hm_clip (), hmesh.h), which binds it as hm_tree_bind () does; -1 if
  .. the tree is too deep.
  ..
  .. hm_tree_cache (M): the cache, made again if the tree changed
  .. (tree_generation) or the points did (M->generation: added, removed,
  .. moved, renumbered; hm_touch () after writing hm_x ()); bound first if
  .. it was not, or if the tree's domain changed (tree_shift_resize ()).
  .. NULL, with a message, if that fails. Made again from the last order:
  .. sorted already if nothing moved across a cell, by insertion if a few
  .. did.
  ..
  .. The points are sorted to the tree's finest cells, cell_bits bits an
  .. axis (tree_get_depth_max () + log2 (tree_get_block_m ())), which is
  .. all the leaves and their cells need, in (D cell_bits) / 11 radix passes:
  .. two in 2D at depth 6 with blocks of 16. Within a cell the points keep
  .. the last order; hm_compact_morton () sorts them by their whole keys.
  ..
  ..   key[a], order[a]       a < n: the points in use, along the curve:
  ..                          order[a] the point, key[a] its whole key,
  ..                          ascending to the cells
  ..   leaf_of[v]             v < hm_end (M, 0): the index in
  ..                          tree_caches ()->leaves of the leaf v is in;
  ..                          HM_NONE if v is not in use or in no leaf
  ..                          (a pruned part of the domain)
  ..   start[q], stop[q]      q < nleaves: leaf q's points are order[a],
  ..                          start[q] <= a < stop[q] (hm_foreach_leaf_point)
  ..   outside                points in no leaf
  ..   duplicates             points with the whole key of another: closer
  ..                          than size / 2^HM_MORTON_BITS on every axis
  ..                          (about 5e-7 size in 3D)
  ..
  .. A point on a face shared by two leaves goes to the upper one (the key
  .. floors), a point on the domain's upper face to the last cell.
  ..
  .. The cache is the mesh's, freed with it; the pointer stays valid until
  .. the next call.
  ..
  .. hm_compact_morton (M, remap): every level renumbered along the curve:
  .. the points by their keys, a k-simplex by the key of its centroid
  .. (hm_compact_order ()). Then the points' ids are the cache's order,
  .. order[a] = a, and a leaf's points are the ids start[q] ... stop[q] - 1:
  .. its fields, read in memory order. -1 if the cache cannot be made.
  */

  /*
  .. Points in the tree
  .. ------------------
  ..
  ..   hm_tree_locate (M, x)      the index in tree_caches ()->leaves of the
  ..                              leaf x is in, HM_NONE outside the domain or
  ..                              in a pruned part
  ..   hm_tree_cell_size (M, x)   the size of the leaf cell x is in, 0 if
  ..                              none: with -DEULERIAN_LAGRANGIAN,
  ..                              hm_tree_bind () makes it the mesh's cell
  ..                              size (hm_set_cell_size ())
  */

  /*
  .. hm_tree_from_function (M, f): the interface f = 0 (hmesh-implicit.h) in
  .. the tree's domain, on a lattice of half its finest cell (2 m 2^depth a
  .. side), so that every edge is shorter than the cells; then bound to the
  .. tree (hm_tree_bind ()). f is evaluated (2 m 2^depth + 1)^D times. The
  .. number of simplices added, or -1.
  */

  #include "hmesh.h"
  #include "hmesh-morton.h"
  #include "hmesh-implicit.h"

  typedef struct {
    unsigned   tree_generation, mesh_generation;
    int        valid;
    uint32_t   n, cap;          /* points sorted, room for */
    uint64_t * key;
    uint32_t * order;
    uint32_t * leaf_of;
    uint32_t   leaf_cap;        /* room in leaf_of */
    uint32_t * start, * stop;
    uint32_t   nleaves, nleaves_cap;
    uint32_t   outside, duplicates;
    int        cell_bits;       /* the sort's bits an axis */
    /* the leaves' keys and levels, for hm_tree_locate () */
    unsigned   table_generation;
    int        table_valid;
    uint32_t   table_n, table_cap;
    uint64_t * table_key;
    int      * table_level;
    /* the sort's scratch (hm_sort_scratch ()) */
    void     * scratch;
    size_t     scratch_bytes;
  } HmTreeCache;

  int                 hm_tree_bind (Hmesh * M);
  int                 hm_tree_clip (Hmesh * M);
  const HmTreeCache * hm_tree_cache (Hmesh * M);
  int                 hm_compact_morton (Hmesh * M, uint32_t ** remap);
  long                hm_tree_from_function (Hmesh * M, HmFunction f);
  uint32_t            hm_tree_locate (Hmesh * M, const double * x);
  double              hm_tree_cell_size (Hmesh * M, const double * x);

  /* the points v of leaf q */
  #define hm_foreach_leaf_point(C, q, v)                                \
    for (uint32_t a_##v = (C)->start[q], v = 0;                         \
         a_##v < (C)->stop[q] && ((v = (C)->order[a_##v]), 1); a_##v++)

#endif
