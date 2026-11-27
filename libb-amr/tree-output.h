#ifndef _BAMR_TREE_OUTPUT_H_
#define _BAMR_TREE_OUTPUT_H_

  /*
  .. The tree's leaves and scalars written to a file, .bamr: every leaf's
  .. index and, for each scalar given, every cell's value, exactly (the
  .. Reals as they are), in the order of the Morton curve; and read back
  .. for a restart.
  ..
  .. tree_output (name, list, n): writes the scalars list [0 ... n - 1] to
  .. name, with ".bamr" added when name has no extension. Under MPI every
  .. process calls it: each writes its own leaves at their place, by prefix
  .. sums of the counts (MPI_Exscan) and MPI-IO, and the file is the same,
  .. byte for byte, whatever the number of processes (the processes' pieces
  .. are runs of the curve, in order). A tree not partitioned yet
  .. (tree_replicated (), whole on every process) is written once, by the
  .. first process. 0, or -1 on every process if one
  .. failed: n < 0, list NULL with n > 0, a scalar not in use, a scalar with
  .. no name (scalar_new ()'s: only those cwrap declares are named), two
  .. with one name, the file not written (then it may be left partly
  .. written).
  ..
  .. The names are the scalars' own (scalar_name (): "h", "u.x"), each
  .. once: a restart finds each field by its name (scalar_find ()).
  ..
  .. tree_restart (path): the tree made again from a .bamr file, and its
  .. scalars filled. The program makes its tree as ever (tree (depthmax,
  .. blockdepth)), then calls it, every process under MPI:
  ..
  ..   1. the file read and checked, before anything changes: the header
  ..      (version, byte order, codec 0), the dimension and m the tree's,
  ..      the size the sections need; the leaves' levels at most the
  ..      tree's depth max, their indices in range, along the curve
  ..      (strictly), covering the domain (an integer sum)
  ..   2. the tree replaced: the same depth max and blocks, the file's
  ..      domain, refined level by level to the file's leaves (the file's
  ..      tree is 2:1 balanced, so the balancing refinements stay within
  ..      it), every block 0; under MPI made whole, then partitioned
  ..      (tree_mpi_partition (): the number of processes need not be the
  ..      writer's)
  ..   3. each name in the file that is a global scalar in use
  ..      (scalar_find (), an id below n_global_scalars) gets its cells,
  ..      converted if the file's Reals are the other precision, its
  ..      ghosts marked stale (ghosts_mark_dirty ()); other names, local
  ..      scalars' among them, are left
  ..
  .. The number of scalars filled, or -1 on every process: no tree, a file
  .. that cannot be read or fails a check of step 1 (the program's tree
  .. then as it was), or one whose leaves the refinement does not give
  .. (not 2:1 balanced: the tree is then replaced, not usable). The
  .. simulation's time and step are not in the file.
  ..
  .. The file, little-endian (the host must be: -1 otherwise):
  ..
  ..   offset  size     what
  ..   0       8        "BAMRGRID"
  ..   8       4        version (1)
  ..   12      4        0x01020304, read as written: the byte order
  ..   16      4        dimension D (1, 2, 3)
  ..   20      4        m, cells a block side (blocks of m^D cells)
  ..   24      4        the size of a Real: 8 (double) or 4 (-DBAMR_SINGLE)
  ..   28      4        n, the scalars
  ..   32      8        N, the leaves
  ..   40      24       the domain's origin x0, y0, z0 (doubles; 0 past D)
  ..   64      8        its size h0 (a double): [x0, x0 + h0]^D
  ..   72      4        codec: 0, none (the field left for compression)
  ..   76      4        0 (reserved)
  ..   80               the scalars' names: n of (u32 length, the bytes,
  ..                    no NUL), then zeros to a multiple of 8: H bytes in
  ..                    all from offset 0
  ..   H       16 N     the leaves along the curve: (level, i, j, k) as
  ..                    int32 each (j, k 0 past D)
  ..   H + 16 N         scalar s, s = 0 ... n - 1: N m^D Reals, leaf by
  ..                    leaf in the order above, a leaf's cells i fastest,
  ..                    then j, then k (its interior: no ghost)
  ..
  .. Scalar s starts at H + 16 N + s N m^D (size of a Real). Leaf (level,
  .. i, j, k) is the box x0 + h0 2^-level ([i, i + 1] x [j, j + 1] x [k, k +
  .. 1]), its cell (a, b, c) of size h0 2^-level / m. libb-amr/tools/
  .. bamr.py reads it in Python (the standard library alone).
  */

  #include "tree.h"

  int tree_output (const char * name, const scalar * list, int n);
  int tree_restart (const char * path);

#endif
