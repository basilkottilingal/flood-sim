#!/bin/sh
# libb-amr's tests, in 1, 2 and 3 dimensions, double and float, under
# AddressSanitizer and UBSan. Run from libb-amr:  sh tests/run.sh
# CWRAP: the cwrap to build with (default ../libc-pp/build/cwrap).
cwrap=${CWRAP:-$(pwd)/../libc-pp/build/cwrap}
lib=$(pwd)
dir=$(mktemp -d)
trap 'rm -rf "$dir"' EXIT
# libutils (bits.h) is a sibling of libb-amr: included as "libutils/bits.h"
src="$lib/grid-scalar.c $lib/morton-curve.c $lib/memindex.c $lib/voxel-memindex.c $lib/nano-vdb.c $lib/pool.c $lib/block-ghosts.c $lib/adapt.c $lib/../libutils/bits.c"
flags="-D_GNU_SOURCE -std=gnu11 -O1 -g -Wall -Wextra -Wno-unused-parameter -fsanitize=address,undefined -I$lib -I$lib/.."
fails=0

run () {   # name, dimension, extra flags, sources...
  name=$1 d=$2 extra=$3; shift 3
  if "$cwrap" -dimension=$d $flags $extra "$@" -lm -o "$dir/t" 2> "$dir/log" &&
     ! grep -v 'mem_destroy\|^In file\|^ *[0-9]* |\|^ *|' "$dir/log" | grep -q 'warning\|error' &&
     out=$("$dir/t" 2>&1); then
    echo "ok   $name, ${d}D ${extra:-double}: $(echo "$out" | tail -1)"
  else
    echo "FAIL $name, ${d}D ${extra:-double}"; cat "$dir/log"; echo "$out" | tail -5
    fails=$((fails + 1))
  fi
}

# the pool, giving memory back either way: MADV_DONTNEED (1) or MADV_FREE (0)
for madv in 1 0; do
  def="-D_BAMR_MADV_DONTNEED_RECOMMENDED_=$madv"
  if gcc -O2 -Wall -Wextra -Werror $def -I"$lib" tests/test_pool.c "$lib/pool.c" -o "$dir/pool" 2> "$dir/pool.log" &&
     "$dir/pool" > "$dir/pool.log" 2>&1; then
    echo "ok   pool $def: $(tail -1 "$dir/pool.log")"; sed 's/^/       /' "$dir/pool.log" | sed -n 2,4p
  else echo "FAIL pool $def"; cat "$dir/pool.log"; fails=$((fails + 1)); fi
done

for prec in "" "-DBAMR_SINGLE"; do
  for d in 1 2 3; do
    run tree   $d "$prec" tests/test_tree.c $src
    run blocks $d "$prec" tests/test_blocks.c "$lib/tree.c" $src
    run ghosts $d "$prec" tests/test_ghosts.c "$lib/tree.c" $src
    run adapt  $d "$prec" tests/test_adapt.c "$lib/tree.c" $src
    run stencil $d "$prec" tests/test_stencil.c "$lib/tree.c" $src
  done
done
# with OpenMP: threads, simd and the offload pragmas (run on the host here)
flags="$flags -fopenmp"
export OMP_NUM_THREADS=4
for d in 2 3; do
  run blocks $d "-fopenmp" tests/test_blocks.c "$lib/tree.c" $src
  run ghosts $d "-fopenmp" tests/test_ghosts.c "$lib/tree.c" $src
  run adapt  $d "-fopenmp" tests/test_adapt.c "$lib/tree.c" $src
  run stencil $d "-fopenmp" tests/test_stencil.c "$lib/tree.c" $src
done
# stored blocks in the aligned layout: rows on 4 Reals
for d in 1 2 3; do
  run blocks $d "-DBAMR_ALIGN=4" tests/test_blocks.c "$lib/tree.c" $src
  run ghosts $d "-DBAMR_ALIGN=4" tests/test_ghosts.c "$lib/tree.c" $src
  run adapt  $d "-DBAMR_ALIGN=4" tests/test_adapt.c "$lib/tree.c" $src
  run stencil $d "-DBAMR_ALIGN=4" tests/test_stencil.c "$lib/tree.c" $src
done
# the level indices side by side: memindex, voxel-memindex, nano-vdb
isrc="$lib/memindex.c $lib/voxel-memindex.c $lib/nano-vdb.c $lib/pool.c"
for d in 1 2 3; do
  run index $d "" tests/test_index.c $isrc
  run index $d "-DNVDB_COMPACT" tests/test_index.c $isrc
done
# the tree on the other level indices (level-index.h): voxel, nano-vdb
for idx in "2" "3" "3 -DNVDB_COMPACT"; do
  for d in 2 3; do
    run tree   $d "-DBAMR_INDEX=$idx" tests/test_tree.c $src
    run blocks $d "-DBAMR_INDEX=$idx" tests/test_blocks.c "$lib/tree.c" $src
    run ghosts $d "-DBAMR_INDEX=$idx" tests/test_ghosts.c "$lib/tree.c" $src
    run adapt  $d "-DBAMR_INDEX=$idx" tests/test_adapt.c "$lib/tree.c" $src
  done
done
# levels are the bits of a uint32_t: BAMR_MAX_DEPTH 32 must not build
if "$cwrap" -dimension=2 $flags -DBAMR_MAX_DEPTH=32 -c "$lib/tree.c" -o "$dir/t.o" 2> "$dir/log"; then
  echo "FAIL tree.c built with BAMR_MAX_DEPTH=32"; fails=$((fails + 1))
elif grep -q "BAMR_MAX_DEPTH must be below 32" "$dir/log"; then
  echo "ok   tree.c refuses BAMR_MAX_DEPTH=32"
else
  echo "FAIL tree.c with BAMR_MAX_DEPTH=32: not the #error"; cat "$dir/log"; fails=$((fails + 1))
fi
# the tree on a pool that gives memory back with MADV_FREE
for d in 2 3; do
  run ghosts $d "-D_BAMR_MADV_DONTNEED_RECOMMENDED_=0" tests/test_ghosts.c "$lib/tree.c" $src
  run adapt  $d "-D_BAMR_MADV_DONTNEED_RECOMMENDED_=0" tests/test_adapt.c "$lib/tree.c" $src
done
[ "$fails" = 0 ] && echo "all passed"
exit "$fails"
