#!/bin/sh
# bench_vof.c in 2D and 3D, optimised, no sanitizers. Run from libvof:
#   sh tests/bench.sh [extra flags]   (e.g. -march=native, -DBAMR_SINGLE)
cwrap=${CWRAP:-$(pwd)/../libc-pp/build/cwrap}
lib=$(pwd)
amr=$lib/../libb-amr
dir=$(mktemp -d)
trap 'rm -rf "$dir"' EXIT
src="$lib/vof.c $lib/vof-ghosts.c $amr/tree.c $amr/grid-scalar.c $amr/memindex.c $amr/voxel-memindex.c $amr/nano-vdb.c $amr/pool.c $amr/adapt.c $lib/../libutils/bits.c"
for d in ${DIMS:-2 3}; do
  "$cwrap" -dimension=$d -D_GNU_SOURCE -std=gnu11 -O3 -DNDEBUG "$@" -I"$lib" -I"$amr" -I"$lib/.." \
    tests/bench_vof.c $src -lm -o "$dir/b" 2> "$dir/log" || { cat "$dir/log"; exit 1; }
  "$dir/b"
done
