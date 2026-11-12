#!/bin/sh
# liblinear-algebra's tests, in 1, 2 and 3 dimensions, double and float,
# under AddressSanitizer and UBSan, failing on any warning. Run from
# liblinear-algebra:  sh tests/run.sh
# CWRAP: the cwrap to build with (default ../libc-pp/build/cwrap).
cwrap=${CWRAP:-$(pwd)/../libc-pp/build/cwrap}
lib=$(pwd)
amr=$lib/../libb-amr
dir=$(mktemp -d)
trap 'rm -rf "$dir"' EXIT
src="$lib/multigrid.c $amr/tree.c $amr/grid-scalar.c $amr/memindex.c $amr/voxel-memindex.c $amr/nano-vdb.c $amr/pool.c $amr/block-ghosts.c $amr/adapt.c $lib/../libutils/bits.c"
flags="-D_GNU_SOURCE -std=gnu11 -O1 -g -Wall -Wextra -Wno-unused-parameter -fsanitize=address,undefined -I$lib -I$amr -I$lib/.."
fails=0

run () {   # name, dimension, extra flags, the test
  name=$1 d=$2 extra=$3 test=$4
  out=
  if "$cwrap" -dimension=$d $flags $extra "$test" $src -lm -o "$dir/t" 2> "$dir/log" &&
     ! grep -v 'mem_destroy\|^In file\|^ *[0-9]* |\|^ *|' "$dir/log" | grep -q 'warning\|error' &&
     out=$("$dir/t" 2>&1); then
    echo "ok   $name, ${d}D ${extra:-double}: $(echo "$out" | tail -1)"
    echo "$out" | sed '$d' | sed 's/^/     /'
  else
    echo "FAIL $name, ${d}D ${extra:-double}"; cat "$dir/log"; echo "$out" | tail -8
    fails=$((fails + 1))
  fi
}

for prec in "" "-DBAMR_SINGLE"; do
  for d in 1 2 3; do
    run poisson $d "$prec" tests/test_poisson.c
  done
done
[ "$fails" = 0 ] && echo "all passed"
exit "$fails"
