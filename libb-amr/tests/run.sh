#!/bin/sh
# libb-amr's tests, in 1, 2 and 3 dimensions, double and float, under
# AddressSanitizer and UBSan. Run from libb-amr:  sh tests/run.sh
# CWRAP: the cwrap to build with (default ../libc-pp/build/cwrap).
cwrap=${CWRAP:-$(pwd)/../libc-pp/build/cwrap}
lib=$(pwd)
dir=$(mktemp -d)
trap 'rm -rf "$dir"' EXIT
src="$lib/grid-scalar.c $lib/morton-curve.c $lib/memindex.c $lib/pool.c $lib/block-ghosts.c"
flags="-D_GNU_SOURCE -std=gnu11 -O1 -g -Wall -Wextra -Wno-unused-parameter -fsanitize=address,undefined -I$lib"
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

gcc -O2 -Wall -Wextra -I"$lib" tests/test_pool.c "$lib/pool.c" -o "$dir/pool" && "$dir/pool" > "$dir/pool.log"
if [ $? = 0 ]; then echo "ok   pool: $(tail -1 "$dir/pool.log")"; sed 's/^/       /' "$dir/pool.log" | head -3
else echo "FAIL pool"; cat "$dir/pool.log"; fails=$((fails + 1)); fi

for prec in "" "-DBAMR_SINGLE"; do
  for d in 1 2 3; do
    run tree   $d "$prec" tests/test_tree.c $src
    run blocks $d "$prec" tests/test_blocks.c "$lib/tree.c" $src
    run ghosts $d "$prec" tests/test_ghosts.c $src
  done
done
# with OpenMP: threads, simd and the offload pragmas (run on the host here)
flags="$flags -fopenmp"
export OMP_NUM_THREADS=4
for d in 2 3; do
  run blocks $d "-fopenmp" tests/test_blocks.c "$lib/tree.c" $src
  run ghosts $d "-fopenmp" tests/test_ghosts.c $src
done
[ "$fails" = 0 ] && echo "all passed"
exit "$fails"
