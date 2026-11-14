#!/bin/sh
# libhmesh's tests, in 1, 2 and 3 dimensions, double and float, under
# AddressSanitizer and UBSan, failing on any warning; then once at -O2 for
# the times (QUICK=1 leaves that out). Run from libhmesh:
#   sh tests/run.sh
#   test_hmesh.c   the mesh alone
#   test_tree.c    the mesh in libb-amr's tree (hmesh-tree.h)
#   test_clip.c    the mesh cut to a box (hm_clip ())
# CWRAP: the cwrap to build with (default ../libc-pp/build/cwrap, else cwrap
# on the PATH).
cwrap=${CWRAP:-$(pwd)/../libc-pp/build/cwrap}
[ -x "$cwrap" ] || cwrap=cwrap
lib=$(pwd)
amr=$lib/../libb-amr
dir=$(mktemp -d)
trap 'rm -rf "$dir"' EXIT
src="$lib/hmesh.c $lib/../libutils/bits.c"
csrc="$src $lib/hmesh-clip.c"
tsrc="$csrc $lib/hmesh-morton.c $lib/hmesh-tree.c $amr/tree.c $amr/grid-scalar.c $amr/memindex.c $amr/voxel-memindex.c $amr/nano-vdb.c $amr/pool.c $amr/block-ghosts.c $amr/adapt.c"
warn="-Wall -Wextra -Wshadow -Wconversion -Wsign-conversion -Wno-float-conversion"
flags="-D_GNU_SOURCE -std=gnu11 -I$lib -I$lib/.."
san="-O1 -g -fsanitize=address,undefined -fno-sanitize-recover=all"
fails=0

run () {   # name, dimension, test, sources, flags...
  name=$1 d=$2 test=$3 sources=$4
  shift 4
  out=
  if "$cwrap" -dimension=$d $flags "$@" "$test" $sources -lm -o "$dir/t" 2> "$dir/log" &&
     ! grep -v 'mem_destroy\|^In file\|^ *[0-9]* |\|^ *|' "$dir/log" | grep -q 'warning\|error' &&
     out=$("$dir/t" 2>&1); then
    echo "ok   $name, ${d}D: $(echo "$out" | tail -1)"
    echo "$out" | sed '$d' | sed 's/^/     /'
  else
    echo "FAIL $name, ${d}D"; cat "$dir/log"; echo "$out" | tail -12
    fails=$((fails + 1))
  fi
}

for d in ${DIMS:-1 2 3}; do
  # the mesh alone, every warning of its own
  run "mesh, double, sanitizers" $d tests/test_hmesh.c "$src" $warn $san
  run "mesh, float, sanitizers" $d tests/test_hmesh.c "$src" $warn $san -DBAMR_SINGLE
  [ -n "$QUICK" ] || run "mesh, double, -O2" $d tests/test_hmesh.c "$src" $warn -O2
  run "clip, double, sanitizers" $d tests/test_clip.c "$csrc" $warn $san
  run "clip, float, sanitizers" $d tests/test_clip.c "$csrc" $warn $san -DBAMR_SINGLE
  # in the tree: libb-amr's sources with its own warnings (-Wall -Wextra)
  run "tree, double, sanitizers" $d tests/test_tree.c "$tsrc" -Wall -Wextra -Wno-unused-parameter -I$amr $san
  run "tree, float, sanitizers" $d tests/test_tree.c "$tsrc" -Wall -Wextra -Wno-unused-parameter -I$amr $san -DBAMR_SINGLE
  [ -n "$QUICK" ] || run "tree, double, -O2" $d tests/test_tree.c "$tsrc" -Wall -Wextra -Wno-unused-parameter -I$amr -O2
done
[ "$fails" = 0 ] && echo "all passed"
exit "$fails"
