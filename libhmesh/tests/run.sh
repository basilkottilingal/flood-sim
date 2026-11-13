#!/bin/sh
# libhmesh's tests, in 1, 2 and 3 dimensions, double and float, under
# AddressSanitizer and UBSan, failing on any warning; then once at -O2 for
# the time of a removal and an addition (QUICK=1 leaves that out). Run from
# libhmesh:
#   sh tests/run.sh
# CWRAP: the cwrap to build with (default ../libc-pp/build/cwrap, else cwrap
# on the PATH).
cwrap=${CWRAP:-$(pwd)/../libc-pp/build/cwrap}
[ -x "$cwrap" ] || cwrap=cwrap
lib=$(pwd)
dir=$(mktemp -d)
trap 'rm -rf "$dir"' EXIT
src="$lib/hmesh.c $lib/../libutils/bits.c"
warn="-Wall -Wextra -Wshadow -Wconversion -Wsign-conversion -Wno-float-conversion"
flags="-D_GNU_SOURCE -std=gnu11 $warn -I$lib -I$lib/.."
san="-O1 -g -fsanitize=address,undefined -fno-sanitize-recover=all"
fails=0

run () {   # name, dimension, flags...
  name=$1 d=$2
  shift 2
  out=
  if "$cwrap" -dimension=$d $flags "$@" tests/test_hmesh.c $src -lm -o "$dir/t" 2> "$dir/log" &&
     ! grep -q 'warning\|error' "$dir/log" &&
     out=$("$dir/t" 2>&1); then
    echo "ok   $name, ${d}D: $(echo "$out" | tail -1)"
    echo "$out" | sed '$d' | sed 's/^/     /'
  else
    echo "FAIL $name, ${d}D"; cat "$dir/log"; echo "$out" | tail -12
    fails=$((fails + 1))
  fi
}

for d in ${DIMS:-1 2 3}; do
  run "double, sanitizers" $d $san
  run "float, sanitizers" $d $san -DBAMR_SINGLE
  [ -n "$QUICK" ] || run "double, -O2" $d -O2
done
[ "$fails" = 0 ] && echo "all passed"
exit "$fails"
