#!/bin/sh
# libpinn's tests, in 1 and 2 dimensions, double and float, with every
# backend Python can import (torch, jax); a missing one is skipped. Run
# from libpinn:  sh tests/run.sh
# CWRAP:      the cwrap to build with (default ../libc-pp/build/cwrap)
# PYTHONPATH: where torch and jax are, if not in the system's site-packages
cwrap=${CWRAP:-$(pwd)/../libc-pp/build/cwrap}
py=${PYTHON:-python3}
lib=$(pwd)
amr=$lib/../libb-amr
dir=$(mktemp -d)
trap 'rm -rf "$dir"' EXIT
src="$amr/tree.c $amr/grid-scalar.c $amr/memindex.c $amr/pool.c $amr/block-ghosts.c $amr/adapt.c $lib/../libutils/bits.c"
flags="-D_GNU_SOURCE -std=gnu11 -O1 -g -Wall -Wextra -Wno-unused-parameter -I$lib -I$amr -I$lib/.."
pyc="$($py-config --includes)"
pyl="$($py-config --ldflags --embed)"
export PYTHONPATH="$lib/tests${PYTHONPATH:+:$PYTHONPATH}"
export TMPDIR="$dir"
fails=0

# the bridge, once, by the C compiler: it includes Python.h
if ! cc $flags $pyc -DPINN_PYDIR="\"$lib/python\"" -c "$lib/pinn.c" -o "$dir/pinn.o" 2> "$dir/log" ||
   grep -q 'warning' "$dir/log"; then
  echo "FAIL pinn.c"; cat "$dir/log"; exit 1
fi

backends=""
for b in torch jax; do
  if $py -c "import $b" 2> /dev/null; then backends="$backends $b"
  else echo "skip $b: python cannot import it"; fi
done
[ -n "$backends" ] || { echo "no backend: install torch or jax"; exit 1; }

for prec in "" "-DBAMR_SINGLE"; do
  for d in 1 2; do
    if ! "$cwrap" -dimension=$d $flags $prec tests/test_pinn.c "$lib/pinn-bamr.c" $src \
         "$dir/pinn.o" $pyl -lm -o "$dir/t" 2> "$dir/log" ||
       grep -v 'mem_destroy\|^In file\|^ *[0-9]* |\|^ *|' "$dir/log" | grep -q 'warning\|error'; then
      echo "FAIL build, ${d}D ${prec:-double}"; cat "$dir/log"; fails=$((fails + 1)); continue
    fi
    for b in $backends; do
      if out=$("$dir/t" $b 2>&1); then
        echo "ok   $b, ${d}D ${prec:-double}: $(echo "$out" | tail -1)"
        echo "$out" | grep '^  ' | sed 's/^/     /'
      else
        echo "FAIL $b, ${d}D ${prec:-double}"; echo "$out" | tail -30
        fails=$((fails + 1))
      fi
    done
  done
done
[ "$fails" = 0 ] && echo "all passed"
exit "$fails"
