#!/bin/sh
# libutils' tests. Run from anywhere:  sh tests/run.sh
# bits.h in every C standard gcc knows, strict ISO (-pedantic-errors) and
# GNU, with the builtins and with the de Bruijn fallback, all warnings as
# errors, and once each under UBSan and ASan. profile.h: on one process
# with and without -DPROFILE_ENABLE, its outputs checked (python3 reads
# the trace's JSON), under ASan and UBSan; then, if mpicc and mpirun are
# here, on 1, 3 and 4 MPI processes. CC: the compiler (default gcc).
cc=${CC:-gcc}
here=$(cd "$(dirname "$0")" && pwd)
lib=$(cd "$here/.." && pwd)
dir=$(mktemp -d)
trap 'rm -rf "$dir"' EXIT
warn="-Wall -Wextra -Werror -Wconversion -Wsign-conversion -Wshadow -Wcast-qual -Wstrict-prototypes"
fails=0

run () {   # label, flags...
  label=$1; shift
  if $cc $warn "$@" -I"$lib" "$here/test_bits.c" "$lib/bits.c" -o "$dir/t" \
       2> "$dir/log" && out=$("$dir/t"); then
    echo "ok   $label: $out"
  else
    echo "FAIL $label"; cat "$dir/log"; echo "$out" | tail -5
    fails=$((fails + 1))
  fi
}

for std in c99 c11 c17 c2x gnu99 gnu11 gnu17 gnu2x c89 gnu89; do
  case $std in
    c89|gnu89) strict="-Wdeclaration-after-statement" ;;   # long long: not ISO C90
    *)         strict="-pedantic-errors" ;;
  esac
  for fb in "" "-DBITS_NO_BUILTINS"; do
    run "-std=$std $fb -O2" -std=$std $strict $fb -O2
  done
  run "-std=$std -O1 UBSan, ASan, de Bruijn" -std=$std $strict -O1 -g \
      -DBITS_NO_BUILTINS -fsanitize=undefined,address -fno-sanitize-recover=all
done
# the instructions, where this machine has them
run "-std=gnu11 -O2 -march=native" -std=gnu11 -O2 -march=native
# a 32 bit target, compiled only (long is 32 bits there): if it can be
for fb in "" "-DBITS_NO_BUILTINS"; do
  if $cc -m32 $warn -std=c11 -pedantic-errors $fb -O2 -c "$lib/bits.c" -o "$dir/b32.o" 2> /dev/null &&
     $cc -m32 $warn -std=c11 -pedantic-errors $fb -O2 -I"$lib" -c "$here/test_bits.c" -o "$dir/t32.o" 2> "$dir/log"; then
    echo "ok   -m32 -std=c11 $fb: compiles"
  else
    echo "skip -m32 -std=c11 $fb: no 32 bit headers here"
  fi
done

# profile.h, one process: the test, then its files
pwarn="-Wall -Wextra -Werror -Wshadow -Wconversion"
prof () {   # label, flags...
  label=$1; shift
  rm -f "$dir"/p.*  "$dir"/p_*
  if $cc -std=gnu11 $pwarn "$@" -I"$lib" "$here/test_profile.c" "$lib/profile.c" \
       -lpthread -lm -o "$dir/p" 2> "$dir/log" &&
     out=$(cd "$dir" && PROFILE=table,flame,trace PROFILE_OUT="$dir/p" ./p) &&
     check_outputs "$@"; then
    echo "ok   profile $label: $(echo "$out" | tail -1)"
  else
    echo "FAIL profile $label"; cat "$dir/log"; echo "$out" | tail -5
    fails=$((fails + 1))
  fi
}
check_outputs () {   # timed: table, folded stacks and a trace that parses
  case "$*" in *PROFILE_ENABLE*) ;; *) [ ! -e "$dir/p.txt" ]; return ;; esac
  grep -q "^outer " "$dir/p.txt" && grep -q "^outer;inner [0-9]" "$dir/p.folded" &&
    python3 -c "import json, sys; e = json.load (open (sys.argv [1]))['traceEvents'];
sys.exit (0 if len (e) > 100 and len ({x.get ('tid') for x in e}) >= 4 else 1)" \
      "$dir/p.trace.json"
}
prof "-O2" -O2 -DPROFILE_ENABLE
prof "-O2, not enabled" -O2
prof "-O1 ASan, UBSan" -O1 -g -DPROFILE_ENABLE -fsanitize=address,undefined \
  -fno-sanitize-recover=all

# profile.h on MPI processes (Open MPI's mpirun takes --oversubscribe)
if command -v mpicc > /dev/null && command -v mpirun > /dev/null; then
  if mpicc -std=gnu11 $pwarn -D_MPI -DPROFILE_ENABLE -I"$lib" \
       "$here/test_profile_mpi.c" "$lib/profile.c" -lpthread -lm -o "$dir/pm" \
       2> "$dir/log"; then
    for np in 1 3 4; do
      if out=$(cd "$dir" && PROFILE=table,flame,trace PROFILE_OUT="$dir/pm" \
                 timeout 120 mpirun --oversubscribe -np $np -x PROFILE \
                 -x PROFILE_OUT ./pm 2>&1); then
        echo "ok   profile MPI np $np: $(echo "$out" | grep checks)"
      else
        echo "FAIL profile MPI np $np"; echo "$out" | tail -8
        fails=$((fails + 1))
      fi
    done
  else
    echo "FAIL profile MPI: build"; cat "$dir/log"; fails=$((fails + 1))
  fi
else
  echo "skip profile MPI: no mpicc or mpirun"
fi
[ "$fails" = 0 ] && echo "all passed"
exit "$fails"
