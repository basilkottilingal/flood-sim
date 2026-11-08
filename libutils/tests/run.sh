#!/bin/sh
# bits.h in every C standard gcc knows, strict ISO (-pedantic-errors) and
# GNU, with the builtins and with the de Bruijn fallback, all warnings as
# errors, and once each under UBSan and ASan. Run from bits:  sh run.sh
# CC: the compiler (default gcc).
cc=${CC:-gcc}
dir=$(mktemp -d)
trap 'rm -rf "$dir"' EXIT
warn="-Wall -Wextra -Werror -Wconversion -Wsign-conversion -Wshadow -Wcast-qual -Wstrict-prototypes"
fails=0

run () {   # label, flags...
  label=$1; shift
  if $cc $warn "$@" test_bits.c bits.c -o "$dir/t" 2> "$dir/log" && out=$("$dir/t"); then
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
  if $cc -m32 $warn -std=c11 -pedantic-errors $fb -O2 -c bits.c -o "$dir/b32.o" 2> /dev/null &&
     $cc -m32 $warn -std=c11 -pedantic-errors $fb -O2 -c test_bits.c -o "$dir/t32.o" 2> "$dir/log"; then
    echo "ok   -m32 -std=c11 $fb: compiles"
  else
    echo "skip -m32 -std=c11 $fb: no 32 bit headers here"
  fi
done
[ "$fails" = 0 ] && echo "all passed"
exit "$fails"
