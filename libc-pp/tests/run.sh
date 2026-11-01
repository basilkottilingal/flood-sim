#!/bin/sh
# End-to-end checks of the wrapper. Run from the libc-pp directory: make test
cwrap=$(pwd)/build/cwrap
tests=$(pwd)/tests
dir=$(mktemp -d)
trap 'rm -rf "$dir"' EXIT
cd "$dir" || exit 1
fails=0

check() {
  if [ "$2" = 0 ]; then
    echo "ok   $1"
  else
    echo "FAIL $1"
    fails=$((fails + 1))
  fi
}

cat > hello.c <<'EOF'
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#ifndef VALUE
#define VALUE 0
#endif
struct point { double x, y; };
int main(void)
{
  struct point p = { .x = 3.0, .y = 4.0 };
  printf("%.1f %d %s\n", sqrt(p.x * p.x + p.y * p.y), VALUE, "a\"b");
  return 0;
}
EOF

"$cwrap" -Wall -DVALUE=42 hello.c -o hello -lm
[ "$(./hello)" = '5.0 42 a"b' ]
check "compile, link and run with stdio, math and -l" $?

CWRAP_KEEP=1 "$cwrap" -c hello.c -o hello.o 2> keep.log
kept=$(sed -n 's/^cwrap: kept //p' keep.log)
cmp -s "$kept/0-hello.i" "$kept/0-hello.t.i"
check "translator output is byte-identical to its input" $?
rm -rf "$kept"

cat > bad.c <<'EOF'
#include <stdio.h>

int main(void)
{
  return missing;
}
EOF
"$cwrap" -c bad.c -o bad.o 2> bad.log
[ $? -ne 0 ] && grep -q '^bad.c:5:' bad.log
check "gcc errors point at the original file and line" $?

printf 'int twice(int x) { return 2 * x; }\n' > a.c
printf 'int twice(int);\nint main(void) { return twice(4) != 8; }\n' > b.c
"$cwrap" a.c b.c -o two && ./two
check "several .c files in one command" $?

"$cwrap" -c a.c && "$cwrap" -c b.c -o b2.o && "$cwrap" a.o b2.o -o two2 && ./two2
check "separate compile, then link objects" $?

printf '#define N 7\n' > dep.h
printf '#include "dep.h"\nint n = N;\n' > dep.c
"$cwrap" -MD -c dep.c -o dep.o
grep -q '^dep.o:' dep.d && grep -q 'dep.h' dep.d
check "-MD writes a dependency file for the object" $?

[ "$("$cwrap" -E -DN=9 -P - <<'EOF'
N
EOF
)" = 9 ]
check "-E passes straight through to gcc" $?

printf 'int x = 1 @ 2;\n' > stray.c
"$cwrap" -c stray.c 2> stray.log
[ $? -ne 0 ] && grep -q '^stray.c:1: error' stray.log
check "translator errors point at the original file and line" $?

cp "$tests/grammar.c" grammar.c
CWRAP_KEEP=1 "$cwrap" -std=gnu11 -w -c grammar.c -o grammar.o 2> keep.log
rc=$?
kept=$(sed -n 's/^cwrap: kept //p' keep.log)
[ "$rc" = 0 ] && cmp -s "$kept/0-grammar.i" "$kept/0-grammar.t.i"
check "the full grammar accepts tests/grammar.c, byte for byte" $?
rm -rf "$kept"

printf 'int in_header;\n' > loc.h
printf '#include "loc.h"\nint\nin_source;\n' > loc.c
CWRAP_DUMP_AST=1 "$cwrap" -c loc.c -o loc.o 2> loc.log
grep -q "'in_header' loc.h:1" loc.log && grep -q "'in_source' loc.c:3" loc.log
check "each leaf records its file and line" $?

printf 'int x = 1 + 2 * 3;\n' > prec.c
CWRAP_DUMP_AST=1 "$cwrap" -c prec.c -o prec.o 2> prec.log
[ "$(grep -nE 'ADDITIVE|MULTIPLICATIVE' prec.log | sed 's/.*[0-9]: *//')" = "ADDITIVE
MULTIPLICATIVE" ]
check "operators nest by precedence" $?

printf 'typedef int T;\nvoid f (void) { int T; T * 2; T = 1; }\n' > hide.c
"$cwrap" -w -c hide.c -o hide.o
check "an inner variable hides a typedef" $?

[ "$fails" = 0 ] && echo "all passed"
exit "$fails"
