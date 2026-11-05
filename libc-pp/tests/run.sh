#!/bin/sh
# End-to-end checks of the wrapper. Run from the libc-pp directory: make test
# The programs they build are in tests/cases, copied into a scratch directory.
cwrap=$(pwd)/build/cwrap
scope_test=$(pwd)/build/scope_test
tests=$(pwd)/tests
dir=$(mktemp -d)
trap 'rm -rf "$dir"' EXIT
cd "$dir" || exit 1
cp "$tests"/cases/* . || exit 1
fails=0

check() {
  if [ "$2" = 0 ]; then
    echo "ok   $1"
  else
    echo "FAIL $1"
    fails=$((fails + 1))
  fi
}

"$cwrap" -Wall -DVALUE=42 hello.c -o hello -lm
[ "$(./hello)" = '5.0 42 a"b' ]
check "compile, link and run with stdio, math and -l" $?

CWRAP_KEEP=1 "$cwrap" -c hello.c -o hello.o 2> keep.log
kept=$(sed -n 's/^cwrap: kept //p' keep.log)
cmp -s "$kept/0-hello.i" "$kept/0-hello.t.i"
check "translator output is byte-identical to its input" $?
rm -rf "$kept"

"$cwrap" -c bad.c -o bad.o 2> bad.log
[ $? -ne 0 ] && grep -q '^bad.c:5:' bad.log
check "gcc errors point at the original file and line" $?

"$cwrap" a.c b.c -o two && ./two
check "several .c files in one command" $?

"$cwrap" -c a.c && "$cwrap" -c b.c -o b2.o && "$cwrap" a.o b2.o -o two2 && ./two2
check "separate compile, then link objects" $?

"$cwrap" -MD -c dep.c -o dep.o
grep -q '^dep.o:' dep.d && grep -q 'dep.h' dep.d
check "-MD writes a dependency file for the object" $?

[ "$("$cwrap" -E -DN=9 -P - <<'EOF'
N
EOF
)" = 9 ]
check "-E passes straight through to gcc" $?

"$cwrap" -c stray.c 2> stray.log
[ $? -ne 0 ] && grep -q '^stray.c:1: warning: not translated' stray.log &&
  grep -q "^stray.c:1:.*error: stray" stray.log
check "a file that does not parse is left to gcc, which reports it" $?

"$cwrap" -c errors.c -o errors.o 2> errors.log
[ $? -ne 0 ] &&
  [ "$(grep -o '^errors.c:[0-9]*: warning' errors.log | cut -d: -f2 | tr '\n' ' ')" = "2 4 6 9 11 " ] &&
  [ "$(grep -c 'not translated' errors.log)" = 1 ] && grep -q '^errors.c:9:19: error: stray' errors.log
check "cwrap lists each syntax error, and none of its own making, before gcc" $?

"$cwrap" -c errors_cap.c -o errors_cap.o 2> errors_cap.log
[ "$(grep -c '^errors_cap.c:[0-9]*: warning' errors_cap.log)" = 20 ] &&
  grep -q '^errors_cap.c:20: note: 20 errors, not looking for more' errors_cap.log
check "cwrap stops listing syntax errors after 20" $?

"$cwrap" -c foreach_twice.c -o foreach_twice.o 2> foreach_twice.log
[ $? -ne 0 ] && [ "$(grep -c ': error: foreach loop' foreach_twice.log)" = 2 ] &&
  [ ! -e foreach_twice.o ]
check "every foreach that clashes with a declared foreach is reported" $?

"$cwrap" -w kr.c -o kr 2> kr.log && ./kr
check "code gcc accepts still builds when the grammar does not cover it" $?

"$cwrap" -c undeclared.c 2> undeclared.log
[ $? -ne 0 ] && grep -q 'unknown type name' undeclared.log
check "gcc reports an undeclared type" $?

"$cwrap" -c conflict.c 2> conflict.log
[ $? -ne 0 ] && grep -q 'conflicting types' conflict.log &&
  ! grep -q 'not translated' conflict.log
check "the scope table accepts a conflict, and gcc reports it" $?

"$scope_test"
check "scope table unit checks" $?

cp "$tests/grammar.c" grammar.c
CWRAP_KEEP=1 "$cwrap" -std=gnu11 -w -c grammar.c -o grammar.o 2> keep.log
rc=$?
kept=$(sed -n 's/^cwrap: kept //p' keep.log)
[ "$rc" = 0 ] && cmp -s "$kept/0-grammar.i" "$kept/0-grammar.t.i"
check "the full grammar accepts tests/grammar.c, byte for byte" $?
rm -rf "$kept"

CWRAP_DUMP_AST=1 "$cwrap" -c loc.c -o loc.o 2> loc.log
grep -q "'in_header' loc.h:1" loc.log && grep -q "'in_source' loc.c:3" loc.log
check "each leaf records its file and line" $?

CWRAP_DUMP_AST=1 "$cwrap" -c prec.c -o prec.o 2> prec.log
[ "$(grep -nE 'ADDITIVE|MULTIPLICATIVE' prec.log | sed 's/.*[0-9]: *//')" = "ADDITIVE
MULTIPLICATIVE" ]
check "operators nest by precedence" $?

"$cwrap" -w -c hide.c -o hide.o
check "an inner variable hides a typedef" $?

CWRAP_DUMP_AST=1 "$cwrap" -translate foreach.c -o foreach.out 2> foreach.log
[ $? -eq 0 ] && grep -q FOR_EACH_STATEMENT foreach.log &&
  ! grep -q 'not translated' foreach.log &&
  grep -q '^  for (int __i_v = 0; __i_v < V.n; ++__i_v$' foreach.out &&
  grep -q '^            ) { __auto_type v = V.a\[__i_v\]; { g (v); } }$' foreach.out
check "foreach (v in V) becomes a for loop, on as many lines" $?

"$cwrap" -Wall foreach_run.c -o foreach_run 2> foreach_run.log &&
  [ "$(./foreach_run)" = 150 ] &&
  grep -q '^foreach_run.c:17:7: warning: unused variable' foreach_run.log
check "translated foreach loops run, nested too, and later lines keep their numbers" $?

CWRAP_DUMP_AST=1 "$cwrap" -translate foreach_gap.c -o foreach_gap.out 2> foreach_gap.log
[ $? -eq 0 ] && grep -q FOR_EACH_STATEMENT foreach_gap.log
check "a foreach header may span linemarkers" $?

CWRAP_DUMP_AST=1 "$cwrap" -translate foreach_scope.c -o foreach_scope.out 2> foreach_scope.log
[ $? -eq 0 ] && [ "$(grep -c MULTIPLICATIVE foreach_scope.log)" = 1 ]
check "the foreach variable hides a typedef until the loop ends" $?

"$cwrap" -translate foreach_names.c -o foreach_names.out 2> foreach_names.log
[ $? -eq 0 ] && ! grep -q 'not translated' foreach_names.log
check "in may also be the names in a foreach" $?

CWRAP_DUMP_AST=1 "$cwrap" -translate foreach_types.c -o foreach_types.out 2> foreach_types.log
[ $? -eq 0 ] && [ "$(grep -c FOR_EACH_STATEMENT foreach_types.log)" = 1 ]
check "a foreach works while in is a type name" $?

"$cwrap" -translate foreach_scoped.c -o foreach_scoped.out 2> foreach_scoped.log
[ $? -eq 0 ] && ! grep -q 'not translated' foreach_scoped.log
check "a member or parameter named foreach leaves the keyword alone" $?

"$cwrap" -c foreach_declared.c -o foreach_declared.o 2> foreach_declared.log
[ $? -ne 0 ] && grep -q '^foreach_declared.c:4: error: ' foreach_declared.log &&
  grep -q '^foreach_declared.c:1: note: ' foreach_declared.log &&
  ! grep -q 'not translated' foreach_declared.log &&
  [ "$(wc -l < foreach_declared.log)" = 2 ] && [ ! -e foreach_declared.o ]
check "foreach (v in V) stops in cwrap if foreach is declared, and says where" $?

"$cwrap" -translate foreach_declared.c -o foreach_declared.out 2> foreach_declared.log
[ $? -eq 1 ] && grep -q ': error: ' foreach_declared.log && [ ! -e foreach_declared.out ]
check "-translate stops there too, and writes nothing" $?

"$cwrap" -c foreach_header.c -o foreach_header.o 2> foreach_header.log
[ $? -ne 0 ] && grep -q '^foreach_taken.h:1: note: ' foreach_header.log
check "the note points into the header that declares foreach" $?

"$cwrap" -c foreach_typedef.c -o foreach_typedef.o 2> foreach_typedef.log
[ $? -ne 0 ] && grep -q '^foreach_typedef.c:2: warning: not translated' foreach_typedef.log
check "a foreach while foreach is a typedef is left to gcc" $?

"$cwrap" -c foreach_gap2.c -o foreach_gap2.o 2> foreach_gap2.log
[ $? -ne 0 ] && grep -q '^foreach_gap2.c:4: error: ' foreach_gap2.log
check "a foreach split by linemarkers still stops if foreach is declared" $?

"$cwrap" -translate foreach_closed.c -o foreach_closed.out 2> foreach_closed.log
[ $? -eq 0 ] && [ ! -s foreach_closed.log ]
check "a foreach whose scope has closed does not stop a foreach" $?

"$cwrap" -w -c foreach_lookalike.c -o foreach_lookalike.o
check "a declaration that looks like a foreach still compiles" $?

for f in foreach_bad*.c; do
  "$cwrap" -translate "$f" > /dev/null 2> foreach_bad.log
  [ $? -eq 1 ] && grep -q 'not translated' foreach_bad.log
  check "$(sed 's/.*{ //; s/ ; }//' "$f") is rejected" $?
done

"$cwrap" -translate scalar.c -o scalar.out 2> scalar.log
[ "$(grep -o '__attribute__ ((section ("scalars"), used)) = &[a-z]*;' scalar.out |
     sed 's/.*&//' | tr '\n' ' ')" = "a; b; d; s; last; " ] &&
  grep -q '^scalar a; static scalar \*const __cwrap_scalar_slot_0 ' scalar.out &&
  grep -q '^extern scalar e;$' scalar.out && grep -q '^typedef scalar t;$' scalar.out &&
  grep -q '^int g (void) { scalar local = 0; return local; }$' scalar.out
check "each global scalar gets a slot in the section scalars" $?

"$cwrap" -w scalar.c -o scalar 2> scalar.log; ./scalar
[ $? = 10 ]
check "the constructor numbers the global scalars 0 to 4" $?

"$cwrap" -translate scalar_local.c -o scalar_local.out &&
  grep -q '^  scalar a __attribute__ ((cleanup (__cwrap_scalar_free))) = scalar_new ();$' scalar_local.out &&
  grep -q '^  static scalar s, \*const __cwrap_local_slot_[0-9]* __attribute__ ((section ("scalars"), used)) = &s;$' scalar_local.out &&
  grep -q '^  scalar kept = scalar_new ();$' scalar_local.out
check "local scalars get an id from scalar_new (), and a cleanup that frees it" $?

"$cwrap" -Wall scalar_local.c -o scalar_local 2> scalar_local.log &&
  [ "$(./scalar_local | tr '\n' ' ')" = "16 in use 3 kept 2 gs 1 globals 2 " ] &&
  grep -q '^scalar_local.c:22:7: warning: unused variable' scalar_local.log &&
  [ "$(grep -c 'warning:' scalar_local.log)" = 1 ]
check "local scalars are freed however their scope ends, in reverse order" $?

"$cwrap" -c scalar_other.c -o scalar_other.o && "$cwrap" -c scalar_main.c -o scalar_main.o &&
  "$cwrap" scalar_main.o scalar_other.o -o scalars && [ "$(./scalars)" = "6 3f" ]
check "scalars of separately compiled files are numbered 0 to n-1 together" $?

"$cwrap" -O2 scalar_main.c scalar_other.c -o scalars2 && [ "$(./scalars2)" = "6 3f" ]
check "and at -O2 in a single command" $?

"$cwrap" -c assign.c -o assign.o 2> assign.log
rc=$?
[ "$rc" = 0 ] && [ -e assign.o ] &&
  [ "$(grep -o '^assign.c:[0-9]*: warning' assign.log | cut -d: -f2 | tr '\n' ' ')" = "7 8 9 10 16 3 " ] &&
  grep -q "^assign.c:9: warning: '++' changes scalar 'a'" assign.log &&
  grep -q "^assign.c:2: note: 'a' is declared here" assign.log &&
  grep -q "^assign.c:5: note: 's' is declared here" assign.log &&
  grep -q "^assign.c:3: warning: scalar 'c' has a value of its own" assign.log
check "changing a scalar is a warning, and the file still compiles" $?

"$cwrap" -Wall scalar_runtime.c -o scalar_runtime 2> scalar_runtime.log &&
  [ ! -s scalar_runtime.log ] &&
  [ "$(./scalar_runtime | tr '\n' ' ')" = "0 1 0 1 0 0,1,2, 0,2, 64 1 64 63 " ]
check "scalar.h comes with every file: ids reused lowest first, all 64 usable, foreach_scalar" $?

./scalar_runtime more > /dev/null 2> scalar_runtime.log
[ $? -ne 0 ] && grep -q 'scalar_new: all 64 scalars are in use' scalar_runtime.log
check "scalar_new () stops the program when no id is left" $?

"$cwrap" -c scalar_runtime.c -o scalar_runtime.o && "$cwrap" scalar_runtime.o -o scalar_runtime2 &&
  ./scalar_runtime2 > /dev/null
check "linking objects alone also brings in the runtime" $?

"$cwrap" -translate -DVALUE=42 hello.c -o hello.tr.c
gcc -w hello.tr.c -o hello.tr -lm && [ "$(./hello.tr)" = '5.0 42 a"b' ]
check "-translate -o writes C that gcc can build" $?

"$cwrap" -translate -DVALUE=42 hello.c > hello.out
cmp -s hello.tr.c hello.out
check "-translate writes to standard output without -o" $?

rm -f flags.o flags
"$cwrap" -translate -c -O2 -Wall -lm missing.o flags.c > flags.out &&
  [ -s flags.out ] && [ ! -e flags.o ] && [ ! -e flags ]
check "-translate disregards compile and link flags" $?

"$cwrap" -translate a.c b.c > both.out
grep -q 'int twice(int x)' both.out && grep -q 'int main(void)' both.out
check "-translate takes several files on standard output" $?

"$cwrap" -translate a.c b.c -o both.c 2> both.log
[ $? -ne 0 ] && grep -q 'multiple files' both.log
check "-translate refuses -o with several files" $?

"$cwrap" -translate -x c a.c 2> x.log
[ $? -ne 0 ] && grep -q 'cannot be used with -x' x.log
check "-translate refuses flags it cannot honour" $?

"$cwrap" -translate stray.c > stray.out 2> stray2.log
[ $? -eq 1 ] && grep -q 'not translated' stray2.log && grep -q '@' stray.out
check "-translate fails on code that does not parse, and still writes it" $?

cp keep.c keep.orig
ln -s keep.c keep.link
for args in '-translate keep.c -o keep.c' '-translate keep.c -o ./keep.c' \
            '-translate keep.c -o keep.link' '-c keep.c -o keep.c' \
            '-S keep.c -o keep.c' 'keep.c -o keep.c'; do
  "$cwrap" $args 2> keep.log
  [ $? -ne 0 ] && grep -q 'is the same as output file' keep.log &&
    cmp -s keep.c keep.orig
  check "cwrap $args leaves the source alone" $?
  cp keep.orig keep.c
done

[ "$fails" = 0 ] && echo "all passed"
exit "$fails"
