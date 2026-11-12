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

"$cwrap" -include cwrap/cache.h -c foreach_twice.c -o foreach_twice.o 2> foreach_twice.log
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

CWRAP_DUMP_AST=1 "$cwrap" -include cwrap/cache.h -translate foreach.c -o foreach.out 2> foreach.log
[ $? -eq 0 ] && grep -q FOR_EACH_STATEMENT foreach.log &&
  ! grep -q 'not translated' foreach.log &&
  grep -q '^  for (int __i_v = 0; __i_v < V.n; ++__i_v$' foreach.out &&
  grep -q '^            ) { __auto_type v = V.a\[__i_v\]; { g (v); } }$' foreach.out
check "foreach (v in V) becomes a for loop, on as many lines" $?

"$cwrap" -include cwrap/cache.h -Wall foreach_run.c -o foreach_run 2> foreach_run.log &&
  [ "$(./foreach_run)" = 150 ] &&
  grep -q '^foreach_run.c:17:7: warning: unused variable' foreach_run.log
check "translated foreach loops run, nested too, and later lines keep their numbers" $?

CWRAP_DUMP_AST=1 "$cwrap" -include cwrap/cache.h -translate foreach_gap.c -o foreach_gap.out 2> foreach_gap.log
[ $? -eq 0 ] && grep -q FOR_EACH_STATEMENT foreach_gap.log
check "a foreach header may span linemarkers" $?

CWRAP_DUMP_AST=1 "$cwrap" -include cwrap/cache.h -translate foreach_scope.c -o foreach_scope.out 2> foreach_scope.log
[ $? -eq 0 ] && [ "$(grep -A1 MULTIPLICATIVE foreach_scope.log | grep -c "foreach_scope.c:")" = 1 ]
check "the foreach variable hides a typedef until the loop ends" $?

"$cwrap" -include cwrap/cache.h -translate foreach_names.c -o foreach_names.out 2> foreach_names.log
[ $? -eq 0 ] && ! grep -q 'not translated' foreach_names.log
check "in may also be the names in a foreach" $?

CWRAP_DUMP_AST=1 "$cwrap" -include cwrap/cache.h -translate foreach_types.c -o foreach_types.out 2> foreach_types.log
[ $? -eq 0 ] && [ "$(grep -c FOR_EACH_STATEMENT foreach_types.log)" = 1 ]
check "a foreach works while in is a type name" $?

"$cwrap" -include cwrap/cache.h -translate foreach_scoped.c -o foreach_scoped.out 2> foreach_scoped.log
[ $? -eq 0 ] && ! grep -q 'not translated' foreach_scoped.log
check "a member or parameter named foreach leaves the keyword alone" $?

"$cwrap" -include cwrap/cache.h -c foreach_declared.c -o foreach_declared.o 2> foreach_declared.log
[ $? -ne 0 ] && grep -q '^foreach_declared.c:4: error: ' foreach_declared.log &&
  grep -q '^foreach_declared.c:1: note: ' foreach_declared.log &&
  ! grep -q 'not translated' foreach_declared.log &&
  [ "$(wc -l < foreach_declared.log)" = 2 ] && [ ! -e foreach_declared.o ]
check "foreach (v in V) stops in cwrap if foreach is declared, and says where" $?

"$cwrap" -include cwrap/cache.h -translate foreach_declared.c -o foreach_declared.out 2> foreach_declared.log
[ $? -eq 1 ] && grep -q ': error: ' foreach_declared.log && [ ! -e foreach_declared.out ]
check "-translate stops there too, and writes nothing" $?

"$cwrap" -include cwrap/cache.h -c foreach_header.c -o foreach_header.o 2> foreach_header.log
[ $? -ne 0 ] && grep -q '^foreach_taken.h:1: note: ' foreach_header.log
check "the note points into the header that declares foreach" $?

"$cwrap" -include cwrap/cache.h -c foreach_typedef.c -o foreach_typedef.o 2> foreach_typedef.log
[ $? -ne 0 ] && grep -q '^foreach_typedef.c:2: warning: not translated' foreach_typedef.log
check "a foreach while foreach is a typedef is left to gcc" $?

"$cwrap" -include cwrap/cache.h -c foreach_gap2.c -o foreach_gap2.o 2> foreach_gap2.log
[ $? -ne 0 ] && grep -q '^foreach_gap2.c:4: error: ' foreach_gap2.log
check "a foreach split by linemarkers still stops if foreach is declared" $?

"$cwrap" -include cwrap/cache.h -translate foreach_closed.c -o foreach_closed.out 2> foreach_closed.log
[ $? -eq 0 ] && [ ! -s foreach_closed.log ]
check "a foreach whose scope has closed does not stop a foreach" $?

"$cwrap" -include cwrap/cache.h -w -c foreach_lookalike.c -o foreach_lookalike.o
check "a declaration that looks like a foreach still compiles" $?

for f in foreach_bad*.c; do
  "$cwrap" -include cwrap/cache.h -translate "$f" > /dev/null 2> foreach_bad.log
  [ $? -eq 1 ] && grep -q 'not translated' foreach_bad.log
  check "$(sed 's/.*{ //; s/ ; }//' "$f") is rejected" $?
done

"$cwrap" -Wall -Wextra cache_run.c -o cache_run 2> cache_run.log &&
  [ ! -s cache_run.log ] && [ "$(./cache_run)" = "7410 1 20 2 hi 3" ]
check "caches of any type grow, update once when dirty, and cache stays a name" $?

"$cwrap" -translate cache_run.c -o cache_run.out &&
  grep -q '^struct cache_int { __typeof__ (int) \*a; int n, cap, is_dirty; void (\*update) (struct cache_int \*self); }; struct cache_int squares = ' cache_run.out &&
  grep -q '^struct cache_point { .* }; struct cache_char_p { __typeof__ (char \*) \*a; .* }; int main (void)$' cache_run.out &&
  grep -q '^  struct cache_point pts __attribute__ ((cleanup (__cwrap_cache_free))) = { 0 };$' cache_run.out &&
  grep -q '^  struct cache_char_p words __attribute__ ((cleanup (__cwrap_cache_free))) = { .update = ' cache_run.out &&
  grep -q '^  for (int __i_s = squares.update && squares.is_dirty ? (squares.is_dirty = 0, squares.update (&squares), 0) : 0; __i_s < squares.n; ++__i_s) { __auto_type s = squares.a\[__i_s\];$' cache_run.out
check "cache (T) becomes a struct defined on the line of its first use" $?

"$cwrap" -c cache_errors.c -o cache_errors.o 2> cache_errors.log
[ $? -ne 0 ] && [ ! -e cache_errors.o ] &&
  [ "$(grep -o '^cache_errors.c:[0-9]*: error' cache_errors.log | cut -d: -f2 | tr '\n' ' ')" = "4 8 9 11 12 15 " ] &&
  grep -q "^cache_errors.c:8: error: '=' copies cache 'w'" cache_errors.log &&
  grep -q "^cache_errors.c:7: note: 'w' is declared here" cache_errors.log &&
  grep -q "^cache_errors.c:15: error: cache (int_p) and cache (int \*) would both be struct cache_int_p" cache_errors.log &&
  ! grep -q 'not translated' cache_errors.log
check "copying a cache, or a typedef of one, stops in cwrap" $?

esc=$(printf '\033')
"$cwrap" -fdiagnostics-color=always -c cache_errors.c -o cache_errors.o 2> color.log
grep -q "^${esc}\[01m${esc}\[Kcache_errors.c:4:${esc}\[m${esc}\[K ${esc}\[01;31m${esc}\[Kerror:${esc}\[m${esc}\[K a typedef" color.log &&
  grep -q "${esc}\[01;36m${esc}\[Knote:" color.log &&
  ! grep -q "$esc" cache_errors.log &&
  "$cwrap" -fdiagnostics-color=always -fno-diagnostics-color -c cache_errors.c -o cache_errors.o 2> color.log
[ $? -ne 0 ] && ! grep -q "$esc" color.log
check "diagnostics are colored as gcc colors them, and only when asked or on a terminal" $?

"$cwrap" -Wall cache_warn.c -o cache_warn 2> cache_warn.log && ./cache_warn &&
  [ "$(grep -c "^cache_warn.c:7: warning: cache 'v' is passed by value" cache_warn.log)" = 2 ]
check "a cache passed by value is a warning" $?

"$cwrap" -Wall cache_off.c -o cache_off 2> cache_off.log && [ "$(./cache_off)" = 15 ] &&
  [ ! -s cache_off.log ]
check "without <cwrap/cache.h>, cache and foreach are ordinary names" $?

"$cwrap" -dimension=2 -Wall tensor.c -o tensor2 2> tensor2.log &&
  [ "$(./tensor2 | tr '\n' ' ')" = "0 1 2 3 4 5 6 9 9 24 9 1406 " ] &&
  [ "$(grep -o '^tensor.c:[0-9]*: warning' tensor2.log | cut -d: -f2 | tr '\n' ' ')" = "19 20 21 " ] &&
  grep -q "^tensor.c:19: warning: '=' changes scalar 't.x.y'" tensor2.log &&
  grep -q "^tensor.c:21: warning: '=' changes vector 't.y'" tensor2.log
check "-dimension=2: a vector has 2 scalars and a tensor 4, numbered or freed each" $?

"$cwrap" -dimension=2 -translate tensor.c -o tensor.out 2> /dev/null &&
  grep -q '^tensor t; static const struct __cwrap_scalar_slot __cwrap_scalar_slot_3 __attribute__ ((section ("scalars"), used)) = { &t.x.x, "t.x.x" }; .*_slot_4 .* = { &t.x.y, "t.x.y" }; .*_slot_5 .* = { &t.y.x, "t.y.x" }; .*_slot_6 .* = { &t.y.y, "t.y.y" };$' tensor.out &&
  grep -q '^  tensor w __attribute__ ((cleanup (__cwrap_tensor_free))) = __cwrap_tensor_new_named ((const char \*const \[\]) { "w.x.x", "w.x.y", "w.y.x", "w.y.y" });$' tensor.out &&
  grep -q '^  static vector kept; static const struct __cwrap_scalar_slot __cwrap_local_slot_[0-9]* .* = { &kept.x, "sum.kept.x" }; static const struct __cwrap_scalar_slot __cwrap_local_slot_[0-9]* .* = { &kept.y, "sum.kept.y" };$' tensor.out
check "a tensor gets a slot for each scalar, t.x.x to t.y.y, on its line, each with its name" $?

"$cwrap" -dimension=3 -Wall tensor.c -o tensor3 2> /dev/null &&
  [ "$(./tensor3 | tr '\n' ' ')" = "0 1 2 4 5 7 8 16 16 40 16 2312 " ]
check "-dimension=3: a vector has 3 scalars and a tensor 9" $?

"$cwrap" -Wall own_scalar.c -o own_scalar 2> own_scalar.log && [ "$(./own_scalar)" = 1.5 ] &&
  [ ! -s own_scalar.log ]
check "without -dimension, scalar and vector are the program's own names" $?

"$cwrap" -dimension=4 -c own_scalar.c 2> dimension.log
[ $? -ne 0 ] && grep -q "error: -dimension must be 1, 2 or 3, not '4'" dimension.log
check "-dimension takes only 1, 2 or 3" $?

"$cwrap" -dimension=1 -translate scalar.c -o scalar.out 2> scalar.log
[ "$(grep -o '__attribute__ ((section ("scalars"), used)) = { &[a-z]*, "[a-z]*" };' scalar.out |
     sed 's/.*&//' | tr '\n' ' ')" = 'a, "a" }; b, "b" }; d, "d" }; s, "s" }; last, "last" }; ' ] &&
  grep -q '^scalar a; static const struct __cwrap_scalar_slot __cwrap_scalar_slot_0 ' scalar.out &&
  grep -q '^extern scalar e;$' scalar.out && grep -q '^typedef scalar t;$' scalar.out &&
  grep -q '^int g (void) { scalar local = 0; return local; }$' scalar.out
check "each global scalar gets a slot in the section scalars, with its name" $?

"$cwrap" -dimension=1 -w scalar.c -o scalar 2> scalar.log; ./scalar
[ $? = 10 ]
check "the constructor numbers the global scalars 0 to 4" $?

"$cwrap" -dimension=1 -translate scalar_local.c -o scalar_local.out &&
  grep -q '^  scalar a __attribute__ ((cleanup (__cwrap_scalar_free))) = __cwrap_scalar_new_named ("a");$' scalar_local.out &&
  grep -q '^  static scalar s; static const struct __cwrap_scalar_slot __cwrap_local_slot_[0-9]* __attribute__ ((section ("scalars"), used)) = { &s, "counted.s" };$' scalar_local.out &&
  grep -q '^  scalar kept = scalar_new ();$' scalar_local.out
check "local scalars get a named id from __cwrap_scalar_new_named (), and a cleanup that frees it" $?

"$cwrap" -dimension=1 -Wall scalar_local.c -o scalar_local 2> scalar_local.log &&
  [ "$(./scalar_local | tr '\n' ' ')" = "16 in use 3 kept 2 gs 1 globals 2 " ] &&
  grep -q '^scalar_local.c:22:7: warning: unused variable' scalar_local.log &&
  [ "$(grep -c 'warning:' scalar_local.log)" = 1 ]
check "local scalars are freed however their scope ends, in reverse order" $?

"$cwrap" -dimension=1 -c scalar_other.c -o scalar_other.o && "$cwrap" -dimension=1 -c scalar_main.c -o scalar_main.o &&
  "$cwrap" -dimension=1 scalar_main.o scalar_other.o -o scalars && [ "$(./scalars)" = "6 3f" ]
check "scalars of separately compiled files are numbered 0 to n-1 together" $?

"$cwrap" -dimension=1 -O2 scalar_main.c scalar_other.c -o scalars2 && [ "$(./scalars2)" = "6 3f" ]
check "and at -O2 in a single command" $?

"$cwrap" -dimension=1 -c assign.c -o assign.o 2> assign.log
rc=$?
[ "$rc" = 0 ] && [ -e assign.o ] &&
  [ "$(grep -o '^assign.c:[0-9]*: warning' assign.log | cut -d: -f2 | tr '\n' ' ')" = "7 8 9 10 16 3 " ] &&
  grep -q "^assign.c:9: warning: '++' changes scalar 'a'" assign.log &&
  grep -q "^assign.c:2: note: 'a' is declared here" assign.log &&
  grep -q "^assign.c:5: note: 's' is declared here" assign.log &&
  grep -q "^assign.c:3: warning: scalar 'c' has a value of its own" assign.log
check "changing a scalar is a warning, and the file still compiles" $?

"$cwrap" -dimension=1 -Wall scalar_runtime.c -o scalar_runtime 2> scalar_runtime.log &&
  [ ! -s scalar_runtime.log ] &&
  [ "$(./scalar_runtime | tr '\n' ' ')" = "0 1 0 1 0 64 1 " ]
check "scalar.h comes with every file: ids reused lowest first, all 64 usable" $?

./scalar_runtime more > /dev/null 2> scalar_runtime.log
[ $? -ne 0 ] && grep -q 'scalar_new: all 64 scalars are in use' scalar_runtime.log
check "scalar_new () stops the program when no id is left" $?

"$cwrap" -dimension=1 -c scalar_runtime.c -o scalar_runtime.o && "$cwrap" -dimension=1 scalar_runtime.o -o scalar_runtime2 &&
  ./scalar_runtime2 > /dev/null
check "linking objects alone also brings in the runtime" $?

"$cwrap" -dimension=1 -Wall loop_run.c loop_other.c -o loop_run 2> loop_run.log &&
  [ "$(./loop_run)" = "46 144127 10 104 -1 7 9 1 0 7" ] && [ ! -s loop_run.log ]
check "@loop: typed arguments, hidden locals, given names, nesting, break, scalars, caches" $?

"$cwrap" -c loop_bad.c -o loop_bad.o 2> loop_bad.log
[ $? -ne 0 ] && [ ! -e loop_bad.o ] &&
  [ "$(grep -o '^loop_bad.c:[0-9]*: error' loop_bad.log | cut -d: -f2 | tr '\n' ' ')" = "1 2 3 4 5 6 7 9 10 11 11 " ] &&
  grep -q "^loop_bad.c:7: error: 'z' is not declared by the loop where {...} is" loop_bad.log &&
  grep -q "^loop_bad.c:8: note: 'ok' is defined here" loop_bad.log &&
  grep -q "^loop_bad.c:11: error: loop 'ok' takes 1 argument, not 0" loop_bad.log
check "@loop: each misuse is an error, and the file never reaches gcc" $?

"$cwrap" -c loop_lines.c -o loop_lines.o 2> loop_lines.log
[ $? -ne 0 ] && grep -q "^loop_lines.c:8:[0-9]*: error: .*no member named .nope." loop_lines.log &&
  grep -q "^loop_lines.c:10:[0-9]*: error: .*no member named .nope2." loop_lines.log &&
  grep -q "^loop_lines.c:11:[0-9]*: warning: initialization of .const Bag \*. from incompatible pointer type" loop_lines.log
check "@loop: errors in the block, in an argument and after the use keep their lines" $?

"$cwrap" -dimension=2 -Wall -Wextra stencil_run.c -o stencil_run 2> stencil_run.log &&
  [ "$(./stencil_run)" = "48 864 24 18 9 336 99" ] && [ ! -s stencil_run.log ]
check "foreach (): every cell of every leaf, s[i,j] offsets, a vector, the body's own names, x, y, Delta" $?

"$cwrap" -dimension=2 -translate stencil_run.c -o stencil_run.out &&
  [ "$(grep -c '[^t] Real \* const __cwrap_block_2_0 = block_data ( __cwrap_foreach_leaf_sblocks_[0-9]* , ( s ) )' stencil_run.out)" = 1 ] &&
  [ "$(grep -c 'const Real \* const __cwrap_block_2_1 = block_data ( __cwrap_foreach_leaf_sblocks_[0-9]* , ( t ) )' stencil_run.out)" = 1 ] &&
  [ "$(grep -c '__cwrap_block_2_[2-9]' stencil_run.out)" = 0 ]
check "foreach (): a block for each scalar used, once, const where only read" $?

"$cwrap" -dimension=2 -c stencil_lines.c -o stencil_lines.o 2> stencil_lines.log
[ $? -ne 0 ] && grep -q "^stencil_lines.c:7:[0-9]*: error: .*nope1" stencil_lines.log &&
  grep -q "^stencil_lines.c:9:[0-9]*: error: .*nope2" stencil_lines.log
check "foreach (): errors in the body and after it keep their lines" $?

"$cwrap" -dimension=2 -c stencil_bad.c -o stencil_bad.o 2> stencil_bad.log
[ $? -ne 0 ] && [ ! -e stencil_bad.o ] &&
  [ "$(grep -o '^stencil_bad.c:[0-9]*: error' stencil_bad.log | cut -d: -f2 | tr '\n' ' ')" = "6 7 8 9 13 17 17 18 " ] &&
  grep -q "^stencil_bad.c:18: error: 's\[1\]' writes a cell other than its own" stencil_bad.log &&
  grep -q "^stencil_bad.c:17: error: 'x' in foreach () is the cell's coordinate, which hides scalar 'x'" stencil_bad.log &&
  grep -q "^stencil_bad.c:14: note: 'x' is declared here" stencil_bad.log &&
  grep -q "^stencil_bad.c:17: error: 'y' in foreach () is the cell's coordinate, which hides vector 'y'" stencil_bad.log &&
  grep -q "^stencil_bad.c:13: error: 'Delta\[...\]': in foreach (), x, y, z and Delta are the cell's" stencil_bad.log &&
  grep -q "^stencil_bad.c:7: error: 'a\[\]', but 'a' is not a scalar" stencil_bad.log &&
  grep -q "^stencil_bad.c:8: error: 's\[1,2,3\]' has 3 offsets, in 2 dimensions" stencil_bad.log
check "foreach (): s[] outside it, x[] of a non-scalar, too many offsets, nesting, a scalar x hidden" $?

"$cwrap" -dimension=2 -c stencil_noloops.c -o stencil_noloops.o 2> stencil_noloops.log
[ $? -ne 0 ] && grep -q "^stencil_noloops.c:4: error: foreach () is not declared: include <block-loops.h>" stencil_noloops.log
check "foreach (): without its @foreach, an error" $?

"$cwrap" -dimension=2 -Wall -Wextra stencil_neighbor.c -o stencil_neighbor 2> stencil_neighbor.log &&
  [ "$(./stencil_neighbor)" = "300 0 | 3 " ] && [ ! -s stencil_neighbor.log ]
check "foreach_neighbor (r): the box around the point, offsets in it, x and Delta the neighbour's, ghosts updated first" $?

"$cwrap" -dimension=2 -translate stencil_neighbor.c -o stencil_neighbor.out &&
  [ "$(grep -c '( ( Real ) __cwrap_block_[0-9]*_[0-9]* \[ [^;]*__cwrap_di_' stencil_neighbor.out)" = 3 ]
check "foreach_neighbor (): its reads are values, not lvalues, for gcc too" $?

"$cwrap" -dimension=2 -c stencil_neighbor_bad.c -o stencil_neighbor_bad.o 2> stencil_neighbor_bad.log
[ $? -ne 0 ] && [ ! -e stencil_neighbor_bad.o ] &&
  grep -q "^stencil_neighbor_bad.c:5: error: foreach_neighbor () outside a stencil loop" stencil_neighbor_bad.log &&
  grep -q "^stencil_neighbor_bad.c:6: error: foreach_neighbor () inside foreach_neighbor ()" stencil_neighbor_bad.log &&
  grep -q "^stencil_neighbor_bad.c:7: error: 's\[\]' in foreach_neighbor () writes a neighbour" stencil_neighbor_bad.log
check "foreach_neighbor (): outside a stencil loop, chained, or writing, an error" $?

"$cwrap" -dimension=2 -Wall -c stencil_reach.c -o stencil_reach.o 2> stencil_reach.log
[ $? -ne 0 ] && [ ! -e stencil_reach.o ] &&
  [ "$(grep -o '^stencil_reach.c:[0-9]*:[0-9]*: error' stencil_reach.log | cut -d: -f2 | tr '\n' ' ')" = "5 6 7 8 10 " ] &&
  grep -q "^stencil_reach.c:10:[0-9]*: error: static assertion failed: \"coarse (t, 2): an offset is not a constant within BAMR_COARSE_REACH" stencil_reach.log &&
  grep -q "^stencil_reach.c:5:[0-9]*: error: static assertion failed: \"t\[3\]: an offset is not a constant within BAMR_GHOSTS" stencil_reach.log &&
  grep -q "^stencil_reach.c:6:[0-9]*: error: expression in static assertion is not constant" stencil_reach.log &&
  grep -q "^stencil_reach.c:7:[0-9]*: error: static assertion failed: \"foreach_neighbor (): its radius is not a constant" stencil_reach.log &&
  grep -q "^stencil_reach.c:8:[0-9]*: error: .*t\[0,1\].*adding the radius of foreach_neighbor" stencil_reach.log
check "offsets and radii past BAMR_GHOSTS, or not constant, are compile errors at their lines" $?

"$cwrap" -dimension=2 -Wall -Wextra stencil_levels.c -o stencil_levels 2> stencil_levels.log &&
  [ "$(./stencil_levels)" = "0 12 | " ] && [ ! -s stencil_levels.log ]
check "two levels: coarse (s, i, j) and child.x against the same by hand, foreach_child () a mean, no ghost tracked" $?

"$cwrap" -dimension=2 -c stencil_levels_bad.c -o stencil_levels_bad.o 2> stencil_levels_bad.log
[ $? -ne 0 ] && [ ! -e stencil_levels_bad.o ] &&
  [ "$(grep -o '^stencil_levels_bad.c:[0-9]*: error' stencil_levels_bad.log | cut -d: -f2 | tr '\n' ' ')" = "5 6 7 8 9 10 " ] &&
  grep -q "^stencil_levels_bad.c:5: error: 't\[1\]' in foreach_up (), a loop over two levels: no offset on the point's own level" stencil_levels_bad.log &&
  grep -q "^stencil_levels_bad.c:6: error: 'coarse (s)' writes the coarse level" stencil_levels_bad.log &&
  grep -q "^stencil_levels_bad.c:7: error: foreach_child () in foreach_up (), whose points have no children" stencil_levels_bad.log &&
  grep -q "^stencil_levels_bad.c:8: error: foreach_neighbor () in foreach_up (), a loop over two levels" stencil_levels_bad.log &&
  grep -q "^stencil_levels_bad.c:9: error: 'coarse (t)' in foreach (), whose points have no coarse cell" stencil_levels_bad.log &&
  grep -q "^stencil_levels_bad.c:10: error: 't\[\]' in foreach_child () writes a child" stencil_levels_bad.log
check "two levels: an own-level offset, writing coarse (), foreach_child () without children, foreach_neighbor (), coarse () in foreach (), writing a child" $?

"$cwrap" -dimension=2 -Wall -Wextra stencil_params.c -o stencil_params 2> stencil_params.log &&
  [ "$(./stencil_params)" = "6 1" ] && [ ! -s stencil_params.log ]
check "@foreach with a parameter: given to its block loop, its argument evaluated once" $?

"$cwrap" -Wall -Wextra ops_run.c -o ops_run 2> ops_run.log &&
  [ "$(./ops_run)" = "1 2 1 | 1 2 -2" ] && [ ! -s ops_run.log ] &&
  "$cwrap" -translate ops_run.c -o ops_run.out &&
  grep -q 'switch ( id ) { case 1 : __cwrap_prolongation_p_twice ( ( pr ) ) ; break ; case 2 : __cwrap_prolongation_p_neg' ops_run.out &&
  ! grep -q '(\*' ops_run.out
check "@prolongation and @restriction: ids by kind, attached by name, a switch with no function pointer" $?

"$cwrap" -c ops_bad.c -o ops_bad.o 2> ops_bad.log
[ $? -ne 0 ] && [ ! -e ops_bad.o ] &&
  grep -q "^ops_bad.c:2: error: @boundary b_one is defined twice" ops_bad.log &&
  grep -q "^ops_bad.c:1: note: 'b_one' is defined here" ops_bad.log &&
  grep -q "^ops_bad.c:5: error: @restriction (...) (...): no @restriction is defined before it" ops_bad.log
check "ghost operators: defined twice, or dispatched with none defined, an error" $?

"$cwrap" -dimension=2 -Wall -Wextra ops_coupled.c -o ops_coupled 2> ops_coupled.log &&
  [ "$(./ops_coupled)" = "34 1 2 0 | 0 -1 2 2 1 2 3 1 0 2 -1" ] && [ ! -s ops_coupled.log ]
check "coupled operators: (scalar a, scalar b) passed from uses, @prolongation uses (id), ghost_ops () with globals, boundary.top, boundary.all" $?

"$cwrap" -dimension=2 -c ops_coupled_bad.c -o ops_coupled_bad.o 2> ops_coupled_bad.log
[ $? -ne 0 ] && [ ! -e ops_coupled_bad.o ] &&
  [ "$(grep -o '^ops_coupled_bad.c:[0-9]*: error' ops_coupled_bad.log | cut -d: -f2 | tr '\n' ' ')" = "5 9 10 11 12 13 14 " ] &&
  grep -q "^ops_coupled_bad.c:5: error: @prolongation wrong: what it depends on is scalars" ops_coupled_bad.log &&
  grep -q "^ops_coupled_bad.c:9: error: @prolongation (id) (...): some @prolongation depends on scalars" ops_coupled_bad.log &&
  grep -q "^ops_coupled_bad.c:10: error: prolongation (needs, ...): 'local', a scalar it depends on, must be global" ops_coupled_bad.log &&
  grep -q "^ops_coupled_bad.c:11: error: ghost_ops (zb, ...): 'smoothing(needs)' is not prolongation" ops_coupled_bad.log &&
  grep -q "^ops_coupled_bad.c:12: error: ghost_ops (zb, ...): 'boundary(b_one)': a boundary condition is bound per face" ops_coupled_bad.log &&
  grep -q "^ops_coupled_bad.c:13: error: ghost_ops (zb, ...): 'boundary.back': no face 'back' in 2 dimensions" ops_coupled_bad.log &&
  grep -q "^ops_coupled_bad.c:14: error: ghost_ops (zb, ...): 'boundary.middle': no face 'middle'; the faces are" ops_coupled_bad.log &&
  grep -q "^ops_coupled_bad.c:15: warning: ghost_ops (zb, ...): 'boundary.all' comes after 'boundary.left'" ops_coupled_bad.log
check "coupled operators: a dependency not a scalar, a dispatch without uses, a local dependency, an item of no kind, boundary with no face, a face beyond the dimension, an unknown face, .all after a face" $?

"$cwrap" -dimension=2 -c scalar_global_free.c -o scalar_global_free.o 2> scalar_global_free.log
[ $? -ne 0 ] && [ ! -e scalar_global_free.o ] &&
  [ "$(grep -o '^scalar_global_free.c:[0-9]*: error' scalar_global_free.log | cut -d: -f2 | tr '\n' ' ')" = "9 10 11 12 13 " ] &&
  grep -q "^scalar_global_free.c:9: error: scalar_free () of 'zb', a global scalar: it is numbered once for the whole run and never freed" scalar_global_free.log &&
  grep -q "^scalar_global_free.c:10: error: scalar_free () of 'u.x', a member of a global vector" scalar_global_free.log &&
  grep -q "^scalar_global_free.c:1: note: 'zb' is declared here" scalar_global_free.log
check "a global scalar is never freed: at file scope, static, extern, a member, a tensor; a local one is" $?

"$cwrap" -dimension=2 -Wall scalar_global_abort.c -o scalar_global_abort 2> scalar_global_abort.log &&
  (ulimit -c 0; ./scalar_global_abort > scalar_global_abort.out 2> scalar_global_abort.err
   echo $? > scalar_global_abort.rc) 2> /dev/null &&
  [ "$(cat scalar_global_abort.rc)" -ne 0 ] && [ ! -s scalar_global_abort.out ] &&
  grep -q "^scalar_free: scalar zb (0) is global, numbered once for the whole run: it is never freed" scalar_global_abort.err
check "a global scalar freed through a pointer stops the run" $?

"$cwrap" -dimension=2 -Wall -Wextra scalar_names.c -o scalar_names 2> scalar_names.log &&
  [ "$(./scalar_names)" = "h zb u.x u.y | inner.h inner.h.2 1 1 inner.kept tmp 1 inner.u.x inner.u.y | - -1 -" ] &&
  [ ! -s scalar_names.log ]
check "scalar names: as declared, a local hiding a scalar qualified, deeper numbered, static its function's, found; scalar_new () unnamed" $?

"$cwrap" -dimension=2 -c stencil_foreach_bad.c -o stencil_foreach_bad.o 2> stencil_foreach_bad.log
[ $? -ne 0 ] && [ ! -e stencil_foreach_bad.o ] &&
  [ "$(grep -o '^stencil_foreach_bad.c:[0-9]*: error' stencil_foreach_bad.log | cut -d: -f2 | tr '\n' ' ')" = "6 7 9 " ] &&
  grep -q "^stencil_foreach_bad.c:6: error: again () inside foreach (): stencil loops do not nest" stencil_foreach_bad.log &&
  grep -q "^stencil_foreach_bad.c:7: error: again () takes 0 arguments, not 1" stencil_foreach_bad.log &&
  grep -q "^stencil_foreach_bad.c:9: error: @foreach foreach is defined again, differently" stencil_foreach_bad.log
check "@foreach: one stencil loop in another, the arguments counted, a different definition again" $?

"$cwrap" -dimension=2 -c stencil_nest.c -o stencil_nest.o 2> stencil_nest.log
[ $? -ne 0 ] && [ ! -e stencil_nest.o ] &&
  [ "$(grep -o '^stencil_nest.c:[0-9]*: error' stencil_nest.log | cut -d: -f2 | tr '\n' ' ')" = "8 9 " ] &&
  grep -q "^stencil_nest.c:8: error: 'fill ()' runs a stencil loop, and is called inside one" stencil_nest.log &&
  grep -q "^stencil_nest.c:3: note: 'fill' runs one here" stencil_nest.log &&
  grep -q "^stencil_nest.c:4: note: 'twice' calls one that runs one" stencil_nest.log
check "@foreach: a call inside a stencil loop to a function that runs one, directly or not" $?

"$cwrap" -dimension=2 -Wall -Wextra stencil_ghosts.c -o stencil_ghosts 2> stencil_ghosts.log &&
  [ "$(./stencil_ghosts)" = "updated 0 4 | stale 14" ] && [ ! -s stencil_ghosts.log ]
check "foreach (): the ghosts of a scalar read at an offset, updated only if written since; A[k] stands for all of A" $?

"$cwrap" -dimension=2 -Wall -Wextra stencil_by_caller.c -o stencil_by_caller 2> stencil_by_caller.log &&
  [ "$(./stencil_by_caller)" = "level: updated | stale 0; then: updated | stale 1" ] &&
  [ ! -s stencil_by_caller.log ]
check "a loop whose block loop gives ghosts_by_caller: no ghost update or mark around it, foreach () still marking" $?

"$cwrap" -dimension=2 -Wall -c stencil_by_caller_bad.c -o stencil_by_caller_bad.o 2> stencil_by_caller_bad.log
[ $? -ne 0 ] && [ ! -e stencil_by_caller_bad.o ] &&
  [ "$(grep -o '^stencil_by_caller_bad.c:[0-9]*: error' stencil_by_caller_bad.log | cut -d: -f2 | tr '\n' ' ')" = "6 " ] &&
  grep -q "^stencil_by_caller_bad.c:6: error: 'p\[1\]' writes a cell other than its own" stencil_by_caller_bad.log &&
  grep -q "^stencil_by_caller_bad.c:9: warning: 'p\[1\]' reads a neighbour of 'p', which this foreach () writes" stencil_by_caller_bad.log
check "ghosts_by_caller keeps the rules: a write at an offset an error, a neighbour of a scalar written a warning" $?

"$cwrap" -dimension=2 -Wall -c stencil_dep.c -o stencil_dep.o 2> stencil_dep.log &&
  [ "$(grep -c 'warning:' stencil_dep.log)" = 1 ] &&
  grep -q "^stencil_dep.c:7: warning: 'p\[1,0\]' reads a neighbour of 'p', which this foreach () writes" stencil_dep.log &&
  grep -q "^stencil_dep.c:6: note: 'p\[\]' is written here" stencil_dep.log
check "foreach (): reading a neighbour of a scalar it writes, a warning" $?

"$cwrap" -dimension=2 -c stencil_nopos.c -o stencil_nopos.o 2> stencil_nopos.log
[ $? -ne 0 ] && grep -q "^stencil_nopos.c:6: error: x, y, z and Delta in foreach () need" stencil_nopos.log
check "foreach (): x without the leaf's and the cell's indices, an error" $?

"$cwrap" -dimension=2 -Wall stencil_call.c -o stencil_call 2> stencil_call.log &&
  [ "$(./stencil_call)" = 2 ] && [ ! -s stencil_call.log ]
check "foreach (): a function named foreach is still called" $?

"$cwrap" -dimension=1 -translate loop_run.c -o loop_run.out 2> /dev/null &&
  ! grep -q '@loop' loop_run.out &&
  grep -q '^  { int __cwrap_repeat_n_[0-9]* = 3 ;$' loop_run.out &&
  grep -q '^  int __cwrap_range_count_[0-9]* = 0;$' loop_run.out &&
  grep -q '^    count++;$' loop_run.out &&
  grep -q '^# [0-9]* "loop.h"$' loop_run.out
check "@loop: -translate writes the expansion with linemarkers, and no definition" $?

"$cwrap" -dimension=1 -Wall -Wextra loop_nest.c -o loop_nest 2> loop_nest.log &&
  [ "$(./loop_nest)" = "24 checks, 0 failed" ] && [ ! -s loop_nest.log ]
check "@loop nested: in blocks, in bodies three levels deep, names passed on, break, statics, caches, scalars" $?

"$cwrap" -c loop_nest_bad.c -o loop_nest_bad.o 2> loop_nest_bad.log
[ $? -ne 0 ] && [ ! -e loop_nest_bad.o ] &&
  [ "$(grep -o '^loop_nest_bad.c:[0-9]*: error' loop_nest_bad.log | cut -d: -f2 | tr '\n' ' ')" = "2 3 4 5 7 " ] &&
  grep -q "^loop_nest_bad.c:2: error: 'k' is not declared by the loop where {...} is" loop_nest_bad.log &&
  grep -q "^loop_nest_bad.c:5: error: loop 'itself' is used in its own definition" loop_nest_bad.log
check "@loop nested: misuse inside bodies and blocks is an error" $?

"$cwrap" -dimension=1 -c loop_nest_lines.c -o loop_nest_lines.o 2> loop_nest_lines.log
[ $? -ne 0 ] && grep -q "^loop_nest_lines.c:9:[0-9]*: error: .nope1. undeclared" loop_nest_lines.log &&
  grep -q "^loop_nest_lines.c:11:[0-9]*: error: .nope2. undeclared" loop_nest_lines.log &&
  grep -q "^loop_nest_lines.h:7:[0-9]*: error: .missing_name. undeclared" loop_nest_lines.log &&
  grep -q "^loop_nest_lines.c:14:[0-9]*: error: .nope3. undeclared" loop_nest_lines.log
check "@loop nested: errors keep their lines, in the block, in a body, and after" $?

"$cwrap" -c loop_pragma.c -o loop_pragma.o 2> loop_pragma.log
[ $? -ne 0 ] && [ "$(grep -c 'error' loop_pragma.log)" = 1 ] &&
  grep -q "^loop_pragma.c:18:[0-9]*: error: .nope. undeclared" loop_pragma.log &&
  sed 's/ + nope//' loop_pragma.c > loop_pragma2.c && "$cwrap" -Wall loop_pragma2.c -o loop_pragma &&
  { ./loop_pragma; [ $? = 12 ]; } &&
  "$cwrap" -translate loop_pragma2.c -o loop_pragma.out &&
  [ "$(grep -c '^#pragma GCC ivdep' loop_pragma.out)" = 2 ]
check "@loop: a pragma in a body goes with each use, and lines stay right" $?

"$cwrap" -D_GNU_SOURCE -Wall floatn.c -o floatn -lm 2> floatn.log && ./floatn &&
  ! grep -q 'not translated' floatn.log
check "_Float32, _Complex _Float64 and the like parse, as <complex.h> declares them" $?

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
