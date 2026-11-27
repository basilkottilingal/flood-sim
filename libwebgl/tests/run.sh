#!/bin/sh
# libwebgl's tests, from libwebgl: sh tests/run.sh
#   test_webgl.c   the scene written and read back: every layer and field
#                  kind, the calls refused, escaping, every blob bit for
#                  bit, three.js and the viewer in the file (plain C);
#                  again with no three.js in the library (WEBGL_ONLINE=1:
#                  the import map's jsDelivr URLs)
#   test_scene.c   the adapters, 2D: a libb-amr tree and a libhmesh circle,
#                  against the cells' function and the mesh's own data,
#                  double and float
#   check_json.py  the manifest and the import map, by Python's parser
#   viewer.js      parsed as a module by gjs (SpiderMonkey, Firefox's
#                  engine) when gjs is there: syntax only
# All under AddressSanitizer and UBSan, failing on any warning. The HTML
# files are left in tests/, to open in a browser.
# CWRAP: the cwrap to build with (default ../libc-pp/build/cwrap, else cwrap
# on the PATH).
cwrap=${CWRAP:-$(pwd)/../libc-pp/build/cwrap}
[ -x "$cwrap" ] || cwrap=cwrap
lib=$(pwd)
amr=$lib/../libb-amr
hm=$lib/../libhmesh
dir=$(mktemp -d)
trap 'rm -rf "$dir"' EXIT
sh gen-assets.sh || exit 1
san="-O1 -g -fsanitize=address,undefined -fno-sanitize-recover=all"
warn="-Wall -Wextra -Wshadow -Wconversion -Wsign-conversion -Wno-float-conversion"
hsrc="$hm/hmesh.c $hm/hmesh-clip.c $hm/hmesh-morton.c $hm/hmesh-tree.c $hm/hmesh-vof.c $hm/hmesh-implicit.c $lib/../libutils/bits.c"
asrc="$amr/tree.c $amr/grid-scalar.c $amr/memindex.c $amr/voxel-memindex.c $amr/nano-vdb.c $amr/pool.c $amr/block-ghosts.c $amr/adapt.c"
fails=0

run () {   # name, compiler, html, sources and flags...
  name=$1 cc=$2 html=$3
  shift 3
  out=
  if $cc "$@" -lm -o "$dir/t" 2> "$dir/log" &&
     ! grep -v 'mem_destroy\|^In file\|^ *[0-9]* |\|^ *|' "$dir/log" | grep -q 'warning\|error' &&
     out=$("$dir/t" "$html" 2>&1) &&
     json=$(python3 -I tests/check_json.py "$html" 2>&1); then
    echo "ok   $name: $(echo "$out" | tail -1); $json"
  else
    echo "FAIL $name"; cat "$dir/log"; echo "$out" | tail -12; echo "$json"
    fails=$((fails + 1))
  fi
}

run "scene, core" gcc tests/test_webgl.html -std=gnu11 -D_GNU_SOURCE $warn $san \
  -I. -Itests tests/test_webgl.c webgl.c webgl-assets.c
# without three.js in the library: the pages load it from jsDelivr
WEBGL_ONLINE=1 sh gen-assets.sh "$dir/assets-online.c" || exit 1
run "scene, core, online" gcc tests/test_webgl_online.html -std=gnu11 \
  -D_GNU_SOURCE $warn $san -I. -Itests tests/test_webgl.c webgl.c \
  "$dir/assets-online.c"
for prec in "" "-DBAMR_SINGLE"; do
  run "scene, tree and hmesh${prec:+, float}" "$cwrap -dimension=2" \
    tests/test_scene${prec:+_float}.html -D_GNU_SOURCE -std=gnu11 -Wall \
    -Wextra -Wno-unused-parameter $san $prec -I. -Itests -I$amr -I$hm \
    -I$lib/.. tests/test_scene.c webgl.c webgl-tree.c webgl-hmesh.c \
    webgl-assets.c $hsrc $asrc
done

if command -v gjs > /dev/null; then
  cp viewer/viewer.js "$dir/viewer.mjs"
  if HOME=$dir gjs -m "$dir/viewer.mjs" 2>&1 | grep -q "SyntaxError"; then
    echo "FAIL viewer.js: a syntax error"
    HOME=$dir gjs -m "$dir/viewer.mjs" 2>&1 | grep SyntaxError
    fails=$((fails + 1))
  else
    echo "ok   viewer.js parses (gjs)"
  fi
fi
[ "$fails" = 0 ] && echo "all passed"
exit "$fails"
