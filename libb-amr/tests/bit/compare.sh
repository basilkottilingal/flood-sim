#!/bin/sh
# Every cell bit for bit: the working tree against a git ref (HEAD by
# default), in 1, 2 and 3 dimensions, double and float. Run from libb-amr:
#
#   sh tests/bit/compare.sh            # against HEAD
#   sh tests/bit/compare.sh main       # against another ref
#   TIMING=1 sh tests/bit/compare.sh   # also the time of ghosts_update_levels
#
# checksum.c refines, updates ghosts with given operators and adapts, then
# hashes every cell of every block, ghosts included. Both versions are built
# from this checksum.c: a ref whose API it does not build against fails to
# build, said so. Exit status: the number of builds that differ or fail.
# CWRAP: the cwrap to build with (default ../libc-pp/build/cwrap).
ref=${1:-HEAD}
cwrap=${CWRAP:-$(pwd)/../libc-pp/build/cwrap}
new=$(pwd)
here=$new/tests/bit
dir=$(mktemp -d)
trap 'rm -rf "$dir"' EXIT
top=$(git rev-parse --show-toplevel)
git -C "$top" archive "$ref" libb-amr libutils | tar -x -C "$dir" || {
  echo "compare: no ref '$ref'"; exit 1; }
old=$dir/libb-amr
bad=0
for d in 1 2 3; do
  for prec in "" -DBAMR_SINGLE; do
    for v in old new; do
      eval L=\$$v
      # the library's sources, those of this version (a ref may have others)
      src=
      for f in tree grid-scalar morton-curve memindex voxel-memindex nano-vdb pool block-ghosts adapt ../libutils/bits; do
        [ -e "$L/$f.c" ] && src="$src $L/$f.c"
      done
      progs=checksum; [ -n "$TIMING" ] && progs="checksum timing"
      for p in $progs; do
        "$cwrap" -dimension=$d -D_GNU_SOURCE -std=gnu11 -O2 $prec -I"$L" -I"$L/.." \
          "$here/$p.c" $src -lm -o "$dir/${p}_$v" 2> "$dir/log" ||
          { n=$ref; [ $v = new ] && n="the working tree"
            echo "${d}D ${prec:-double}: $p failed to build against $n"; cat "$dir/log"; }
      done
      eval "h_$v=\"\$(\"$dir/checksum_$v\" 2>/dev/null)\""
      [ -n "$TIMING" ] && eval "t_$v=\"\$(\"$dir/timing_$v\" 2>/dev/null)\""
      rm -f "$dir/checksum_$v" "$dir/timing_$v"
    done
    if [ -n "$h_old" ] && [ "$h_old" = "$h_new" ]; then r=same
    else r=DIFFERENT; bad=$((bad + 1)); fi
    echo "${d}D ${prec:-double}: $r  $ref $h_old  new $h_new"
    [ -n "$TIMING" ] && echo "    $ref: $t_old" && echo "    new:  $t_new"
  done
done
exit "$bad"
