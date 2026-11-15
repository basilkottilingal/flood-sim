#!/bin/sh
# Profiling, not testing: a ghost update piece by piece (restriction, the
# ghosts' copies, level jumps and walls, prolongation into children, a
# stencil) on the CPU and on an NVIDIA GPU, each as fast as it goes
# (tests/bench_gpu.c says how). Run from libb-amr:  sh tests/run2.sh
# The results are checked against the library's on the host, to a
# tolerance, before any timing: not bit for bit (tests/run.sh, GPU=1).
#
# CWRAP: the cwrap to build with (default ../libc-pp/build/cwrap).
# DIMS (2 3), PRECS (double float): what to run.
# ARGS2, ARGS3: base depth blockdepth of the tree (default 4 8 5, 2 5 4):
#   uniform at level base, refined to depth on a shell; blocks of
#   2^blockdepth cells a side. Several trees: ';' between them, as
#   ARGS2="5 9 4; 4 8 5; 3 7 6" (the same finest cell, m 16, 32, 64).
# GPU=0: the CPU only. GPU_ARCH: the PTX target (default sm_80).
# RAW: a file to keep the raw lines in.
# Times are medians of 7 or more calls; the machine should be idle.
cwrap=${CWRAP:-$(pwd)/../libc-pp/build/cwrap}
lib=$(pwd)
dir=$(mktemp -d)
trap 'rm -rf "$dir"' EXIT
raw=${RAW:-$dir/raw}
: > "$raw"
src="$lib/tree.c $lib/grid-scalar.c $lib/memindex.c $lib/voxel-memindex.c $lib/nano-vdb.c $lib/pool.c $lib/adapt.c $lib/../libutils/bits.c"
base="-D_GNU_SOURCE -std=gnu11 -O3 -DNDEBUG -Wall -Wextra -Wno-unused-parameter -I$lib -I$lib/.."
# the CPU: every thread, the widest vectors; offload off (a gcc built with
# --enable-offload-defaulted offloads to every offload compiler installed)
cpu="$base -march=native -fopenmp -foffload=disable"
# the GPU: as tests/run.sh GPU=1, but fma allowed and the PTX JIT at its
# default optimisation: no bit identity asked for
gpu="$base -march=native -fopenmp -foffload=nvptx-none -foffload-options=nvptx-none=-misa=${GPU_ARCH:-sm_80}"
gpu="$gpu -fcf-protection=none -fno-stack-protector -no-pie"
dims=${DIMS:-2 3}
precs=${PRECS:-double float}

variants=cpu
[ "${GPU:-1}" != 0 ] && variants="$variants gpu"

# every binary at once
fails=0 pids=
for d in $dims; do
  for p in $precs; do
    pf=; [ "$p" = float ] && pf=-DBAMR_SINGLE
    for v in $variants; do
      case $v in
        cpu)  f="$cpu" ;;
        gpu)  f="$gpu" ;;
      esac
      ( "$cwrap" -dimension=$d $f $pf tests/bench_gpu.c $src -lm -o "$dir/$v$d$p" \
          > "$dir/$v$d$p.log" 2>&1 || echo failed >> "$dir/$v$d$p.log" ) &
      pids="$pids $!"
    done
  done
done
wait $pids
for d in $dims; do
  for p in $precs; do
    for v in $variants; do
      if [ ! -x "$dir/$v$d$p" ] || grep -q 'warning\|error\|failed' "$dir/$v$d$p.log"; then
        echo "FAIL build $v ${d}D $p"; cat "$dir/$v$d$p.log"; fails=$((fails + 1))
      fi
    done
  done
done
[ "$fails" = 0 ] || exit "$fails"

# one at a time: they would compete. ARGSd may list several trees, ';'
# between them: each a case of its own, tagged base-depth-blockdepth
for d in $dims; do
  eval list=\${ARGS$d:-}
  [ -z "$list" ] && { [ $d = 2 ] && list="4 8 5" || list="2 5 4"; }
  echo "$list" | tr ';' '\n' > "$dir/cases"
  while read -r args; do
    [ -z "$args" ] && continue
    tag=$(echo $args | tr ' ' -)
    for p in $precs; do
      for v in $variants; do
        mode=dev; [ $v = cpu ] && mode=cpu
        if [ $v = cpu ]; then
          "$dir/$v$d$p" $mode $args > "$dir/out" 2>&1 < /dev/null
        else
          OMP_TARGET_OFFLOAD=MANDATORY "$dir/$v$d$p" $mode $args > "$dir/out" 2>&1 < /dev/null
        fi
        st=$?
        grep '^T\|^C' "$dir/out" | sort -u | sed "s/^/  [$v $tag] /"
        sed -n "s/^R \([^ ]*\) \([^ ]*\) [^ ]* /R \1 \2,$tag $v /p" "$dir/out" >> "$raw"
        if [ $st != 0 ]; then
          echo "FAIL $v ${d}D $p $tag (exit $st)"; grep -v '^[RTC] ' "$dir/out" | tail -5
          fails=$((fails + 1))
        fi
      done
    done
  done < "$dir/cases"
done

echo
awk '
$1 == "R" {
  k = $2 "D " $3 " | " $5
  if (!(k in seen)) { seen [k] = 1; order [n++] = k }
  ms [k, $4] = $6; cells [k] = $7
}
function col(k, v) { return ((k, v) in ms) ? sprintf ("%.3f", ms [k, v]) : "-" }
function ns(k, v) { return ((k, v) in ms) && cells [k] > 0 ? sprintf ("%.2f", 1e6 * ms [k, v] / cells [k]) : "-" }
END {
  printf "| case | phase | cells | CPU ms | CPU ns/cell | GPU ms | GPU ns/cell | CPU / GPU |\n"
  printf "|---|---|--:|--:|--:|--:|--:|--:|\n"
  for (i = 0; i < n; i++) {
    k = order [i]
    r = ((k, "cpu") in ms) && ((k, "gpu") in ms) && ms [k, "gpu"] > 0 ? sprintf ("%.1f", ms [k, "cpu"] / ms [k, "gpu"]) : "-"
    printf "| %s | %d | %s | %s | %s | %s | %s |\n", k, cells [k], col(k, "cpu"),
      ns(k, "cpu"), col(k, "gpu"), ns(k, "gpu"), r
  }
}' "$raw"
cat <<EOF

ghosts:jumps and ghosts:walls on the CPU are differences (the fill with
them, less without): noisy where small. ghosts = the fill (copies, jumps,
walls); ghosts_update = restriction + ghosts; step = ghosts_update +
stencil. upload, download: every block, host to device and back, as a
ghost update on the host would need. CPU: $(nproc) threads.
EOF
[ "$fails" = 0 ] || exit "$fails"
