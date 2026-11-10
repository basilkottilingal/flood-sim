#!/bin/sh
# memindex against voxel-memindex and nano-vdb (dense and compact):
# run from libb-amr:  sh tests/bench_memindex.sh [RAW]
# Each case runs once per implementation, in its own process. All must
# agree (checksum). RAW, if given, keeps the key=value lines.
# Cache and dTLB misses need kernel.perf_event_paranoid <= 2, e.g.
#   sudo sysctl kernel.perf_event_paranoid=1
# The sparse cases cost memindex a page and a page-table page per row
# touched: about 0.5-1.5 GB at the sizes below.
lib=$(pwd)
dir=$(mktemp -d)
trap 'rm -rf "$dir"' EXIT
raw=${1:-$dir/raw}
: > "$raw"

# twice: nano-vdb's dense nodes, and its compact ones (-DNVDB_COMPACT)
for d in 2 3; do
  for c in "" "-DNVDB_COMPACT"; do
    gcc -O2 -mpopcnt -std=gnu11 -Wall -Wextra -Werror -Ddimension=$d -D_GNU_SOURCE \
      -I"$lib" $c tests/bench_memindex.c memindex.c voxel-memindex.c nano-vdb.c pool.c -lm \
      -o "$dir/bench$d${c:+c}" || exit 1
  done
done

# dimension pattern level [N]
cases="
2 dense  11
2 shell  16
2 river  16 32
2 random 14 1000000
2 random 20 50000
3 dense  7
3 shell  9
3 shell  10
3 river  12 16
3 random 12 100000
3 random 20 30000
"

echo "$cases" | while read d pat lev n; do
  [ -z "$d" ] && continue
  for impl in mem voxel nvdb nvdb-compact; do
    bin="$dir/bench$d" name=$impl
    [ $impl = nvdb-compact ] && bin="$dir/bench${d}c" name=nvdb
    # nano-vdb's dense nodes on scattered 3D points: about 300 KB a point
    if [ $impl = nvdb ] && [ $d = 3 ] && [ $pat = random ]; then
      echo "skip nvdb ${d}D $pat L$lev: dense VDB nodes, ~300 KB per isolated point"
      continue
    fi
    if ! "$bin" $name $pat $lev $n >> "$raw" 2> "$dir/err"; then
      # memindex runs out of memory first on the large sparse cases
      echo "FAIL $impl ${d}D $pat L$lev $n (killed: out of memory?)"; cat "$dir/err"
    fi
  done
done

awk '
function kv(   i, a) { delete v; for (i = 1; i <= NF; i++) { split ($i, a, "="); v[a[1]] = a[2] } }
BEGIN {
  printf "| case | impl | points | B/point | build ns | build faults | ordered ns | shuffled ns | nbr ns | nbr+cursor ns | uniform ns (hits) | free ns | kB left after free |\n"
  printf "|---|---|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|\n"
}
{
  kv()
  c = v["dim"] "D " v["pattern"] " L" v["level"]
  printf "| %s | %s | %d | %.0f | %.1f | %d | %.1f | %.1f | %.1f | %s | %.1f (%.0f%%) | %.1f | %d |\n",
    c, v["impl"], v["n"], v["bytes_per_pt"], v["build_ns"], v["build_flt"],
    v["ordered_ns"], v["shuffled_ns"], v["nbr_ns"],
    ("nbrcur_ns" in v) ? sprintf ("%.1f", v["nbrcur_ns"]) : "-",
    v["uniform_ns"], 100 * v["uniform_hits"], v["free_ns"],
    v["rss_after_free_kb"] + v["pte_after_free_kb"]
  if (v["impl"] == "mem") sum[c] = v["checksum"]
  else if ((c in sum) && sum[c] != v["checksum"]) bad = bad " " c
  if (v["ordered_llc"] != "n/a") perf = 1
  line[NR] = $0
}
END {
  if (bad != "") print "\nCHECKSUMS DIFFER:" bad
  else print "\nchecksums agree in every case both finished"
  if (!perf) { print "cache/dTLB misses: not available (kernel.perf_event_paranoid > 2)"; exit }
  printf "\nmisses per lookup, LLC / dTLB\n\n| case | impl | ordered | shuffled | nbr | uniform |\n|---|---|--:|--:|--:|--:|\n"
  for (r = 1; r <= NR; r++) {
    $0 = line[r]; kv()
    printf "| %sD %s L%s | %s | %s / %s | %s / %s | %s / %s | %s / %s |\n",
      v["dim"], v["pattern"], v["level"], v["impl"],
      v["ordered_llc"], v["ordered_tlb"], v["shuffled_llc"], v["shuffled_tlb"],
      v["nbr_llc"], v["nbr_tlb"], v["uniform_llc"], v["uniform_tlb"]
  }
}' "$raw"
