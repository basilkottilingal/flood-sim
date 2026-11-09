#!/bin/sh
# memindex against voxel-memindex on 2:1 balanced trees: run from libb-amr:
#   sh tests/bench_balanced.sh [RAW]
# Each case runs once per implementation, in its own process; the two must
# agree (checksum). RAW, if given, keeps the key=value lines. memindex
# needs up to about 1 GB for the cases below.
lib=$(pwd)
dir=$(mktemp -d)
trap 'rm -rf "$dir"' EXIT
raw=${1:-$dir/raw}
: > "$raw"

for d in 2 3; do
  gcc -O2 -std=gnu11 -Wall -Wextra -Werror -Ddimension=$d -D_GNU_SOURCE -I"$lib" \
    tests/bench_balanced.c memindex.c voxel-memindex.c pool.c -lm \
    -o "$dir/bench$d" || exit 1
done

# dimension pattern depth
cases="
2 dense 10
2 shell 15
2 river 15
3 dense 6
3 shell 8
3 river 10
3 shell 9
"

echo "$cases" | while read d pat depth; do
  [ -z "$d" ] && continue
  for impl in mem voxel; do
    if ! "$dir/bench$d" $impl $pat $depth >> "$raw" 2> "$dir/err"; then
      echo "FAIL $impl ${d}D $pat depth $depth (killed: out of memory?)"; cat "$dir/err"
    fi
  done
done

awk '
function kv(   i, a) { delete v; for (i = 1; i <= NF; i++) { split ($i, a, "="); v[a[1]] = a[2] } }
function c(x) { return (x in v) ? sprintf ("%.1f", v[x]) : "-" }
BEGIN {
  printf "| case | impl | nodes (leaves) | B/node | build ns | dfs ns | dfs+cur ns | list ns | list+cur ns | nbr ns/leaf | nbr+cur ns/leaf | ns/lookup (+cur) | free ns |\n"
  printf "|---|---|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|\n"
}
{
  kv()
  k = v["dim"] "D " v["pattern"] " depth " v["depth"]
  printf "| %s | %s | %d (%d) | %.0f | %.1f | %.1f | %s | %.1f | %s | %.0f | %s | %.1f (%s) | %.1f |\n",
    k, v["impl"], v["nodes"], v["leaves"], v["bytes_per_node"], v["build_ns"],
    v["dfs_ns"], c("dfscur_ns"), v["list_ns"], c("listcur_ns"),
    v["nbr_leaf_ns"], c("nbrcur_leaf_ns"), v["nbr_lookup_ns"], c("nbrcur_lookup_ns"),
    v["free_ns"]
  if (v["impl"] == "mem") { sum[k] = v["checksum"]; look[k] = v["lookups_per_leaf"]; co[k] = v["coarser_share"] }
  else if ((k in sum) && sum[k] != v["checksum"]) bad = bad " " k
}
END {
  if (bad != "") print "\nCHECKSUMS DIFFER:" bad
  else print "\nchecksums agree in every case both finished"
  printf "\nneighbour lookups per leaf (coarser fall-backs share):"
  for (k in look) printf " %s: %.1f (%.0f%%);", k, look[k], 100 * co[k]
  print ""
}' "$raw"
