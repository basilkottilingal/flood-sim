#!/usr/bin/env bash
#
# depgraph.sh - draw the dependency graph of lib*/ directories from the
#               DEPS variable set in each lib*/Makefile.
#
# Usage:
#   ./depgraph.sh                  # DOT to stdout
#   ./depgraph.sh -o deps.png      # render (format from extension: png, svg, pdf)
#   ./depgraph.sh -o deps.dot      # write DOT file
#   ./depgraph.sh -r path/to/root  # project root (default: current dir)
#   ./depgraph.sh -l               # print link order (dependents first) and exit
#
# An edge A -> B means "A depends on B".
# A name in DEPS is resolved as lib<name>/ first, then <name>/.

set -euo pipefail

root=.
out=
link_order=0

while getopts "o:r:lh" opt; do
    case $opt in
        o) out=$OPTARG ;;
        r) root=$OPTARG ;;
        l) link_order=1 ;;
        h|*) sed -n '2,13p' "$0" | sed 's/^# \{0,1\}//'; exit 0 ;;
    esac
done

cd "$root"

# Ask make for DEPS of one directory (handles +=, :=, continuations, ifdefs).
get_deps() {
    local dir=$1
    make -s --no-print-directory -C "$dir" \
        -f Makefile -f <(printf '__print_deps:\n\t@echo $(DEPS)\n') \
        __print_deps 2>/dev/null || true
}

# Resolve a DEPS entry to an existing directory name.
resolve() {
    local d=$1
    if   [[ -d lib$d ]]; then echo "lib$d"
    elif [[ -d $d    ]]; then echo "$d"
    else
        echo "depgraph: warning: dependency '$d' has no directory" >&2
        echo "lib$d"
    fi
}

nodes=()
edges=()   # "from to" pairs

for mk in lib*/Makefile; do
    [[ -e $mk ]] || continue
    dir=${mk%/Makefile}
    nodes+=("$dir")
    for d in $(get_deps "$dir"); do
        edges+=("$dir $(resolve "$d")")
    done
done

if ((${#nodes[@]} == 0)); then
    echo "depgraph: no lib*/Makefile found under $(pwd)" >&2
    exit 1
fi

# Link order via tsort: dependents come before the libraries they use.
# tsort also detects cycles (prints "input contains a loop" on stderr).
link_order_list() {
    {
        for n in "${nodes[@]}"; do echo "$n $n"; done      # keep isolated nodes
        for e in "${edges[@]}"; do echo "$e"; done
    } | tsort
}

if ((link_order)); then
    link_order_list
    exit 0
fi

# Cycle check (warn only; the graph is still drawn).
if ! link_order_list >/dev/null 2>/tmp/depgraph_tsort.$$; then
    echo "depgraph: WARNING: circular dependency detected:" >&2
    cat /tmp/depgraph_tsort.$$ >&2
fi
rm -f /tmp/depgraph_tsort.$$

emit_dot() {
    echo 'digraph deps {'
    echo '    rankdir=LR;'
    echo '    node [shape=box, style=rounded, fontname="monospace"];'
    for n in "${nodes[@]}"; do
        echo "    \"$n\";"
    done
    for e in "${edges[@]}"; do
        read -r from to <<<"$e"
        echo "    \"$from\" -> \"$to\";"
    done
    echo '}'
}

if [[ -z $out ]]; then
    emit_dot
elif [[ $out == *.dot ]]; then
    emit_dot >"$out"
    echo "wrote $out" >&2
else
    command -v dot >/dev/null || {
        echo "depgraph: 'dot' (graphviz) not installed; writing DOT to stdout instead" >&2
        emit_dot
        exit 0
    }
    emit_dot | dot -T"${out##*.}" -o "$out"
    echo "wrote $out" >&2
fi
