#!/usr/bin/env bash

set -euo pipefail

help="
Convert a newick tree with PhylTree notations (pipes and spaces in labels, multiline, internal node names) to a 'standard' newick format

    ${0} [-s space_repl] [-b] inputtree

-s  character to replace spaces in labels
-b  even more basic: no internal node labels, no root length.
"

space_repl='.'
basic=0

while getopts 's:b' opt; do
    case "$opt" in
        s) space_repl="$OPTARG" ;;
        b) basic=1 ;;
        *) echo "$help"; exit 2 ;;
    esac
done
shift $((OPTIND-1))

inputtree=${1:-}

[ -n "$inputtree" ] || { echo "$help"; exit 2; }


if [ "$basic" -eq 0 ]; then

sed -r 's/([A-Za-z0-9]+) /\1'"${space_repl}"'/g; s/\|[A-Za-z0-9_.|-]+:/:/' \
    "${inputtree}" |\
    sed ':a; N;s/\n\s*//; ta'

else
# Additionally, remove internal node labels, and the root branch length.
sed -r -e 's/([A-Za-z0-9]+) /\1./g' \
    -e 's/\|[A-Za-z0-9_.|-]+:/:/' \
    -e 's/\)[A-Za-z0-9_.-]+([:;])/)\1/' \
    -e 's/:[0-9.-];/;/' \
    "${inputtree}" |\
    sed ':a; N;s/\n\s*//; ta'
done
