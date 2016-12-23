#!/bin/bash

set -euo pipefail
trap exit EXIT

help="USAGE:
    ./reshape_tree.sh [-h] <input_tree>

OPTIONS:
    -h: print this help and exit

DESCRIPTION:
    prints reshaped newick tree to stdout.
    Reshaping means removing the inner node labels, and removing nodes that
    start from the root and have only one child (codeml fails on such a tree).
    Needs newicktools installed.
"

input_tree=${1:-}

if [ -z "$input_tree" ] ; then
	echo "$help"
	exit 1
elif [ "$input_tree" = "-h" ]; then
	echo "$help"
	exit
fi

input_labels="$(nw_labels -I $input_tree)"

nw_topology -bI "$input_tree" | nw_clade - $input_labels

