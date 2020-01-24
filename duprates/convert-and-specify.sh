#!/usr/bin/env bash

set -euo pipefail
IFS=$'\t\n'


genetree=$1
subtree=$2
subtreedir=${3:-subtreesGoodQualO2}


SCRIPTS=$HOME/scripts


subtree_input="./$genetree/$subtreedir/${subtree}_codeml.nwk"
subtree_species="./$genetree/$subtreedir/${subtree}_codeml.species.nwk"
$SCRIPTS/seqtools/specify.py -e 93 -f nwk -l '{sp}_{gene}' "$subtree_input" "$subtree_species"

al_root="./$genetree/$subtreedir/realign/${subtree}"
al_input="${al_root}_protfsa.fa"
al_phy="${al_root}_protfsa.species.phy"
# Prefix the species names + convert the AA alignment to phylip
# 'phylip-relaxed' works with pb too.
$SCRIPTS/seqtools/specify.py -e 93 -f fasta -l '{sp}_{gene}' "$al_input" |\
    $SCRIPTS/seqtools/seqconv.py -t 'phylip-sequential-relaxed' > "$al_phy"


