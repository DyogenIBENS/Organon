#!/usr/bin/env bash

set -euo pipefail
IFS=$'\t\n'


genetree=$1
subtree=$2
suffix=${3:-species}
subtreedir=${4:-subtreesGoodQualO2}

SCRIPTS=$HOME/scripts


subtree_input="./$genetree/$subtreedir/${subtree}_codeml.nwk"
subtree_species="./$genetree/$subtreedir/${subtree}_codeml.${suffix}.nwk"
$SCRIPTS/seqtools/specify.py -e 93 -f nwk -l '{shortsp}_{gene}' \
    -t 'gene=s/^ENS(G|[A-Z]{3}G)0{5}/G/' "$subtree_input" "$subtree_species"

# Use this option to significantly shorten most sequence names: -t 'gene=s/^ENS(G|[A-Z]{3}G)0{5}/G/'

al_root="./$genetree/$subtreedir/realign/${subtree}"
al_input="${al_root}_protfsa.fa"
al_phy="${al_root}_protfsa.${suffix}.phy"
# Prefix the species names + convert the AA alignment to phylip
# 'phylip-relaxed' works with pb too.
$SCRIPTS/seqtools/specify.py -e 93 -f fasta -l '{shortsp}_{gene}' \
    -t 'gene=s/^ENS(G|[A-Z]{3}G)0{5}/G/' "$al_input" |\
    $SCRIPTS/seqtools/seqconv.py -t 'phylip-sequential-relaxed' > "$al_phy"


