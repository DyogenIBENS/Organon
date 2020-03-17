#!/usr/bin/env bash

set -euo pipefail
IFS=$'\t\n'

help="
USAGE

./convert-and-specify.sh <genetree> <subtree> [ <suffix> <altype> <subtreedir> ]

ARGUMENTS

  genetree: e.g. ENSGT00390000000002
  subtree: e.g. RodentiaENSGT00390000000002.a.a.a.a.a.a
  suffix: for the 'specified' alignment [species]
  altype: 'prot' or 'cds' [prot]
  subtreedir: [subtreesGoodQualO2]
"

genetree=$1
subtree=$2
suffix=${3:-species}
altype=${4:-prot}  # "prot" or "cds"
subtreedir=${5:-subtreesGoodQualO2}

SCRIPTS=$HOME/scripts

[[ "$altype" =~ ^prot$|^cds$ ]] || { echo "$help">&2 ; exit 2; }

print_variables() {
    echo "    genetree=${genetree}
    subtree=${subtree}
    suffix=${suffix}
    altype=${altype}
    subtreedir=${subtreedir}
    SCRIPTS=${SCRIPTS}
    ---
    subtree_input=${subtree_input:-}
    subtree_species=${subtree_species:-}
    al_root=${al_root:-}
    al_suffix=${al_suffix:-}
    al_input=${al_input:-}
    al_phy=${al_phy:-}
    " >&2
}
trap print_variables INT ERR
transform_keys='gene=s/^ENS(G|[A-Z]{3}G)0{5}|^MGP_(SPRET|CAROLI|Pahari)EiJ_G0|^WBGene00/G/'

subtree_input="./$genetree/$subtreedir/${subtree}_codeml.nwk"
subtree_species="./$genetree/$subtreedir/${subtree}_codeml.${suffix}.nwk"
$SCRIPTS/seqtools/specify.py -e 93 -f nwk -l '{shortsp}_{gene}' \
    -t "$transform_keys" "$subtree_input" "$subtree_species"

# Use this option to significantly shorten most sequence names: -t 'gene=s/^ENS(G|[A-Z]{3}G)0{5}/G/'

if [[ "$altype" = 'prot' ]]; then
    al_root="./$genetree/$subtreedir/realign/${subtree}"
    al_suffix="protfsa"
elif [[ "$altype" = 'cds' ]]; then
    al_root="./$genetree/$subtreedir/${subtree}"
    al_suffix="fsa"
fi

al_input="${al_root}_${al_suffix}.fa"
al_phy="${al_root}_${al_suffix}.${suffix}.phy"
# Prefix the species names + convert the AA alignment to phylip
# 'phylip-relaxed' works with pb too.
$SCRIPTS/seqtools/specify.py -e 93 -f fasta -l '{shortsp}_{gene}' \
    -t "$transform_keys" "$al_input" |\
    $SCRIPTS/seqtools/seqconv.py -t 'phylip-sequential-relaxed' > "$al_phy"

echo "Done." >&2
