#!/usr/bin/env bash

set -euo pipefail
IFS=$'\t\n'


echo '[FAMILIES]'

i=0 || true
while read -r genetree subtree; do
    #echo "- $genetree - $subtree"
    if (( i < 6000 )); then suffix='shortsp'; else suffix='species'; fi
    (( ++i ))
    echo -ne "\r% ${i}" >&2
    alfile="${genetree}/subtreesGoodQualO2/realign/${subtree}_protfsa.${suffix}.fa"
    [[ -f "$alfile" ]] || $HOME/scripts/seqtools/seqconv.py -f phylip-relaxed -t fasta "${alfile%.fa}.phy" "$alfile"
    echo -e "- ${subtree}
starting_gene_tree = ${genetree}/subtreesGoodQualO2/${subtree}_codeml.${suffix}.nwk
alignment = ${alfile}
subst_model = LG+G"
done < "$1"
