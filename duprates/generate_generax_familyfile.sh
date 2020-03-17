#!/usr/bin/env bash

set -euo pipefail
IFS=$'\t\n'

altype=${2:-prot}
suffix='shortsp'

[[ "$altype" =~ ^prot$|^cds$ ]] || { echo 'altype != "prot|cds"' >&2; exit 2; }
if [[ "$altype" = 'prot' ]]; then
    subdir='subtreesGoodQualO2/realign'
    al_suffix='protfsa'
    subst_model='LG+G'
else
    subdir='subtreesGoodQualO2'
    al_suffix='fsa'
    subst_model='GTR+G'
fi

echo '[FAMILIES]'

i=0 || true
while read -r genetree subtree; do
    #echo "- $genetree - $subtree"
    #if (( i < 6000 )); then suffix='shortsp'; else suffix='species'; fi
    (( ++i ))
    echo -ne "\r% ${i}" >&2
    alfile="${genetree}/${subdir}/${subtree}_${al_suffix}.${suffix}.fa"
    [[ -f "$alfile" ]] || $HOME/scripts/seqtools/seqconv.py -f phylip-relaxed -t fasta "${alfile%.fa}.phy" "$alfile"
    echo -e "- ${subtree}
starting_gene_tree = ${genetree}/subtreesGoodQualO2/${subtree}_codeml.${suffix}.nwk
alignment = ${alfile}
subst_model = ${subst_model}"
done < "$1"
