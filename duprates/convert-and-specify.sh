#!/usr/bin/env bash

set -euo pipefail
IFS=$'\t\n'

help="
USAGE

./convert-and-specify.sh <genetree> <subtree> [ <suffix> <alinputtype> <aloutputext> <drop_outgroup> <binarize> <subtreedir>]

ARGUMENTS

  genetree:      e.g. ENSGT00390000000002
  subtree:       e.g. RodentiaENSGT00390000000002.a.a.a.a.a.a
  suffix:        for the output 'specified' alignment [species]
  alinputtype:   'prot'/'cds' [prot]
  aloutputext:   'phy'/'fa' [phy]
  drop_outgroup: 1/0 [0]
  binarize:      1/0 [0]
  subtreedir:    [subtreesGoodQualO2]
"

(( $# > 1 )) || { echo -e "ERROR:Not enough arguments.\n$help">&2 ; exit 2; }

genetree=$1
subtree=$2
suffix=${3:-species}
alinputtype=${4:-prot}  # "prot" or "cds"
aloutputext=${5:-phy}
drop_outgroup=${6:-0}
binarize=${7:-0}
subtreedir=${8:-subtreesGoodQualO2}

SCRIPTS=$HOME/scripts

[[ "$alinputtype" =~ ^prot$|^cds$ ]] || {
    echo -e "ERROR:Wrong alinputtype='${alinputtype}'\n$help">&2 ; exit 2; }
[[ "$aloutputext" =~ ^phy$|^fa$ ]] || {
    echo -e "ERROR:Wrong aloutputext='${aloutputext}'\n$help">&2 ; exit 2; }
[[ "$drop_outgroup" =~ ^0|1$ ]] || {
    echo -e "ERROR:Wrong drop_outgroup='${drop_outgroup}'\n$help">&2 ; exit 2; }
[[ "$binarize" =~ ^0|1$ ]] || {
    echo -e "ERROR:Wrong binarize='${binarize}'\n$help">&2 ; exit 2; }

print_variables() {
    echo "    genetree=${genetree}
    subtree=${subtree}
    suffix=${suffix}
    alinputtype=${alinputtype}
    aloutputext=${aloutputext}
    drop_outgroup=${drop_outgroup}
    binarize=${binarize}
    subtreedir=${subtreedir}
    SCRIPTS=${SCRIPTS}
    ---
    subtree_input=${subtree_input:-}
    subtree_species=${subtree_species:-}
    al_root=${al_root:-}
    al_suffix=${al_suffix:-}
    al_input=${al_input:-}
    al_output=${al_output:-}
    transform_keys=${transform_keys:-}
    outleaves=${outleaves-undefined}
    " >&2
}
trap print_variables INT ERR
# Use this option to significantly shorten most sequence names: -t 'gene=s/^ENS(G|[A-Z]{3}G)0{5}/G/'
transform_keys='gene=s/^ENS(G|[A-Z]{3}G)0{5}|^MGP_(SPRET|CAROLI|Pahari)EiJ_G0|^WBGene00/G/'

subtree_input="./$genetree/$subtreedir/${subtree}_codeml.nwk"
subtree_species="./$genetree/$subtreedir/${subtree}_codeml.${suffix}.nwk"

if (( ${drop_outgroup} )); then
    subtree_input="${subtree_input%_codeml.nwk}.nwk"
    outleaves=$(
        if (( ${binarize} )); then
            $SCRIPTS/duprates/drop_outgroup.py "$subtree_input" \
            | $SCRIPTS/seqtools/specify.py -e 93 -f nwk -l '{shortsp}_{gene}' \
                -t "$transform_keys" - \
            | $SCRIPTS/dendro/unroot_binarise.py -U >"$subtree_species"
        else
            $SCRIPTS/duprates/drop_outgroup.py "$subtree_input" \
            | $SCRIPTS/seqtools/specify.py -e 93 -f nwk -l '{shortsp}_{gene}' \
                -t "$transform_keys" - "$subtree_species"
        fi 2>&1 | sed -n 's/^Outgroup leaves://; s/,/|/gp'
        )
else
    if (( ${binarize} )); then
        $SCRIPTS/seqtools/specify.py -e 93 -f nwk -l '{shortsp}_{gene}' \
            -t "$transform_keys" "$subtree_input" \
        | $SCRIPTS/dendro/unroot_binarise.py -U >"$subtree_species"
    else
        $SCRIPTS/seqtools/specify.py -e 93 -f nwk -l '{shortsp}_{gene}' \
            -t "$transform_keys" "$subtree_input" "$subtree_species"
    fi
fi


if [[ "$alinputtype" = 'prot' ]]; then
    al_root="./$genetree/$subtreedir/realign/${subtree}"
    al_suffix="protfsa"
elif [[ "$alinputtype" = 'cds' ]]; then
    al_root="./$genetree/$subtreedir/${subtree}"
    al_suffix="fsa"
fi

al_input="${al_root}_${al_suffix}.fa"
al_output="${al_root}_${al_suffix}.${suffix}.${aloutputext}"
# Prefix the species names + convert the AA alignment to phylip
# 'phylip-relaxed' works with pb too.
if (( ${drop_outgroup} )); then
    if [[ "aloutputext" = "phy" ]]; then
        $SCRIPTS/seqtools/seqname_grep.py -v "$outleaves" "$al_input" \
        | $SCRIPTS/seqtools/specify.py -e 93 -f fasta -l '{shortsp}_{gene}' \
            -t "$transform_keys" - \
        | $SCRIPTS/seqtools/seqconv.py -t "$al_format_option" > "$al_output"
    else
        $SCRIPTS/seqtools/seqname_grep.py -v "$outleaves" "$al_input" \
        | $SCRIPTS/seqtools/specify.py -e 93 -f fasta -l '{shortsp}_{gene}' \
            -t "$transform_keys" - "$al_output"
    fi
else
    if [[ "aloutputext" = "phy" ]]; then
        $SCRIPTS/seqtools/specify.py -e 93 -f fasta -l '{shortsp}_{gene}' \
            -t "$transform_keys" "$al_input" \
        | $SCRIPTS/seqtools/seqconv.py -t 'phylip-sequential-relaxed' > "$al_output"
    else
        $SCRIPTS/seqtools/specify.py -e 93 -f fasta -l '{shortsp}_{gene}' \
            -t "$transform_keys" "$al_input" "$al_output"
    fi
fi

echo "Done." >&2
