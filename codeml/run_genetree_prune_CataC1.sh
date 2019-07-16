#!/usr/bin/env bash

set -euo pipefail
IFS=$'\t\n'

while IFS= read -r treefile; do
    $HOME/scripts/dendro/genetree_prune_species.py \
        -e 93 -o "${treefile/subtreesGoodQualO2/subtreesCataC1}" \
        $treefile \
        <(echo -e "Macaca fascicularis\nMacaca nemestrina\nPapio anubis\nMandrillus leucophaeus\nChlorocebus sabaeus\nCercocebus atys\nColobus angolensis palliatus")
 done <&0
