#!/usr/bin/env bash

set -euo pipefail
IFS=$'\t\n'

while IFS= read -r treefile; do
    $HOME/scripts/dendro/genetree_prune_species.py \
        -e 93 -o "${treefile/subtreesGoodQualO2/subtreesCataH1}" \
        $treefile \
        <(echo -e "Pan troglodytes\nPan paniscus\nGorilla gorilla gorilla\nPongo abelii\nNomascus leucogenys")
 done <&0
