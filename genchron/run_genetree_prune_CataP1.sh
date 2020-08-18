#!/usr/bin/env bash

set -euo pipefail
IFS=$'\t\n'

while IFS= read -r treefile; do
    $HOME/scripts/dendro/genetree_prune_species.py \
        -e 93 -o "${treefile/subtreesGoodQualO2/subtreesCataP1}" \
        $treefile \
        <(echo -e "Saimiri boliviensis boliviensis\nAotus nancymaae\nCebus capucinus imitator")
 done <&0
