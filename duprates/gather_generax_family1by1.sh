#!/usr/bin/env bash

set -euo pipefail
IFS=$'\t\n'

# Collect DTL rates into a single file.

generaxdir="generax_family-1by1"
for (( i=0; i<24562; i++ )); do
    printf -v nb "%05d" $i
    subdir="$generaxdir/${nb::2}/${nb:2}"
    famname="$(ls -U $subdir/results/)" || continue
    echo -ne "${nb}\t${famname}\t"
    cat "$subdir/gene_optimization_5/dtl_rates.txt" \
        "$subdir/results/$famname/stats.txt" \
        || echo -ne '\t\t\t'
    #cut -d ':' -f2 "$subdir/reconciliations/${famname}_eventCounts.txt" | tr '\n' '\t' \
    sed -nr 's/^.*:/\t/; H; ${x;s/\n//g;p}' "$subdir/reconciliations/${famname}_eventCounts.txt" \
        || echo -e '\t\t\t\t\t\t\t'
done | sed 's/ \+/\t/g'
