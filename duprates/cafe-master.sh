#!/usr/bin/env bash

set -euo pipefail
IFS=$'\t\n'

SCRIPTS=$HOME/scripts
CAFEDIR=$HOME/install/CAFE5/
CAFEBIN=$CAFEDIR/bin
subtreesdir="subtreesGoodQualO2"

# Step 1: Generate the tabulated file of gene counts per family and species.

cd $HOME/ws7/DUPLI_data93/alignments/
species_tree=../PhylTree.TimeTree201901.Ensembl93-like.goodQual.nwk
family_file=../alignments_analysis/subtrees_stats/"$subtreesdir"_familystats-Simiiformes.tsv

#while IFS= read -r family subtree; do
#    treefile="${family}/${subtreesdir}/${subtree}.nwk"
#    #$SCRIPTS/genomicustools/identify.py -e93 -f "txt" $treefile | sort | uniq -c
#done <<< all_Simiiformes.subtreesGoodQualO2.Queue.txt

# This was slow: I used splitjob.sh.
$SCRIPTS/codeml/subtrees_stats.py family -s "$subtreesdir" -I -e 93 -i \
    all_gene_trees.shuffled.txt Simiiformes "$species_tree" \
    > "$family_file" \
    2> "${family_file%.tsv}.log"
# Replace spaces in species names by underscores
sed -i '1s/ /_/g' "$family_file"

simii_tree=../taxon_data/PhylTree-Simii.TimeTree201901.Ensembl93-like.goodQual.binary.nonull.basic.nwk

# Check that species tree is valid

Rscript --vanilla $SCRIPTS/duprates/check_newick.R "$simii_tree"

# Check that the set of leaf names is *exactly* identical:

comm -3 <(nw_labels -I "$simii_tree" | sort) \
    <(head -1 subtrees_stats/subtreesGoodQualO2_familystats-Simiiformes.tsv | cut -f3- | tr '\t' '\n' | sort)

# From CAFE tutorial, it is recommended to filter lines with >100 genes in one species.

cd subtrees_stats/
python3 $CAFEDIR/tutorial/clade_and_size_filter.py -s -i "$family_file" -o "${family_file%.tsv}.filtered.tsv"


