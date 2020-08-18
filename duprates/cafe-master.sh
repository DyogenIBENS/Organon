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
$SCRIPTS/genchron/subtrees_stats.py family -s "$subtreesdir" -I -e 93 -i \
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

cd ~/ws7/DUPLI_data93/alignments_analysis/
# Edit Condor submission files:
# - cafe_onelambda.condor.txt
# - cafe_onelambda-zero_root.condor.txt
# - cafe_lambda-gammadist4.condor.txt
# - cafe_lambda-per-lineage.condor.txt  -> output to cafe/Simiiformes_lambda-per-lineage/

# Run with condor (lambda-per-lineage takes max 4 hours)

# Plot number of families under expansion/contraction

# This script is not in my CAFE directory...
#$CAFEDIR/python_scripts/cafetutorial_draw_tree.py -i reports/summary_run1_node.txt -t <tree> -d <?> -o reports/summary_run1_tree_rapid.png -y Rapid

# The 34 lambdas are in `cafe/Simiiformes_lambda-per-lineage/Base_results.txt`, line 2.

# Now, correlate those with the diversification rates.
