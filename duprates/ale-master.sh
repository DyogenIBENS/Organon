#!/usr/bin/env bash

# Just a backbone/README for the analysis. Everything was run in parallel with Condor.
# See:
# - fa2phy-for-pb_Simiiformes.subtreesGoodQualO2.condor.txt
# - pb_Simiiformes.subtreesGoodQualO2.condor.txt


set -euo pipefail
IFS=$'\t\n'

SCRIPTS=$HOME/scripts
PBDIR=$HOME/install/phylobayes4.1c/data
ALEDIR=$HOME/install/ALE/build/bin
NCORES=4

# cd ~/ws7/DUPLI_data93/alignments/
# tail all_Simiiformes-robusts_fsa.subtreesGoodQualO2.txt 
species_tree="../PhylTree.TimeTree201901-Ensembl93-like.goodQual.basic.nwk"

subtreedir='subtreesGoodQualO2'

# Example data.
genetree='ENSGT00390000001638'
subtree='SimiiformesENSGT00390000001638'

subtree_input="./$genetree/$subtreedir/${subtree}_codeml.nwk"
subtree_species="./$genetree/$subtreedir/${subtree}_codeml.species.nwk"
al_root="./$genetree/$subtreedir/realign/${subtree}"
al_input="${al_root}_protfsa.fa"
al_phy="${al_root}_protfsa.species.phy"
dataroot="${al_root}_pb"

# If needed, translate because we will align aa.
#$SCRIPTS/seqtools/fasta_translate.py "$subtree_phy" "${subtree_prot_phy}"

$SCRIPTS/seqtools/specify.py -e 93 -f nwk -l '{sp}_{gene}' "$subtree_input" "$subtree_species"

# Prefix the species names + convert the AA alignment to phylip
# 'phylip-relaxed' works with pb too.
$SCRIPTS/seqtools/specify.py -e 93 -f fasta -l '{sp}_{gene}' "$al_input" |\
    $SCRIPTS/seqtools/seqconv.py -t 'phylip-sequential-relaxed' > "$al_phy"


#mpirun -np "$NCORES" pbmpi -d "$al_phy" 
# pb options
# -x <every> [<until>]
# -dgam <n> : Number of categories for discrete gamma.
# -ratecat : CAT model
# -lg -wag -jtt -mtrev -mtzoa -mtart : empirical exchangeabilities
$PBDIR/pb -x 1 2000 -dgam 1 -lg -t "$subtree_species" -d "$al_phy" "${dataroot}1" &
mcmc_pid1=$!
$PBDIR/pb -x 1 2000 -dgam 1 -lg -t "$subtree_species" -d "$al_phy" "${dataroot}2" &
mcmc_pid2=$!
# TODO: Run 2 parallel chains and estimate their discrepancy with bpcomp and tracecomp
wait $mcmc_pid1 $mcmc_pid2

mcmctrees="${dataroot}.treelist"
# Check out output: Tracer.
aletrees="${al_root}.treelist"

alefile="${aletrees}.ale"

# Autocorrelation time is typically 2 or 3 for the posterior proba.
# Discard burnin. 300 should be enough.
#bp > "$aletrees"
sed -n '300~2p' "$mcmctrees" > "$aletrees"

# Prepare the forest of trees:
$ALEDIR/ALEobserve "$aletrees"
# Or directly
#ALEobserve "$mcmctrees" 300

# ALEmcmc exist only for the "undated" algo
# The other is called ALEsample.
#
# Arguments of ALEml:
# burnin=0 is invalid.
# tau=0
# output_species_tree=<y/n>
# fraction_missing=<file with fraction of missing genes per species>
#
# The output file is saved in the species tree folder... --'
$ALEDIR/ALEsample "$species_tree" "$alefile"

$ALEDIR/ALEml "$species_tree" "$alefile" tau=0 sample=200

#ALEml failed with:
#
#    terminate called after throwing an instance of 'std::out_of_range'
#      what():  vector::_M_range_check: __n (which is 18446744073709551615) >= this->size() (which is 512)
#    Abandon (core dumped)

# For data:
# - ENSGT00390000001638/subtreesGoodQualO2/realign/SimiiformesENSGT00390000001638_pb.treelist.shortnames.ale
# 1. Test shorter leaf names without dots;
# 2. test by putting only trees without null branch lengths (`head -3`/`tail -2`)
# Still fails.
