#!/usr/bin/env Rscript

# RUN THIS ANALYSIS FROM: ~glouvel/ws7/DUPLI_data93/alignments_analysis/


library(ape)
library(geiger)
library(RPANDA)


treefile_Opistho <- '~/ws7/databases/timetree/Opisthokonta_species.taxa.nwk'

tree_Opistho <- read.tree(treefile_Opistho)
Nopistho <- Ntip(tree_Opistho)

Primate_species <- tips(tree_Opistho, Nopistho + which(tree_Opistho$node.label == "Primates"))
Simii_species <- tips(tree_Opistho, Nopistho + which(tree_Opistho$node.label == "Simiiformes"))

tree_Primates <- extract.clade(tree_Opistho, 'Primates')
tree_Simii <- extract.clade(tree_Opistho, 'Simiiformes')

# From ~/ws7/DUPLI_data93/alignments_analysis/
# fit_ClaDS(tree, samplefraction, iterations, ...)
primates_sampling <- 0.8
# Hedges & Kumar 2015 (Supp table 4) sampling effort for Primates: 0.73.
# I suppose there has been species added since then.

clads.div.result <- fit_ClaDS(tree_Simii, primates_sampling, 500000, nCPU=3, file_name="div_clads_mcmc.txt", it_save=1000, model_id="ClaDS2")
