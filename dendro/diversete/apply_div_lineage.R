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

## From RPANDA doc ##
## # 0.85
## t_ClaDS_chains(Caprimulgidae_ClaDS2$sampler)
## # extract the Maxima A Posteriori for each parameter
## maps = getMAPS_ClaDS(Caprimulgidae_ClaDS2$sampler, thin = 1)
## print(paste0("sigma = ", maps[1], " ; alpha = ",
## maps[2], " ; epsilon = ", maps[3], " ; l_0 = ", maps[4] ))
## # plot the infered branch specific speciation rates
## plot_ClaDS_phylo(Caprimulgidae_ClaDS2$tree, maps[-(1:4)])

# clads.div.result has a $chains attribute (mcmc.list object)
claDS_plot_chains(clads.div.result)

library(coda)  # Analysis of MCMC
effectiveSize(clads.div.result$chains)

# After 65000 iterations:
#    sigma     alpha        mu
# 9.949311 10.186789 47.354408
