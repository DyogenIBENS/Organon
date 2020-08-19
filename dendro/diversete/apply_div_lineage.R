#!/usr/bin/env Rscript

# RUN THIS ANALYSIS FROM: ~/ws7/DUPLI_data93/alignments_analysis/


library(ape)
library(geiger)
library(RPANDA)
library(coda)  # Analysis of MCMC


treefile_Opistho <- '~/ws7/databases/timetree/Opisthokonta_species.taxa.nwk'
ClaDS_filename <- 'div_clads_mcmc.RData'  #TODO: .RData

run <- function() {
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

  #clads.div.sampler <- fit_ClaDS(tree_Simii, primates_sampling, 500000, nCPU=3, file_name=ClaDS_filename, it_save=1000, model_id="ClaDS2")
  prev.run <- new.env()
  load('div_clads_mcmc1.RData', prev.run)
  clads.div.sampler <- fit_ClaDS(tree_Simii, primates_sampling, 500000, nCPU=3, file_name=ClaDS_filename,
                                 it_save=1000, model_id="ClaDS2", mcmcSampler=prev.run$mcmcSampler)

  ## From RPANDA doc ##
  ## # 0.85
  ## t_ClaDS_chains(Caprimulgidae_ClaDS2$sampler)
  ## # extract the Maxima A Posteriori for each parameter
  ## maps = getMAPS_ClaDS(Caprimulgidae_ClaDS2$sampler, thin = 1)
  ## print(paste0("sigma = ", maps[1], " ; alpha = ",
  ## maps[2], " ; epsilon = ", maps[3], " ; l_0 = ", maps[4] ))
  ## # plot the infered branch specific speciation rates
  ## plot_ClaDS_phylo(Caprimulgidae_ClaDS2$tree, maps[-(1:4)])

}

check_running <- function() {
  check.env <- new.env()
  load(ClaDS_filename, check.env)
  nparams <- ncol(check.env$mcmcSampler$chains[[1]])
  # clads.div.sampler has a $chains attribute (mcmc.list object)
  cat('\nChains:', check.env$mcmcSampler$Nchain,
      ' Parameters:', nparams,
      ' Thinning:', check.env$mcmcSampler$thin,
      ' Saved states:', nrow(check.env$mcmcSampler$chains[[1]]),
      ' Iterations:', nrow(check.env$mcmcSampler$chains[[1]])*check.env$mcmcSampler$thin,
      '\n')
  ESS <- effectiveSize(check.env$mcmcSampler$chains)
  cat('check\n')
  cat(sprintf('\nESS (%s elements):\n',
              length(ESS)))#, data.class(ESS), dim(ESS)))
  cat('check2\n')
  print(ESS[c(1:5, length(ESS)-1, length(ESS))])
  cat(sprintf('min ESS = %s (%s)\n', min(ESS), names(ESS)[which.min(ESS)]))
  cat(sprintf('max ESS = %s (%s)\n', max(ESS), names(ESS)[which.max(ESS)]))

  X11()
  #dev.new()
  plot_ClaDS_chains(check.env$mcmcSampler)
  cat("\nClic on the window to close it.")
  locator(1) # Wait for the user to click
  cat("\n")
}

# After 65000 iterations:
#    sigma     alpha        mu
# 9.949311 10.186789 47.354408

if(!interactive()) check_running()

