#!/usr/bin/env Rscript

library(ape)

agefile <- "ages/Rodentia_M1_agesdNdSdist.subtreesNoEdit-um2.tsv"
treefile <- "trees/Rodentia_M1_dist.subtreesNoEdit.fulltree.nwk"

load_calibration  <- function(agefile) {
    ages <- read.delim(agefile)
    calib <- ages[ages$type != "dup", c("name", "age_dist", "type")]
    names(calib) <- c("name", "age.min", "type")
    calib$age.max <- calib$age.min
    # cbind
    # calib$soft.bounds
    # calib$node
    return(calib)
}

subset_calibration <- function(calib, tree) {
  ntips <- length(tree$tip.label)
  subcalib_names <- tree$node.label
  subcalib_nodes <- sapply(tree$node.label, function (nl) {grep(nl, subcalib_names) + ntips})
  subcalib_agemin <- sapply(subcalib_names, function(nl) {
                              calib$age.min[calib$name == nl]})
  subcalib <- data.frame(node=subcalib_nodes,
                         age.min=subcalib_agemin,
                         age.max=subcalib_agemin,
                         soft.bounds=FALSE)
}


date_line <- function(line, calib) {
	tree <- read.tree(text=line)
}


date_all <- function(treefile, agefile) {
  calib <- load_calibration(agefile)
	f <- file(treefile, "r")
	while ( TRUE ) {
		line <- readLines(f, n=1)
		if ( length(line) == 0 ) {
			break
		}
		date_line(line, calib)
	}
	close(f)
}
