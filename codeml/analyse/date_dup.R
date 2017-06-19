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
  subcalib_nodes <- sapply(tree$node.label, function (nl) {
                           grep(nl, subcalib_names) + ntips})
  subcalib_agemin <- sapply(subcalib_names, function(nl) {
                              calib$age.min[calib$name == nl]})
  subcalib <- data.frame(node=subcalib_nodes,
                         age.min=subcalib_agemin,
                         age.max=subcalib_agemin,
                         soft.bounds=FALSE)
}


# Available dating methods:
#  - chronos: penalized likelihood
#  - chronos(phy, lambda=0): NPRS
#  - chronosMPL: Mean Path Length, Britton et al

date_line <- function(line) {
  func_params <- list(...)
  tree <- read.tree(text=line)
}

date_tree <- function(tree, dating_func=chronos, ...) {
  func_params <- list(...)

}

date_PL <- function(tree, calibration) {
  return(chronos(tree, calibration=calibration))
}

date_NPRS <- function(tree, calibration) {
  return(chronos(tree, lambda=0, calibration=calibration))
}

date_MPL <- function(tree, calibration) {
  return(chronos(tree, lambda=0, calibration=calibration))
}

date_PL <- function(tree, calibration) {
  return(chronos(tree, model='relaxed', calibration=calibration))
}


# Function that fills an age at one row (i.e one edge). To use in `apply`.
fill.dates <- function(edgedatum, datation.env=new.env(), scale=1) {
  edge    <- as.character(edgedatum[1:2])
  elength <- edgedatum[3]
  datation.env$dated[edge[2], "age"] <- datation.env$dated[edge[1], "age"] - scale * elength
}

chronogram2table <- function(chronogram,
                             root.age,
                             leaf.age=0,
                             datation.env=new.env())
{
  ntips <- length(chronogram$tip.label)
  rootnb <- ntips + 1
  nnodes <- chronogram$Nnode
  tot_nodes <- ntips + nnodes
  # either write calibrated nodes too, or remove them from the edges beforehand
  #dated <- data.frame(row.names=seq(1, ntips+nnodes),
  #                    name=c(chronogram$tip.label, chronogram$node.label), 
  #                    age=numeric(ntips+nnodes))
  dated <- data.frame(row.names=seq(rootnb, tot_nodes),
                      name=chronogram$node.label, 
                      age=numeric(nnodes))
  dated[as.character(rootnb), "age"] <- root.age
  
  assign("dated", dated, envir=datation.env)
  edge_data <- cbind(chronogram$edge, chronogram$edge.length)
  edge_data <- edge_data[edge_data[,2] > ntips,]

  apply(edge_data, 1, fill.dates, datation.env, scale=(root.age - leaf.age))
  return(datation.env$dated)
}


date_all <- function(treefile, agefile) {
  calib <- load_calibration(agefile)
  f <- file(treefile, "r")
  while ( TRUE ) {
    line <- readLines(f, n=1)
    if ( length(line) == 0 ) {
      break
    }
    tree <- read.tree(text=line)
    subcalib <- subset_calibration(calib, tree)
    if (nrow(subcalib) < tree$Nnode) {
      date_tree(line, calibration=tree)
    }
  }
  close(f)
}

