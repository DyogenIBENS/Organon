#!/usr/bin/env Rscript

library(ape)

agefile <- "ages/Rodentia_M1_agesdNdSdist.subtreesNoEdit-um2.tsv"
treefile <- "trees/Rodentia_M1_dist.subtreesNoEdit.fulltree.nwk"

load_calibration <- function(agefile) {
  ages <- read.delim(agefile, row.names=1)
  #calib <- ages[ages$type == "spe", c("age_dist", "type")]
  #names(calib) <- c("name", "age.min", "type")
  #calib$age.max <- calib$age.min
  # cbind
  # calib$soft.bounds
  # calib$node
  #return(calib)
  return(ages)
}

subset_calibration <- function(ages, tree, age_col="age_dist", leaf_age=0) {
  ntips <- Ntip(tree)
  subcalib_names <- tree$node.label
  subcalib_nodes <- seq_along(subcalib_names) + ntips
  subcalib_agemin <- ages[subcalib_names, age_col] - leaf_age
  subcalib_types <- ages[subcalib_names, "type"]
  subcalib <- data.frame(node=subcalib_nodes,
                         age.min=subcalib_agemin,
                         age.max=subcalib_agemin,
                         soft.bounds=FALSE,
                         row.names=subcalib_names)[subcalib_types == "spe",]
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

date_PL <- function(tree, calibration, lambda=1) {
  return(chronos(tree, lambda=lambda, calibration=calibration))
}

date_NPRS <- function(tree, calibration) {
  return(chronos(tree, lambda=0, calibration=calibration))
}

date_MPL <- function(tree, calibration) {
  return(chronos(tree, lambda=0, calibration=calibration))
}

date_ML <- function(tree, calibration) {
  return(chronos(tree, model='relaxed', calibration=calibration))
}

date_bestMPL <- function(tree, calibration, minl=0, maxl=2, lstep=0.1) {
  all_chro <- lapply(seq(minl, maxl, lstep), function(lambda) {
                        chronos(tree, lambda=lambda, calibration=calibration)})
  #all_lambda <- sapply(all_chro, function(chro) {attr(chro, "PHIIC", exact=T)$lambda})
  all_ploglik <- sapply(all_chro, function(chro) {attr(chro, "ploglik", exact=T)})

  return(all_chro[which.min(all_ploglik)])

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
  ntips <- Ntip(chronogram)
  rootnb <- ntips + 1
  nnodes <- Nnode(chronogram)
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

  apply(edge_data, 1, fill.dates, datation.env) #, scale=(root.age - leaf.age))
  return(datation.env$dated)
}


#date_all <- function(treefile, agefile) {
#  calib <- load_calibration(agefile)
#  f <- file(treefile, "r")
#  while ( TRUE ) {
#    line <- readLines(f, n=1)
#    if ( length(line) == 0 ) {
#      break
#    }
#    tree <- read.tree(text=line)
#    subcalib <- subset_calibration(calib, tree)
#    if (nrow(subcalib) < tree$Nnode) {
#      date_tree(line, calibration=tree)
#    }
#  }
#  close(f)
#}

process_line <- function(line, calibration, outfile, date_func=date_PL,
                         age_col="age_dist", datation.env=new.env()) {
  tree <- read.tree(text=line)
  leaf_ages <- calibration[tree$tip.label, age_col]
  test1 <- all(leaf_ages == leaf_ages[1])
  subcalib <- subset_calibration(calibration, tree, age_col=age_col, leaf_age=leaf_ages[1])
  test2 <- nrow(subcalib) < tree$Nnode
  if (test1 & test2) {
    return(c("calibrate", tree$node.label[1], test1, test2))
    chronogram <- date_func(tree, calibration=subcalib)
    root.age <- subcalib[tree$node.label[1], "age"]
    dated <- chronogram2table(chronogram, root.age, datation.env=datation.env)
    dated <- cbind(dated, ifelse(dated$name %in% subcalib), "spe", "dup")
    write.table(dated, outfile, append=TRUE, sep='\t')
  } else {
    return(c("skip", tree$node.label[1], test1, test2))
  }

}


date_all <- function(treefile, agefile, outfile, date_func=date_PL,
                     age_col="age_dist", test=-1) {
  calib <- load_calibration(agefile)
  #alltrees <- read.tree(treefile)
  allnewick <- scan(treefile, what=character(), n=test)
  datation.env <- new.env()
  sapply(allnewick, process_line, calib, outfile, date_func, age_col, datation.env,
         USE.NAMES=FALSE)
}

