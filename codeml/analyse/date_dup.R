#!/usr/bin/env Rscript

usage="
date_forest.R <datasetname> <calibfile> <treefile>"


library(ape)

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
  # return a "calibration" dataframe formatted for ape::chronos.
  # Does not retrieve leaf ages.
  ntips <- Ntip(tree)
  subcalib_names <- tree$node.label
  subcalib_nodes <- seq_along(subcalib_names) + ntips
  subcalib_agemin <- ages[subcalib_names, age_col] - leaf_age
  
  # handle NA type nodes (not found in ages: duplication as root)
  calibrated <- ages[subcalib_names, "calibrated"] == 1
  calibrated[is.na(calibrated)] <- TRUE
  max_age <- max(subcalib_agemin, na.rm=TRUE)
  subcalib_agemin[is.na(calibrated)] <- max_age

  subcalib <- data.frame(node=subcalib_nodes,
                         age.min=subcalib_agemin,
                         age.max=subcalib_agemin,
                         soft.bounds=FALSE,
                         row.names=subcalib_names)[calibrated,]
}


# Available dating methods:
#  - chronos: penalized likelihood
#  - chronos(phy, lambda=Inf): NPRS
#  - chronosMPL: Mean Path Length, Britton et al

date_line <- function(line) {
  func_params <- list(...)
  tree <- read.tree(text=line)
}

date_tree <- function(tree, dating_func=chronos, ...) {
  func_params <- list(...)

}

date_PL1 <- function(tree, calibration, lambda=1, ...) {
  return(chronos(tree, lambda=lambda, calibration=calibration, ...))
}
date_PL100 <- function(tree, calibration, lambda=100, ...) {
  return(chronos(tree, lambda=lambda, calibration=calibration, ...))
}
date_PL10000 <- function(tree, calibration, lambda=10000, ...) {
  return(chronos(tree, lambda=lambda, calibration=calibration, ...))
}

date_NPRS <- function(tree, calibration, ...) {
  return(chronos(tree, lambda=1000, calibration=calibration, ...))
}

date_MPL <- function(subtree, calibration, ...) {
  # TODO: split the input tree into subtrees with calibrated leaves
  if ( !is.binary(subtree) ) subtree <- multi2di(subtree)
  return(chronoMPL(subtree, ...))
}

#date_ML <- function(tree, calibration) {
#  return(chronos(tree, model='relaxed', calibration=calibration))
#}


### NOTE: the best lambda value (penalty roughness) must be chosen by
###       cross-validation. Testing multiple subtrees to check the prediction on
###       the removed nodes
#date_bestPL <- function(tree, calibration, minl=0, maxl=2, lstep=0.1) {
#  all_chro <- lapply(seq(minl, maxl, lstep), function(lambda) {
#                        chronos(tree, lambda=lambda, calibration=calibration)})
#  #all_lambda <- sapply(all_chro, function(chro) {attr(chro, "PHIIC", exact=T)$lambda})
#  all_ploglik <- sapply(all_chro, function(chro) {attr(chro, "ploglik", exact=T)})
#
#  return(all_chro[which.min(all_ploglik)])
#
#}

# Function that fills an age at one row (i.e one edge). To use in `apply`.
fill.dates <- function(edgedatum, datation.env=new.env(), scale=1) {
  edge    <- as.character(edgedatum[1:2])
  elength <- edgedatum[3]
  datation.env$dated[edge[2], "age"] <- datation.env$dated[edge[1], "age"] - scale * elength
}

chronogram2table <- function(chronogram,
                             root.age,
                             leaf.age=0,
                             datation.env=new.env(), 
                             rescale=FALSE)
{
  # also see ape::branching.times
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
  scale <- ifelse(rescale, root.age - leaf.age, 1)
  edge_data <- edge_data[edge_data[,2] > ntips,,drop=FALSE]
  apply(edge_data, 1, fill.dates, datation.env, scale)
  return(datation.env$dated)
}


process_line <- function(line, calibration, outfile,
                         date_func=date_PL,
                         #date_func_list=list(
                         age_col="age_dist", datation.env=new.env(),
                         col.names=TRUE, ...) {
  count_iter <<- ifelse(exists("count_iter"), count_iter + 1, 0)
  rescale <- identical(date_func, date_MPL)
  tree <- read.tree(text=line)
  leaf_ages <- calibration[tree$tip.label, age_col]
  # In case tips are not in the table, ignore (assume leaf age is zero)
  leaf_ages[is.na(leaf_ages)] <- 0
  test1 <- all(leaf_ages == leaf_ages[1])
  subcalib <- subset_calibration(calibration, tree, age_col=age_col, leaf_age=leaf_ages[1])
  test2 <- nrow(subcalib) < tree$Nnode  # There are remaining nodes to date.
  test3 <- nrow(subcalib) > 0
  test4 <- all(tree$edge.length > 0)
  cat("\r", count_iter, tree$node.label[1], "     ")
  if (test1 && test2 && test3 && test4) {
    chronogram <- tryCatch(date_func(tree, calibration=subcalib, ...),
                            error=function(e) {
                                    cat(" Got Error:")
                                    if (e$message == "NA/NaN gradient evaluation") {
                                      cat("known error", e$message)
                                      warning(e)
                                      return(e)
                                    } else {
                                      cat("unknown")
                                      stop(e)
                                    }})
    if (data.class(chronogram) != "chronos") {
      # return the error
      return(c(as.character(as.vector(chronogram)), tree$node.label[1], test1, test2, test3))
    }
    root.age <- subcalib[tree$node.label[1], "age.min"]
    
    dated <- chronogram2table(chronogram, root.age, leaf_ages[1],
                              datation.env=datation.env,
                              rescale=rescale)
    dated <- cbind(dated,
                   calibration[as.character(dated$name),
                               c("calibrated", "type", "parent", "taxon", "subgenetree")])
    write.table(dated, outfile, sep='\t', quote=FALSE,
                row.names=FALSE, col.names=col.names)
    return(c("calibrate", tree$node.label[1], test1, test2, test3))
  } else {
    return(c("skip", tree$node.label[1], test1, test2, test3))
  }

}


#date_all <- function(treefile, agefile, outfilename, date_func=date_PL,
#                     age_col="age_dist", n=-1, ...) {
#  require(ape)
#  calib <- load_calibration(agefile)
#  #alltrees <- read.tree(treefile)
#  allnewick <- scan(treefile, what=character(), n=n, sep="\n")
#}

date_all <- function(allnewicks, calib, outfilename, date_func=date_PL,
                     age_col="age_dist", n=-1, ...) {
  datation.env <- new.env()
  count_iter  <<- 0
  outfile <- file(outfilename, "w")
  cat("name\tage\tcalibrated\ttype\tparent\ttaxon\tsubgenetree\n", file=outfile)
  #cat("\n")
  r <- sapply(allnewicks, process_line, calib, outfile, date_func, age_col,
              datation.env, col.names=FALSE, ..., USE.NAMES=FALSE)
  close(outfile)
  skipped <- sum(r[1,] == "skip")
  calibrated <- sum(r[1,] == "calibrate")
  failed <- ncol(r) - skipped - calibrated
  cat("Calibrated:", calibrated, "Skipped:", skipped, "Failed:", failed, "\n")
  return(t(r))
}


run <- function(datasetname, agefile, treefile) {
                #outdir="/users/ldog/glouvel/ws2/DUPLI_data85/alignments_analysis/ages",
                #subtrees_src=".dSsubtrees.nwk") {
  ## Penalized likelihood
  #agefile <- "ages/Rodentia_M1_agesdNdSdist.subtreesNoEdit-um2.tsv"
  #treefile <- "trees/Rodentia_M1_dist.subtreesNoEdit.fulltree.nwk"
  #datasetbase <- paste0(outdir, "/", dataset)
  #agefile <- paste0(datasetbase, ".tsv")
  #treefile <- paste0(outdir_trees, datasetbase, subtrees_src)

  calib <- load_calibration(agefile)
  allnewicks <- scan(treefile, what=character(), n=n, sep="\n")

  outfile <- paste0(datasetname, "-C.tsv")
  x_c <- date_all(treefile, agefile, outfile, model='discrete',
                  control=chronos.control(nb.rate.cat=1), quiet=TRUE)

  outfile <- paste0(datasetname, "-PL1.tsv")
  x1 <- date_all(treefile, agefile, outfile, quiet=TRUE)

  outfile <- paste0(datasetname, "-PL100.tsv")
  x100 <- date_all(treefile, agefile, outfile, lambda=100, quiet=TRUE)

  outfile <- paste0(datasetname, "-PL10000.tsv")
  x10000 <- date_all(treefile, agefile, outfile, lambda=10000, quiet=TRUE)
  
  # Penalized relaxed
  outfile <- paste0(datasetname, "-PLR.tsv")
  x_r <- date_all(treefile, agefile, outfile, model='relaxed', quiet=TRUE)

  ## MPL
  outfile <- paste0(datasetname, "-MPL.tsv")
  xMPL <- date_all(treefile, agefile, outfile, date_MPL,
                   se=FALSE, test=FALSE)
}

if(!interactive()) {
  args <- commandArgs(trailingOnly=TRUE)
  if( length(args) != 3 ) {
    cat(usage)
    stop()
  }
  dataset <- "Simiiformes_m1w04_ages.subtreesCleanO2-um2-withSG"
  run(dataset)
}
