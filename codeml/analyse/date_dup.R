#!/usr/bin/env Rscript

usage="
date_forest.R <datasetname> <calibfile> <treefile> <ncores>
"


library(ape)
library(parallel)


load_calibration <- function(agefile) {
  #ages <- read.delim(agefile, row.names=1)
  #calib <- ages[ages$type == "spe", c("age_dist", "type")]
  #names(calib) <- c("name", "age.min", "type")
  #calib$age.max <- calib$age.min
  # cbind
  # calib$soft.bounds
  # calib$node
  #return(calib)

  # Load without row.names, then drop duplicate row names.
  ages <- read.delim(agefile, row.names=NULL, header=TRUE)
  dup_rownames <- duplicated(ages[,'name'])
  if( any(dup_rownames) ) {
    cat('WARNING: drop', sum(dup_rownames), 'duplicate row.names from ages:\n',
        paste(head(ages[dup_rownames,'name'], 20), collapse='\n '), '...\n', file=stderr())
    ages <- ages[!dup_rownames,]
  }
  rownames(ages) <- ages[,'name']
  ages <- ages[,-1]  # Assuming 'name' is the first column
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

# the MPL algo works only with purely dichotomic trees. So I need to randomly
# convert multifurcations to bifurcations, but then convert them back with the
# correct length
restore_multi_custom <- function(tree, max_number) {
  # Remove all nodes that are above the given number, and transfer its
  # edge length to its descendants.
  ntips <- Ntip(tree)

  e  <- tree$edge
  el <- tree$edge.length
  
  # Update the edge matrix
  #removed_base <- which(tree$edge[,1] > max_number)
  removed_ends <- which(e[,2] > max_number)
  
  # Put the dist into descendants branches and set dist to zero
  for (i in removed_ends) {
    removed_i_base <- which(e[,1] == e[i, 2])
    tree$edge.length[removed_i_base] <- el[removed_i_base] + el[i]
    
    # update the node indexing of the edges (short-circuit the node to delete)
    e[removed_i_base,1] <- e[i,2]
  }
  tree$edge.length[removed_ends] <- 0

  # Update the edge matrix
  tree$edge <- e[-removed_ends,]

  # Delete nodes
  tree$node.label <- tree$node.label[1:(max_number - ntips)]
  return(tree)
}

mark_extra_edges <- function(tree, max_number) {
  e <- tree$edge
  el <- tree$edge.length

  # In case there are negative values
  min_el <- min(0, el)
  
  removed_edges <- which(e[,2] > max_number)
  
  #print(length(removed_edges))
  for(r in removed_edges) {
    # extend the lengths of edges **starting** with that node.
    descendants_of_removed <- which(e[,1] == e[r,2])
    tree$edge.length[descendants_of_removed] <- el[descendants_of_removed] + el[r]
    tree$edge.length[r] <- min_el - 1
  }
  return(tree)
}

restore_multi <- function(multi_tree, dicho_tree) {
  dicho_tree_marked <- mark_extra_edges(dicho_tree,
                                        Ntip(multi_tree) + Nnode(multi_tree))
  return(di2multi(dicho_tree_marked, tol=min(0, dicho_tree$edge.length)))
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

chronos_PL1 <- function(tree, calibration, lambda=1, ...) {
  return(chronos(tree, lambda=lambda, calibration=calibration, ...))
}
chronos_PL100 <- function(tree, calibration, lambda=100, ...) {
  return(chronos(tree, lambda=lambda, calibration=calibration, ...))
}
chronos_PL10000 <- function(tree, calibration, lambda=10000, ...) {
  return(chronos(tree, lambda=lambda, calibration=calibration, ...))
}

chronos_NPRS <- function(tree, calibration, ...) {
  return(chronos(tree, lambda=10^16, calibration=calibration, ...))
}

chronos_MPL <- function(subtree, calibration, ...) {
  # TODO: split the input tree into subtrees with calibrated leaves
  if ( !is.binary(subtree) ) subtree <- multi2di(subtree)
  return(chronoMPL(subtree, ...))
}

date_PL0 <- function(subtree, subcalib) {
  return(chronos(subtree, lambda=0, calibration=subcalib, quiet=TRUE))
}
date_PL1 <- function(subtree, subcalib) {
  return(chronos(subtree, lambda=1, calibration=subcalib, quiet=TRUE))
}
date_PL100 <- function(subtree, subcalib) {
  return(chronos(subtree, lambda=100, calibration=subcalib, quiet=TRUE))
}
date_PL10000 <- function(subtree, subcalib) {
  return(chronos(subtree, lambda=10000, calibration=subcalib, quiet=TRUE))
}
date_NPRS <- function(subtree, subcalib) {
  return(chronos(subtree, lambda=1e16, calibration=subcalib, quiet=TRUE))
}
date_R1 <- function(subtree, subcalib) {
  return(chronos(subtree, model='relaxed', calibration=subcalib, quiet=TRUE))
}
date_R100 <- function(subtree, subcalib) {
  return(chronos(subtree, lambda=100, model='relaxed', calibration=subcalib, quiet=TRUE))
}
date_R10000 <- function(subtree, subcalib) {
  return(chronos(subtree, lambda=10000, model='relaxed', calibration=subcalib, quiet=TRUE))
}
date_C <- function(subtree, subcalib) {
  return(chronos(subtree, model='discrete', calibration=subcalib,
                 control=chronos.control(nb.rate.cat=1), quiet=TRUE))
}

date_MPL <- function(subtree, subcalib) {
  is_multi <- !is.binary(subtree)
  disubtree <- if(is_multi) {multi2di(subtree)} else {subtree}
  chronogram <- chronoMPL(disubtree, se=TRUE, test=TRUE)
  #neg_branchlen <- chronogram$edge.length < 0
  if( is_multi ) {
    chronogram <- restore_multi(subtree, chronogram)
    attr(chronogram, "stderr") <- attr(chronogram, "stderr")[1:Nnode(subtree)]
    attr(chronogram, "Pval")   <- attr(chronogram, "Pval")[1:Nnode(subtree)]
  }
  # And rescale edge lengths
  MPL_depth <- max(node.depth.edgelength(chronogram))
  real_depth <- subcalib[chronogram$node.label[1], "age.min"]
  # If multiple calibration, do that separately for each which.edge(rowname of subcalib)
  chronogram$edge.length <- chronogram$edge.length * real_depth/MPL_depth
  attr(chronogram, "stderr") <- attr(chronogram, "stderr") * real_depth/MPL_depth
  return(chronogram)
}



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

    dtree <- tryCatch(date_func(tree, calibration=subcalib, ...),
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
    if (data.class(dtree) != "chronos") {
      # return the error
      return(c(as.character(as.vector(dtree)), tree$node.label[1], test1, test2, test3))
    }
    root.age <- subcalib[tree$node.label[1], "age.min"]
    
    dated <- chronogram2table(dtree, root.age, leaf_ages[1],
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

extract_MPL_nodeinfo <- function(dtree) {
  # For outputs of the `chronoMPL` function
  nodeinfo <- do.call(cbind, attributes(dtree)[c("stderr", "Pval")])
  rownames(nodeinfo) <- dtree$node.label
  return(nodeinfo)
}

extract_rate_edgeinfo <- function(dtree) {
  # For outputs of the `chronos` function
  edgeinfo <- attr(dtree, "rates")
  # Add terminal node name as edge label
  #if( length(edgeinfo)>1 )
  names(edgeinfo) <- c(dtree$tip.label, dtree$node.label)[dtree$edge[,2]]
  return(edgeinfo)
}

extract_runinfo <- function(dtree) {
  # For outputs of the `chronos` function
  dattr <- attributes(dtree)
  return(list(message=dattr$message,
              ploglik=dattr$ploglik,
              loglik=dattr$PHIIC$logLik,
              PHIIC=dattr$PHIIC$PHIIC))
}

date_methods <- list(PL1=date_PL1, PL100=date_PL100, PL10000=date_PL10000, NPRS=date_NPRS, C=date_C, R1=date_R1, R100=date_R100, R10000=date_R10000, MPL=date_MPL)

date_all_methods <- function(line, calibration,
                         date_methods=list(
                                           PL1=date_PL1,
                                           PL100=date_PL100,
                                           PL10000=date_PL10000,
                                           NPRS=date_NPRS,
                                           C=date_C,
                                           R1=date_R1,
                                           R100=date_R100,
                                           R10000=date_R10000,
                                           MPL=date_MPL),
                         age_col="age_dist") {
  count_iter <<- count_iter + 1
  tree <- read.tree(text=line)
  cat("\n", count_iter, tree$node.label[1], "     ")
  leaf_ages <- calibration[tree$tip.label, age_col]
  # In case tips are not in the table, ignore (assume leaf age is zero)
  leaf_ages[is.na(leaf_ages)] <- 0
  test1 <- all(leaf_ages == leaf_ages[1])
  subcalib <- subset_calibration(calibration, tree, age_col=age_col, leaf_age=leaf_ages[1])
  test2 <- nrow(subcalib) < tree$Nnode  # There are remaining nodes to date.
  test3 <- nrow(subcalib) > 0
  test4 <- all(tree$edge.length > 0)
  if (test1 && test2 && test3 && test4) {

    dtrees <- lapply(date_methods, do.call, list(tree, subcalib))

    dated <- sapply(dtrees, branching.times)
    # Standard error of node ages and pvalue
    dating_nodeinfo <- extract_MPL_nodeinfo(dtrees[["MPL"]])
    # Extend dated with the tip data.
    dating_nodeinfo <- rbind(dating_nodeinfo, matrix(NA, Ntip(tree), ncol(dating_nodeinfo)))

    dated <- rbind(dated, matrix(0, Ntip(tree), ncol(dated)))
    # Extend dated with the tip data.
    rownames(dated)[(Nnode(tree)+1):(Ntip(tree)+Nnode(tree))] <- tree$tip.label
    
    # Molecular rates
    dating_edgeinfo <- sapply(dtrees[c("PL1", "PL100", "PL10000", "NPRS", "R1",
                                       "R100", "R10000")],
                              extract_rate_edgeinfo)
    dating_edgeinfo <- cbind(dating_edgeinfo, C=attr(dtrees[["C"]], "rates"))
    dating_edgeinfo <- rbind(NA, dating_edgeinfo)
    rownames(dating_edgeinfo)[1] <- tree$node.label[1]
    # Reorder based on the nodes then the tips, as in the other data.frames
    dating_edgeinfo <- dating_edgeinfo[c(tree$node.label, tree$tip.label), ]

    info <- calibration[c(tree$node.label, tree$tip.label), c("calibrated", "type", "parent", "taxon", "subgenetree")]
    rownames(info) <- c(tree$node.label, tree$tip.label)
    info[tree$tip.label,"subgenetree"] <- info[tree$node.label[1],"subgenetree"]
    #print(info[tree$node.label[1],"subgenetree"])

    #if ( !all(sapply(list(dating_nodeinfo, dating_edgeinfo, info),
    #                 function(d1, d2) identical(rownames(d1), rownames(d2)),
    #                 dated)
    #         ) )
    #  stop("NOT equal rownames\n")

    dated_info <- data.frame(age=dated, age.MPL=dating_nodeinfo, rate=dating_edgeinfo, info)

    runinfo <- unlist(lapply(dtrees[c("PL1", "PL100", "PL10000", "R1", "R100",
                                      "R10000", "C")],
                             extract_runinfo),
                      recursive=FALSE)

    loginfo <- c(tree$node.label[1], "calibrate", test1, test2, test3, test4)
    return(list(log=loginfo, run=runinfo, ages=dated_info))
  } else {
    loginfo <- c(tree$node.label[1], "skip", test1, test2, test3, test4)
    return(list(log=loginfo, run=NULL, ages=NULL))
  }

}

date_all_trees <- function(allnewicks, calib, outfilename, date_func=date_PL,
                           age_col="age_dist", ...) {
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
  x_c <- date_all_trees(treefile, agefile, outfile, model='discrete',
                  control=chronos.control(nb.rate.cat=1), quiet=TRUE)

  outfile <- paste0(datasetname, "-PL1.tsv")
  x1 <- date_all_trees(treefile, agefile, outfile, quiet=TRUE)

  outfile <- paste0(datasetname, "-PL100.tsv")
  x100 <- date_all_trees(treefile, agefile, outfile, lambda=100, quiet=TRUE)

  outfile <- paste0(datasetname, "-PL10000.tsv")
  x10000 <- date_all_trees(treefile, agefile, outfile, lambda=10000, quiet=TRUE)
  
  # Penalized relaxed
  outfile <- paste0(datasetname, "-PLR.tsv")
  x_r <- date_all_trees(treefile, agefile, outfile, model='relaxed', quiet=TRUE)

  ## MPL
  outfile <- paste0(datasetname, "-MPL.tsv")
  xMPL <- date_all_trees(treefile, agefile, outfile, date_MPL,
                         se=FALSE, test=FALSE)
}


date_all_trees_all_methods <- function(datasetname, agefile, treefile, ncores=6, n=-1) {
  calib <- load_calibration(agefile)
  allnewicks <- scan(treefile, what=character(), n=n, sep="\n")
  outfile_ages <- paste0(datasetname, "_chronos-ages.tsv")
  outfile_runs <- paste0(datasetname, "_chronos-runs.tsv")
  outfile_logs <- paste0(datasetname, "_chronos-logs.tsv")
  count_iter <<- 0

  # Name each tree by its root node name.
  names(allnewicks) <- strcapture("^.*\\)([A-Za-z.0-9]+):[0-9.-][^(),]+;$",
                                  allnewicks,
                                  proto=data.frame(root=character()))$root

  if( ncores == 1 ) {
    out <- lapply(allnewicks, date_all_methods, calib)
  } else {
    cl <- makeForkCluster(ncores, outfile="")
    out <- parLapply(cl, allnewicks, date_all_methods, calib)
    stopCluster(cl)
  }

  out_ages <- do.call(rbind, lapply(out, `[[`, "ages"))  # Do not check row names
  out_runs <- do.call(rbind, lapply(out, `[[`, "run"))
  out_logs <- data.frame(do.call(rbind, lapply(out, `[[`, "log")))
  colnames(out_logs) <- c("name", "action", "test1", "test2", "test3","test4")
  #out_logs$test1 <- as.logical(out_logs$test1)
  #out_logs$test2 <- as.logical(out_logs$test2)
  #out_logs$test3 <- as.logical(out_logs$test3)

  #return(list(out_ages, out_runs, out_logs))
  write.table(out_ages, outfile_ages, quote=FALSE, sep="\t", na="")
  write.table(out_runs, outfile_runs, quote=FALSE, sep="\t", na="")
  write.table(out_logs, outfile_logs, quote=FALSE, sep="\t", na="", row.names=FALSE)
}

if(!interactive()) {
  args <- commandArgs(trailingOnly=TRUE)
  if( length(args) != 5 ) {
    stop("Wrong number of arguments", usage)
  }
  #datasetname <- "Simiiformes_m1w04_ages.subtreesCleanO2-um2-withSG"
  #run(datasetname)
  datasetname <- args[1]
  agefile <- args[2]
  treefile <- args[3]
  ncores <- as.integer(args[4])
  nlines <- ifelse(length(args) >= 5, as.integer(args[5]), -1)

  date_all_trees_all_methods(datasetname, agefile, treefile, ncores, nlines)
}
