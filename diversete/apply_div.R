#!/usr/bin/env Rscript

library(ape)
library(geiger)
library(RPANDA)
library(diversitree)

if( interactive() ) {
  subtreefile <- "/users/ldog/glouvel/ws2/databases/timetree/Opisthokonta-listens90size10.subtrees.nwk"
  tablefile <- "/users/ldog/glouvel/ws2/databases/timetree/Opisthokonta-listens90size10.tsv"
  args <- c(subtreefile, tablefile)
} else {
  args <- commandArgs(trailingOnly=TRUE)
  subtreefile <- args[1]
}

if(length(args)>1) {
  tablefile  <- args[2]
  clade.data <- read.delim(tablefile, row.names=1, header=TRUE)
} else {
  #tablefile <- NULL
  clade.data <- NULL
}

fix_single_nodes <- function(lines){
  # newick line with a single node and no parenthesis can't be read by read.tree
}

leaf_heights <- function(tree) {
  cum_edge_length <- tree$edge.length
  for(e in seq_along(cum_edge_length)){
      edge <- tree$edge[e,]
      cur_dist <- tree$edge.length[e]
      prev_edge_pos <- tree$edge[,2] == edge[1]
      prev_height <- ifelse(any(prev_edge_pos),
                            cum_edge_length[prev_edge_pos],
                            0)
      cum_edge_length[e] <- prev_height + cur_dist
  }
  leaf_edge_pos <- which(tree$edge[,2] %in% 1:Ntip(tree))
  return(cum_edge_length[leaf_edge_pos])
  # similar to phylobase::nodeHeights:
  #return(phylobase::nodeHeights(tree)[tree$edge[,2] %in% 1:Ntip(tree), 2])
}

extend_to_ultrametric <- function(tree) {
  # Extend tree leaves to the maximum of the leaf-to-root distance.
  # Return the maximum extension performed
  lheights <- leaf_heights(tree)
  maxh <- max(lheights)
  cat(paste0("Max extension: ", maxh - min(lheights), "\n"), file=stderr())
  leaf_edge_pos <- which(tree$edge[,2] %in% 1:Ntip(tree))
  tree$edge.length[leaf_edge_pos] <- tree$edge.length[leaf_edge_pos] + (maxh - lheights)
  return(tree)
}

summary_posterior <- function(samples, paramname, probaname="p",
                              CIrange=c(0.05, 0.95), burnin_prop=0.1){
  # Computes the cumulative distribution function of the parameter, given a
  # mcmc run ('samples'). 'samples' must contain the columns 'paramname'
  # (parameter) and 'probaname' (the posterior probability of the parameter)
  # Return the value of maximum posterior proba and the 5% interval.
  samples <- samples[-(1:floor(nrow(samples) * burnin_prop)),]
  #max_posterior_pos <- which(samples[,probaname] == max(samples[,probaname]))
  #o_samples <- samples[order(samples[,paramname]),]
  
  #tot_area <- sum(diff(o_samples[,paramname]) * o_samples[-1, probaname])
  ## cumulative distribution function
  # cdf <- cumsum(diff(o_samples[,paramname]) * o_samples[-1, probaname]) / tot_area
  #CI_lower_pos <- which(cdf >= 0.05)[1]
  #CI_upper_pos <- which(cdf >= 0.95)[1]
  #CI <- o_samples[c(CI_lower_pos, CI_upper_pos)+1, paramname]
  #best <- samples[max_posterior_pos, paramname]
  post_density <- density(samples[,paramname])
  best <- post_density$x[which.max(dlam$y)]
  CI <- quantile(samples[,paramname], CIrange)
  return(list(best=best, CI=CI))
}


get_div_stats <- function(tree, metadata) {
  stderrf <- stderr()
  if(!is.ultrametric(tree)) tree <- extend_to_ultrametric(tree)
  treename <- tree$node.label[1]
  cat("- Processing", treename, file=stderrf)
  if(Ntip(tree) <= 2) {
    cat("Ntip <= 2. Skip.\n", file=stderrf)
    return()
  }
  sampling_f <- metadata[treename, "ncbi_sp_sampling"]

  out <- list()
  # In Ape (maximum likelihood)
  cat(" Ape... ", file=stderrf)
  BD_0 <- birthdeath(tree)
  out$b0 <- as.numeric(BD_0$para["b-d"] / (1 - BD_0$para["d/b"]))
  out$d0 <- as.numeric(BD_0$para["d/b"] * out$b0)
  out[["b0-d0"]] <- as.numeric(BD_0$para["b-d"])
  out[["d0/b0"]] <- as.numeric(BD_0$para["d/b"])
  out[["CI_b0-d0"]] <- BD_0$CI["b-d",]
  out[["CI_d0/b0"]] <- BD_0$CI["d/b",]

  out$gamma <- gammaStat(tree)
  # In Geiger (moment estimators)
  cat("Geiger... ", file=stderrf)
  tree$root.edge <- 0         # otherwise will compute the stem estimate
  out$net_d_ms <- bd.ms(tree)
  out$net_d_ms_missing <- bd.ms(tree, missing=round(Ntip(tree) * (1 - sampling_f)))
  out$net_d_km <- bd.km(tree) # missing not implemented for Kendall-Moran
  # In RPANDA
  cat("RPANDA... ", file=stderrf)
  cst.rate <- function(t,r){r[1]}
  tot_time <- max(leaf_heights(tree))
  BD_1 <- fit_bd(tree, tot_time, cst.rate, cst.rate, out$b0, out$d0,
                 cst.lamb=TRUE, cst.mu=TRUE)
  # With incomplete sampling
  BD_2 <- fit_bd(tree, tot_time, cst.rate, cst.rate, out$b0, out$d0, sampling_f,
                 cst.lamb=TRUE, cst.mu=TRUE)
  out$b1 <- BD_1$lamb_par
  out$d1 <- BD_1$mu_par
  out$b2 <- BD_2$lamb_par
  out$d2 <- BD_2$mu_par
  # In diversitree (MCMC)
  cat("Diversitree..", file=stderrf)
  bd_func <- make.bd(tree)
  bd_func_2 <- make.bd(tree, sampling.f=sampling_f)
  # Max-likelihood
  # TODO: catch warnings (improper likelihood convergence)
  cat("ML...", file=stderrf)
  BD_3 <- find.mle(bd_func, c(out$b0, out$d0), method="optim") #condition.surv=F
  BD_4 <- find.mle(bd_func_2, c(out$b0, out$d0), method="optim")
  # MCMC
  cat("MCMC...\n", file=stderrf)
  samples_5 <- mcmc(bd_func, as.numeric(c(out$b0, out$d0)), nsteps=10000, print.every=0,
                    w=c(0.1, 0.1))
  samples_6 <- mcmc(bd_func_2, as.numeric(c(out$b0, out$d0)), nsteps=10000, print.every=0,
                    w=c(0.1, 0.1))

  out <- c(out, coef(BD_3), coef(BD_4),
           summary_posterior(samples_5, 'lambda'),
           summary_posterior(samples_5, 'mu'),
           summary_posterior(samples_6, 'lambda'),
           summary_posterior(samples_6, 'mu'))

  return(unlist(out))
}


if( !interactive() ) {
  subtreelines  <- scan(subtreefile, what=character(), sep="\n")
  not_single_node_lines <- grepl('(', subtreelines, fixed=TRUE)
  subtreelines <- subtreelines[not_single_node_lines]

  subtrees <- read.tree(text=subtreelines, keep.multi=TRUE)
  names(subtrees) <- sapply(subtrees, function(tree){tree$node.label[1]})

  div_stats <- sapply(subtrees, get_div_stats, clade.data)
}



