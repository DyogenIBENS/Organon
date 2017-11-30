#!/usr/bin/env Rscript

library(ape)

args <- commandArgs(trailingOnly=TRUE)

subtreefile <- args[1]
subtreefile <- "/users/ldog/glouvel/ws2/databases/timetree/Opisthokonta-listens90size10.subtrees.nwk"

if(len(args)>1) {
  tablefile <- args[2]
  tablefile <- "/users/ldog/glouvel/ws2/databases/timetree/Opisthokonta-listens90size10.tsv"
} else {
  tablefile <- NULL
}

fix_single_nodes <- function(lines){
  # newick line with a single node and no parenthesis can't be read by read.tree
  
}

subtreelines  <- scan(subtreefile, what=character(), sep="\n")
not_single_node_lines <- grepl('(', subtreelines, fixed=TRUE)
subtreelines <- subtreelines[not_single_node_lines]

subtrees <- read.tree(text=subtreelines, keep.multi=TRUE)
names(subtrees) <- sapply(subtrees, function(tree){tree$node.label[1]})

get_div_stats <- function(tree) {
  bd_out <- birthdeath(tree)
  return(c(bd_out$para, bd_out$se, bd_out$CI, gammaStat(tree)))
}

div_stats <-



