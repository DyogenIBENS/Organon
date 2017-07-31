#!/usr/bin/env Rscript

library(ape)
library(scales)

infer_anc <- function(treefile, valuefile) {
  values <- read.delim2(valuefile, header=F, row.names=1)
  tree <- read.tree(treefile)
  anc <- ace(gl[,1], tree)
  return(list(tree=tree, values=values, anc=anc))
}

#plot_anc <- function(anc, colours=c('DarkGreen', 'yellow', 'DarkRed')) {
plot_anc <- function(tree, nodevalues, colours=c('#006837', '#FFFFBF', '#A50026')) {
  cramp <- colour_ramp(colours)
  minval <- min(nodevalues)
  maxval <- max(nodevalues)

  edge.color <- if(maxval == minval) {'black'} else {
                    cramp((nodevalues - minval)/(maxval - minval))}
  plot(tree, edge.color=edge.color)
}

rodents_params <- list(treefile="~/ws2/databases/GL_rodentsALL.nwk",
                       valuefile="~/ws2/databases/GL_rodentsALL-values.tsv")

