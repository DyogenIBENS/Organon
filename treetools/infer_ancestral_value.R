#!/usr/bin/env Rscript

library(ape)
library(scales)

infer_anc <- function(treefile, valuefile) {
  values <- read.delim2(valuefile, header=F, row.names=1)
  tree <- read.tree(treefile)
  anc <- ace(gl[,1], tree)
  return(list(tree=tree, values=values, anc=anc))
}

colorizevalues <- function(values, colours=c('#006837', '#FFFFBF', '#A50026')) {
  cramp <- colour_ramp(colours)
  minval <- min(values)
  maxval <- max(values)

  if(maxval == minval) {
    return('black')
  } else {
    return(cramp((values - minval)/(maxval - minval)))
  }
}

#plot_anc <- function(anc, colours=c('DarkGreen', 'yellow', 'DarkRed')) {
plot_anc <- function(tree, nodevalues, colours=c('#006837', '#FFFFBF', '#A50026'), ...) {
  edge.color <- colorizevalues(nodevalues, colours)
  plot(tree, edge.color=edge.color, ...)
}

prune <- function(tree, nodes) {
  ntips = Ntip(tree)
  for (node in nodes) {
    node <- ifelse(is.character(node), which(tree$node.label == node), node)
    edges <- which.edge(tree, tips(tree, node))
    tree$edge.length[edges] <- -1
  }
  return(di2multi(tree, tol=-0.5))
}


hide_edges <- function(tree, by) {
  edge.width <- rep(1, nrow(tree$edge))
  for(group in unique(by)){
    group_tips <- which(group == by)
    if ( length(group_tips) > 1 )
      edge.width[which.edge(tree, group_tips)] <- 0 
  }
  return(edge.width)
}

renamebyclade <- function(tree, by) {
  # Rename tips by the clade name (only the first tip belonging has the name)
  # Useful for plotting
  dupfam <- duplicated(by)
  tree$tip.label <- as.character(by)
  #by[dupfam]
  tree$tip.label[dupfam] <- ""
  return(tree)
}


plot_triangle <- function(tree, by, ...) {
  # Adapted from https://stackoverflow.com/a/34405198/4614641
  #groups <- c("A", "B", "C", "D")

  edge.width <- hide_edges(tree, by)
  
  # Step 2: plot the tree with the hidden edges
  treecopy <- renamebyclade(tree)
  plot(treecopy, show.tip.label=TRUE, edge.width=edge.width, ...)
  last.phylo.plot <- get("last_plot.phylo", envir=.PlotPhyloEnv)

  # Step 3: add triangles
  add_polytomy_triangle <- function(phy, tip_group, last.phylo.plot){
    # tip_group: a list of tip numbers
    root <- Ntip(phy)+1
    #group_node_labels <- phy$tip.label[group]
    #group_nodes <- which(phy$tip.label %in% group_node_labels)
    group_mrca <- getMRCA(phy,tip_group) #group_nodes)

    x0 <- last.phylo.plot$xx[group_mrca]
    y0 <- last.phylo.plot$yy[group_mrca]
    x1 <- max(last.phylo.plot$xx[tip_group])
    y1 <- tip_group[1]
    y2 <- tip_group[length(tip_group)]

    xcoords <- c(x0, x1, x1)
    ycoords <- c(y0, y1, y2)
    polygon(xcoords, ycoords, col="gray", border="gray")
    return(list(x=xcoords, y=ycoords))
  }
  for(group in unique(by)){
    group_tips <- which(by == group)
    if (length(group_tips) > 1) {
      coords <- add_polytomy_triangle(tree, group_tips, last.phylo.plot)
      text(coords$x[2], mean(coords$y[2:3]), group, adj=0)
    }
  }
}

collapse_groups <- function(tree, by) {
  ### `clades` is a column whose rownames match exactly the tip.labels
  ### FIXME
  ntips = Ntip(tree)
  for (group in unique(by)) {
    edges <- which.edge(tree, which(group == by))
    tree$edge.length[edges] <- -1
  }
  polytomictree <- di2multi(tree, tol=-0.5)
  polytomictree$edge.length[polytomictree$edge.length < 0] <- 0
  return(polytomictree)
}


rodents_params <- list(treefile="~/ws2/databases/GL_rodentsALL.nwk",
                       valuefile="~/ws2/databases/GL_rodentsALL-values.tsv",
                       cladefile="~/ws2/databases/GL_rodentsALL-clades.tsv")

