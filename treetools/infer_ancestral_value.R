#!/usr/bin/env Rscript

library(ape)
library(geiger)
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

get_descendant_groups <- function(tree, ancestors) {
  # Label leaves with the name of the chosen ancestral node.
  # ancestors: a list of ancestral nodes.
  anc_tip_labels <- tree$tip.label
  for(anc in ancestors) {
    node_num <- which(anc == tree$node.label) + Ntip(tree)
    if(length(node_num)) {
      tip_pos <- which(tree$tip.label %in% tips(tree, node_num))
      anc_tip_labels[tip_pos] <- anc
    } else {
      warning("Node not found: ", anc, call.=FALSE)
    }
  }
  return(anc_tip_labels)
}


hide_edges <- function(tree, by) {
  # by: vector/factor of the same lengths as Ntip(tree), grouping leaves.
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
  dupfam <- rev(duplicated(rev(by))) # reverse because tips are plotted bottom to top
  tree$tip.label <- as.character(by)
  #by[dupfam]
  tree$tip.label[dupfam] <- ""
  return(tree)
}


plot_triangle <- function(tree, by, ...) {
  # 
  # by: tip labels (groupings) of the same length as Ntip(tree)
  # ... Arguments to be passed to `plot.phylo`
  
  # Adapted from https://stackoverflow.com/a/34405198/4614641
  #groups <- c("A", "B", "C", "D")

  edge.width <- hide_edges(tree, by)
  
  # Step 2: plot the tree with the hidden edges
  treecopy <- renamebyclade(tree, by)
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

polytomize_by <- function(tree, by) {
  ### `clades` is a column whose rownames match exactly the tip.labels
  ntips = Ntip(tree)
  for (group in unique(by)) {
    edges <- which.edge(tree, which(group == by))
    tree$edge.length[edges] <- -1
  }
  polytomictree <- di2multi(tree, tol=-0.5)
  polytomictree$edge.length[polytomictree$edge.length < 0] <- 0
  return(renamebyclade(polytomictree, by))
}

collapse_by <- function(tree, by) {
  dup <- duplicated(by)
  todelete <- which(dup)
  #for (group in unique(by)) {
  #  group_tips <- which(by == group)
  #}
  tree <- drop.tip(tree, todelete)
  #tree$tip.label <- as.character(by[!dup])
  return(renamebyclade(tree, by))
}

ancestors <- function(tree, node, to=(Ntip(tree)+1), names=TRUE) {
  # List the ancestors (names or numbers) between `node` and `to`.
  alllabels <- c(tree$tip.label, tree$node.label)
  print(c(node, to))
  print(which(node == alllabels))
  if (!is.numeric(node)) {node <- which(node == alllabels)[0]}
  print(which(node == alllabels))
  print(c(node, to))
  if (!is.numeric(to)) {to <- which(to == alllabels)[0]}
  print(c(node, to))

  parents <- integer(0)

  #print(to)
  while (node != to) {
    parents <- c(parents, node)
    node <- tree$edge[tree$edge[,2] == node, 1]
  }
  parents <- c(parents, node)
  if (names) {parents <- alllabels[parents]}
  return(parents)
}


format_stabletraitsout <- function(tree, brlensfile) {
  brlens <- read.delim(brlensfile, header=T)
  #branches <- strsplit(as.character(brlens$Branch), "->")
  # edges should be in the same order
}


rodents_params <- list(treefile="~/ws2/anctraits/GL_rodentsALL.nwk",
                       valuefile="~/ws2/anctraits/GL_rodentsALL-values.tsv",
                       cladefile="~/ws2/anctraits/GL_rodentsALL-clades.tsv",
                       brlensfile="~/ws2/anctraits/stabletraitsout/GL_rodents.brlens",
                       brlensfileBM="~/ws2/anctraits/stabletraitsout/GL_rodents_BM.brlens",
                       ancstatesfile="~/ws2/anctraits/stabletraitsout/GL_rodents.ancstates",
                       ancstatesfileBM="~/ws2/anctraits/stabletraitsout/GL_rodents_BM.ancstates",
                       treefileST="~/ws2/anctraits/stabletraitsout/GL_rodents.tree")

# Example use:
if (F) {
  w <- 12 # plot width
  h <- 36 #      height

  # Estimation using ape::ace (phylogenetic independent contrasts)
  ttree <- read.tree("GL_rodentsALL.taxa.nwk")
  ttreefam <- renamebyclade(ttree, clade.data$family)

  ancGL.pic <- ace(gl[,1], ttree, method="pic")
  
  edge.color.pic <- colorizevalues(ancGL.pic$ace)
  pdf("~/ws2/anctraits/GL_rodents_treeplot-pic.pdf", width=w, height=h)
  plot(ttreefam, edge.color=edge.color.pic, main="PIC estimate", cex=0.5)
  dev.off()

  # Estimation from StableTraits
  treest <- read.tree(rodents_params$treefileST)
  treestfam <- renamebyclade(treest, clade.data$family)

  st_rates <- read.delim(rodents_params$brlensfile, header=T)
  st_anc <- read.delim(rodents_params$ancstatesfile, header=T)
  st_ancBM <- read.delim(rodents_params$ancstatesfileBM, header=T)

  edge.color.anc <- colorizevalues(st_anc$Median[4:2752])
  edge.color.ancBM <- colorizevalues(st_ancBM$Median[4:2752])

  pdf("~/ws2/anctraits/GL_rodents_treeplot-STBM.pdf", width=w, height=h)
  plot(treest, edge.color=edge.color.ancBM, main="Brownian with StableTraits",
       cex=0.5)
  dev.off()
  pdf("~/ws2/anctraits/GL_rodents_treeplot-ST.pdf", width=w, height=h)
  plot(treest, edge.color=edge.color.anc, main="StableTraits", cex=0.5)
  dev.off()

  # or display with phytools::contMap
}
