#!/usr/bin/env python3
# -*- coding: utf-8 -*-


"""
Summarise tree topology and clade ages from MCMC output in nexus format.

Also see Beast2:TreeAnnotator
"""


from io import StringIO
from functools import partial
from collections import Counter, defaultdict
from Bio import Phylo
from Bio.Nexus import Nexus #, Nodes, Trees

from dendro.sorter import leaf_sort
from dendro.any import BioPhylo

from dendro.bates import dfw_descendants_generalized, rev_dfw_descendants


def node_topology(tree, node, get_children, get_label):
    children = get_children(tree, node)
    if children:
        return '(' + ','.join(node_topology(tree, ch, get_children, get_label)
                              for ch in children) + ')'
    else:
        return get_label(tree, node)


biophylo_topology = partial(node_topology,
                            get_children=BioPhylo.get_children,
                            get_label=BioPhylo.get_label)
biophylo_leaf_sort = partial(leaf_sort,
                             get_children=BioPhylo.get_children,
                             set_children=BioPhylo.set_children,
                             get_attribute=BioPhylo.get_label)

def represent_clade(tree, node, children):
    return tuple(sorted([tuple(sorted(iter_leaves(tree,
                                           BioPhylo.get_children,
                                           queue=[ch])))
                  for ch in children], key=min))


def represent_clades(tree, get_children, get_label):
    clades = {}
    for clade, children in rev_dfw_descendants(tree,
                                               get_children,
                                               include_leaves=True,
                                               queue=[tree.root]):
        cladename = get_label(tree, clade)

        if not children:
            clades[cladename] = ((cladename,),)
            continue

        # Make unique name
        rawcladename = cladename
        seen = 0
        namefmt = '%s.%d'
        while cladename in clades:
            seen += 1
            cladename = namefmt % (rawcladename, seen)
        clade.name = cladename

        # Concatenate each children subclades (flatten grand-children)
        clades[cladename] = []
        for child in children:
            childleaves = ()
            for subclade in clades[get_label(tree, child)]:
                childleaves += subclade
            clades[cladename].append(tuple(sorted(childleaves)))
        clades[cladename] = tuple(sorted(clades[cladename], key=min))

    return clades


def main(nexusfile, reftree, burnin=10):

    # Using the Nexus module
    data = Nexus.Nexus(nexusfile)
    taxlabels = data.structured[1].commandlines[1].options.split()
    nb2taxlabels = data.translate
    trees = data.trees
    # Using the Phylo module
    trees = list(Phylo.parse(nexusfile, 'nexus'))

    N0 = len(trees)

    trees = trees[N0*burnin/100+1:]
    N = N0*(100 - burnin)/100

    topologies = Counter()
    topo_groups = defaultdict(list)

    for tree in trees:
        # Ensure all equivalent topologies will be represented the same way
        biophylo_leaf_sort(tree, tree.root)
        topo = biophylo_topology(tree, tree.root)
        topologies[topo] += 1
        topo_groups.append(tree)

    MAP_topology, MAP_count = topologies.most_common(1)[0]
    MAP_proba = float(MAP_count) / sum(topologies.values())

    clades = represent_clades(reftree, BioPhylo.get_children, BioPhylo.get_label)

