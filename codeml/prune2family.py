#!/usr/bin/env python2.7

#from Bio import Phylo
#
#def load_tree(filename):
#    tree = Phylo.read(filename, "newick")
#    return tree

import os.path
import argparse

import ete3
import LibsDyogen.utils.myPhylTree as PhylTree
from select_leaves_from_specieslist import convert_gene2species


ENSEMBL_VERSION = 85
PHYLTREE = PhylTree.PhylogeneticTree("/users/ldog/alouis/ws2/GENOMICUS_SVN/data{0}/PhylTree.Ensembl.{0}.conf".format(ENSEMBL_VERSION))


#def prune_tree(treefile, specieslist):
#    """Isolate gene subtrees contained in monophyletic species clades.
#    Generate several substrees per original file"""
#    root = ete3.Tree(treefile, format=1)

def split_species_gene(nodename):
    """When genename is the concatenation of the species and gene names"""
    idx = nodename.index('ENSGT')
    return nodename[:idx].replace('.', ' '), nodename[idx:]

def add_species_nodes_back(tree, phyltree):
    """Add missing species ancestors in gene tree"""
    # Iterate from leaves to root
    for node in tree.iter_descendants("postorder"):
        if node.is_leaf():
            ancestor = convert_gene2species(node.name)
            genename = node.name
        else:
            ancestor, _ = split_species_gene(node.name)
        
        parent_node = node.up
        parent_ancestor, genename = split_species_gene(parent_node.name)

        ancestor_lineage = phyltree.dicLinks[parent_ancestor][ancestor] 

        # If doesn't match species tree
        # same ancestor as child is possible for duplication node
        # So check for length 1 or 2
        if len(ancestor_lineage) > 2:
            # Add missing links
            for link in ancestor_lineage[1:-1]:
                parent_node = parent_node.add_child(name=(link + genename))
            # Move the node on top of the created intermediate links
            parent_node.add_child(child=node.detach())


def search_by_ancestorspecies(tree, ancestor):
    for node in tree.traverse():
        if node.name.startswith(ancestor):
            yield node

def save_subtrees_byspecieslist(tree, specieslist, outdir='.'):
    ancestor = PHYLTREE.lastCommonAncestor(specieslist)
    for node in search_by_ancestorspecies(tree, ancestor):
        outfile = os.path.join(outdir, node.name + '.nwk')
        node.write(format=1, outfile=outfile)


def save_subtrees(treefile, ancestors, outdir='.'):
    tree = ete3.Tree(treefile, format=1)
    add_species_nodes_back(tree, PHYLTREE)
    for ancestor in ancestors:
        ancestor = ancestor.capitalize()
        for node in search_by_ancestorspecies(tree, ancestor):
            if len(node.get_leaves()) > 1:
                outfile = os.path.join(outdir, node.name + '.nwk')
                node.write(format=1, outfile=outfile)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("treefile")
    parser.add_argument("ancestors", nargs='+')
    parser.add_argument("-o", "--outdir", default='.')
    args = parser.parse_args()
    save_subtrees(**vars(args))


