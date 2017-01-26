#!/usr/bin/env python2.7

import re
import sys
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
    # TODO: conserve branch length
    # Iterate from leaves to root
    for node in tree.iter_descendants("postorder"):
        if node.is_leaf():
            ancestor = convert_gene2species(node.name)
            genename = node.name
        else:
            ancestor, _ = split_species_gene(node.name)
        
        parent_node = node.up
        parent_ancestor, genename = split_species_gene(parent_node.name)

        try:
            ancestor_lineage = phyltree.dicLinks[parent_ancestor][ancestor] 
        except KeyError:
            print >>sys.stderr, "node       : %s (%r)" % (node.name, ancestor)
            print >>sys.stderr, "parent_node: %s (%r)" % (parent_node.name,
                                                          parent_ancestor)
            raise


        # If doesn't match species tree
        # same ancestor as child is possible for duplication node
        # So check for length 1 or 2
        if len(ancestor_lineage) > 2:
            # conserve original distance
            n_new_branches = len(ancestor_lineage) - 1
            dist_new_branches = float(node.dist) / n_new_branches
            # Add missing links
            for link in ancestor_lineage[1:-1]:
                parent_node = parent_node.add_child(name=(link + genename),
                                                    dist=dist_new_branches)
            # Move the node on top of the created intermediate links
            parent_node.add_child(child=node.detach(),
                                  dist=dist_new_branches)


def search_by_ancestorlist(tree, ancestorlist):
    def stop_at_any_ancestor(node):
        return any(node.name.startswith(anc) for anc in ancestorlist)
    return tree.iter_leaves(is_leaf_fn=stop_at_any_ancestor)


def save_subtrees_byspecieslist(tree, specieslist, outdir='.'):
    ancestor = PHYLTREE.lastCommonAncestor(specieslist)
    for node in search_by_ancestorspecies(tree, ancestor):
        outfile = os.path.join(outdir, node.name + '.nwk')
        node.write(format=1, outfile=outfile)

#def stop_at_duplicates(values):
#    for values in
def with_dup(leafnames):
    leafspecies = [convert_gene2species(leaf) for leaf in leafnames]
    return (len(leafspecies) > len(set(leafspecies)))


def save_subtrees(treefile, ancestors, outdir='.', only_dup=False,
                  dry_run=False):
    try:
        tree = ete3.Tree(treefile, format=1)
    except ete3.parser.newick.NewickError as e:
        print >>sys.stderr, 'ERROR with treefile %r' % treefile
        raise
    add_species_nodes_back(tree, PHYLTREE)
    for ancestor in ancestors:
        ancestor = ancestor.capitalize()
        ancestorlist = sorted(PHYLTREE.allDescendants[ancestor],
                              key=lambda anc: -PHYLTREE.ages[anc])
        ancestor_regex = re.compile('^(%s)(?=ENSGT)' % \
                                    '|'.join(ancestorlist).replace(' ', '.'))
        for node in search_by_ancestorlist(tree, ancestorlist):
            leafnames = node.get_leaf_names()
            #print node.name
            #print node.get_ascii()
            if len(leafnames) > 1:
                # check that there is at least one duplication
                if not only_dup or with_dup(leafnames):
                    outname = ancestor_regex.sub(node.name, ancestor)
                    outfile = os.path.join(outdir, node.name + '.nwk')
                    if dry_run:
                        print outfile
                    else:
                        node.write(format=1,
                                   format_root_node=True,
                                   outfile=outfile)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("treefile")
    parser.add_argument("ancestors", nargs='+')
    parser.add_argument("-o", "--outdir", default='.')
    parser.add_argument("--only-dup", action="store_true",
                        help="do not extract trees that don't have at least "\
                             "one duplication")
    parser.add_argument("-n", "--dry-run", action="store_true",
                        help="only print out the output files it would produce")
    #parser.add_argument("-e", "--ensembl-version", type=int, default=85)
    #parser.add_argument("-p", "--phyltree", help="")
    args = parser.parse_args()
    save_subtrees(**vars(args))

