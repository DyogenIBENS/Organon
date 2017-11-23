#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Tool for phylogenetic comparative testing: Clades-collapse-summarize.

1 - Compute statistics in a phylogeny for each clade of a given rank
    (e.g, each family).

    Statistics are:
        - clade age
        - number of extant species
        - diversification rate
    
    Return a table with values corresponding to each clade.

2 - Return the tree with collapsed clades.
"""

import sys
import argparse
import ete3
# Will use PhyloTree, NCBITaxa

from itertools import chain


def match_duplicate_taxid(taxids, node, taxid2name):
    """Find the correct taxid corresponding to node name, by comparing the node
    ancestors with the taxid lineage."""
    lineages = ncbi.get_lineage_translator(taxids)
    named_lin = {t: set(taxid2name.get(lt) for lt in lin) for t,lin \
                                        in lineages.items()}
    ancestors = [n.name for n in node.get_ancestors()]

    # Find the correct taxid: the one whose lineage (ncbi) matches
    # first to the most recent ancestor (in tree).

    match = None
    while not match:
        anc = ancestors.pop()
        for taxid, nl in named_lin.items():
            if anc in nl:
                match = taxid
                print('Choose %s' % taxid, file=sys.stderr)
                break

    return match


def main(inputtree, outbase, rank='family'):
    tree = ete3.PhyloTree(inputtree, format=1, quoted_node_names=False)

    print("Loading taxonomy", file=sys.stderr)
    ncbi = ete3.NCBITaxa()

    name2taxid = ncbi.get_name_translator([node.name for node in \
                                                    tree.traverse()])
    # Won't return anything for names not found

    taxid2rank = ncbi.get_rank(chain(*name2taxid.values()))
    taxid2name = ncbi.get_taxid_translator(chain(*name2taxid.values()))
    
    # Could also simply annotate the tree using ncbi.annotate_tree
    ncbi.annotate_tree(tree, 'taxid')

    def is_leaf_fn(node):
        return node.rank == rank

    #def is_leaf_fn(node):
    #    """Stop at given rank"""
    #    taxids = name2taxid.get(node.name)
    #    if taxids is not None:
    #        if len(taxids) > 1:
    #            print('WARNING: non unique name %r: %s.' % (node.name, taxids),
    #                  end=' ', file=sys.stderr)
    #            taxids = [match_duplicate_taxid(taxids, node, name2taxid)]

    #        noderank = taxid2rank[taxids[0]]
    #        return noderank == rank
    #    else:
    #        return False
    
    #subtrees = []

    with open(outbase + '-%s.tsv' % rank, 'w') as outtsv, \
         open(outbase + '-%s.subtrees.nwk' % rank, 'w') as outsub:

        outtsv.write('%s\tsize\tage\tdiv_rate\n' % rank)
        
        for node in tree.iter_leaves(is_leaf_fn):
            #subtrees.append(node.copy())
            outsub.write(node.write(format=1, format_root_node=True) + '\n')
            
            # Collapse
            #for child in node.children:
            #    node.remove_child(child)
            size = len(node)
            ### NOTE: divide by stem or crown age??
            _, age = node.get_farthest_leaf()
            div_rate = float(size) / age if age else ''
            
            outtsv.write('\t'.join(
                         [node.name, str(size), str(age), str(div_rate)])
                         + '\n')


    tree.write(outfile=(outbase + '-%s.nwk' % rank), format=1, is_leaf_fn=is_leaf_fn,
                format_root_node=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__,
                        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('inputtree')
    parser.add_argument('outbase',
                        help='basename for the 2 output files (.tsv, .nwk)')
    parser.add_argument('-r', '--rank', default='family',
                        choices=['phylum', 'class', 'family', 'order', 'genus'])
    
    args = parser.parse_args()
    main(**vars(args))
