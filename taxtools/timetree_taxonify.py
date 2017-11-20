#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Annotate a time tree with the taxonomic names at internal nodes.
Assumes the tree topology matches the NCBI taxonomy topology.

Output to stdout."""


import sys
import argparse
from ete3 import PhyloTree, NCBITaxa


def get_last_common_element(*iterables):
    """Used to find the last common ancestor, given two lineages."""
    last = None
    for grouped in zip(*iterables):
        if len(set(grouped)) > 1:
            return last
        else:
            last = grouped[0]

    return last


def matchrename_ncbitax(timetree):
    print('Discarding duplicate common ancestors', file=sys.stderr)
    # Apply NCBI name only when it corresponds to a NCBI taxonomy node.
    for node in timetree.traverse():
        node.add_feature('oldname', node.name)
        try:
            if not node.is_leaf() and \
                    (node.is_root() or node.up.taxid != node.taxid):
                node.name = node.sci_name
            else:
                node.sci_name = ''
        except:
            print(node.name, node.features, 'leaf:', node.is_leaf(), file=sys.stderr)
            if node.up:
                print(node.up.name, node.up.features, file=sys.stderr)
            raise
    #return timetree


def myannotate(timetree, ncbi):
    name2taxid = ncbi.get_name_translator([sp.replace('_', ' ') for sp in \
                                                    timetree.get_leaf_names()])
    
    seen = set()
    for node in timetree.traverse('postorder'):
        if node.is_leaf():
            node.add_feature('taxid', name2taxid[node.name.replace('_', ' ')][0])
        else:
            try:
                lineages = [ncbi.get_lineage(ch.taxid) for ch in node.children]
            except:
                print(repr(node.name), [ch.name for ch in node.children], file=sys.stderr)
                print([(ch.name in seen) for ch in node.children], 'in seen:',
                      len(seen), file=sys.stderr)
                raise

            lctaxid = get_last_common_element(*lineages)
            try:
                lca = ncbi.get_taxid_translator([lctaxid])[lctaxid]
            except:
                print(', '.join(ch.name for ch in node.children), file=sys.stderr)
                print(lineages, file=sys.stderr)
                raise

            node.add_feature('sci_name', lca)
            node.add_feature('taxid', lctaxid)

        seen.add(node.name)

    #return timetree


def name_ancestors(timetreefile, to_table=False):
    print('Loading data', file=sys.stderr)
    ### /!\ quoted_node_names only from ete3 v3.1.1
    timetree = PhyloTree(timetreefile, format=1,
                         quoted_node_names=True)
    ncbi = NCBITaxa()
    

    name2taxid = ncbi.get_name_translator([sp.replace('_', ' ') for sp in \
                                                    timetree.get_leaf_names()])
    
    for leaf in timetree.get_leaves():
        try:
            leaf.add_feature('taxid', name2taxid[leaf.name.replace('_', ' ')][0])
        except KeyError:
            print('WARNING: species %r not found' % leaf.name,
                  file=sys.stderr)
            leaf.delete(prevent_nondicotomic=True,
                        preserve_branch_length=True)


    # TODO: use ncbi.annotate_tree
    print('Placing common ancestors', file=sys.stderr)
    ncbi.annotate_tree(timetree, 'taxid')
    #myannotate(timetree, ncbi)
    matchrename_ncbitax(timetree)
    
    #print({ft:getattr(timetree, ft) for ft in timetree.features},
    #        file=sys.stderr)

    if not to_table:
        print(timetree.write(format=1, format_root_node=True))
    else:
        for node in timetree.traverse():
            if not node.is_leaf():
                print(node.oldname + '\t' + getattr(node, 'sci_name', ''))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('timetreefile')
    parser.add_argument('--to-table', action='store_true')
    
    args = parser.parse_args()
    name_ancestors(**vars(args))
    
