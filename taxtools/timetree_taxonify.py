#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Annotate a time tree with the taxonomic names at internal nodes.
Assumes the tree topology matches the NCBI taxonomy topology.

Output to stdout."""


import argparse
from ete3 import PhyloTree, NCBITaxa
import logging
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
#logging.basicConfig(format='%(levelname)s:l.%(lineno)d:%(message)s', level=logging.INFO)


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
    logger.info('Discarding duplicate common ancestors')
    # Apply NCBI name only when it corresponds to a NCBI taxonomy node.
    for node in timetree.traverse():
        node.add_feature('oldname', node.name)
        try:
            if not node.is_leaf() and \
                    (node.is_root() or node.up.taxid != node.taxid):
                node.name = node.sci_name
            else:
                node.sci_name = ''
        except BaseException as err:
            logger.error('%s %s %s leaf: %s', err, node.name, node.features, node.is_leaf())
            if node.up:
                logger.error('%s %s %s', err, node.up.name, node.up.features)
            raise
    #return timetree


def myannotate(timetree, ncbi):
    name2taxid = ncbi.get_name_translator([sp.replace('_', ' ') for sp in \
                                                    timetree.get_leaf_names()])
    
    seen = set()
    for node in timetree.traverse('postorder'):
        if node.is_leaf():
            taxids = name2taxid[node.name.replace('_', ' ')]
            assert len(taxids) == 1
            node.add_feature('taxid', taxids[0])
        else:
            try:
                lineages = [ncbi.get_lineage(ch.taxid) for ch in node.children]
            except BaseException as err:
                logger.error('%s %r %s\n%s in seen: %d', err, node.name,
                              [ch.name for ch in node.children],
                              [(ch.name in seen) for ch in node.children],
                              len(seen))
                raise

            lctaxid = get_last_common_element(*lineages)
            try:
                lca = ncbi.get_taxid_translator([lctaxid])[lctaxid]
            except BaseException as err:
                logger.error('%s: %s\n%s', err,
                              ', '.join(ch.name for ch in node.children),
                              lineages)
                raise

            node.add_feature('sci_name', lca)
            node.add_feature('taxid', lctaxid)

        seen.add(node.name)

    #return timetree


def name_ancestors(timetreefile, to_table=False):
    logger.info('Loading data')
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
            logger.warning('Species %r not found', leaf.name)
            leaf.delete(prevent_nondicotomic=True,
                        preserve_branch_length=True)

    logger.info('Placing common ancestors')
    #ncbi.annotate_tree(timetree, 'taxid')
    myannotate(timetree, ncbi)
    matchrename_ncbitax(timetree)
    
    #logger.debug({ft:getattr(timetree, ft) for ft in timetree.features})

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

