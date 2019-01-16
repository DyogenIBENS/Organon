#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Annotate a time tree with the taxonomic names at internal nodes.
Assumes the tree topology matches the NCBI taxonomy topology.

Output to stdout."""


import argparse
from ete3 import PhyloTree, NCBITaxa
import logging
logger = logging.getLogger(__name__)
logging.basicConfig(format='%(levelname)s:%(message)s')
logger.setLevel(logging.INFO)


def get_last_common_element(*iterables):
    """Used to find the last common ancestor, given two lineages."""
    last = None
    for grouped in zip(*iterables):
        if len(set(grouped)) > 1:
            return last
        else:
            last = grouped[0]

    return last


def matchrename_ncbitax(timetree, uniq=True):
    """Apply scientific name as the node name.
    If a child has the same name as its parent and `uniq` is True, do not apply
    the scientific name."""

    if uniq:
        logger.info('Discarding duplicate common ancestors')

    # Apply NCBI name only when it corresponds to a NCBI taxonomy node.
    for node in timetree.traverse():
        node.add_feature('oldname', node.name)
        if not node.is_leaf():
            if not node.is_root():
                if node.up.taxid == node.taxid:
                    logger.warning("Node %r is a duplicated taxon %d (%r).",
                                   node.name, node.taxid, node.sci_name)
                    # Indicates a multifurcation in the NCBI taxonomy.
                    if uniq:
                        node.sci_name = ''

            if node.sci_name:
                node.name = node.sci_name


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


def name_ancestors(timetreefile, to_table=False, ete3_algo=False, uniq=True):
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
    if ete3_algo:
        ncbi.annotate_tree(timetree, 'taxid')
    else:
        myannotate(timetree, ncbi)
    matchrename_ncbitax(timetree, uniq)
    
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
    parser.add_argument('--ete3-algo', '--ete3', action='store_true',
                        help='Use the annotate function from ete3 instead of mine.')
    parser.add_argument('-d', '--duplicate', action='store_false', dest='uniq',
                        help='Keep duplicated consecutive annotations')
    
    args = parser.parse_args()
    name_ancestors(**vars(args))

