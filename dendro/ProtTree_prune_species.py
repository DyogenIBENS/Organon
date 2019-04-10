#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from sys import stdout, stdin, setrecursionlimit
import argparse as ap
from LibsDyogen import myPhylTree, myProteinTree
from dendro.trimmer import thin_prottree
import logging
logger = logging.getLogger(__name__)


def fix_thinned_dups(phyltree, tree, node):
    logger.debug('Entering fix_thinned_dup at %d', node)
    for child,_ in tree.data.get(node, []):
        fix_thinned_dups(phyltree, tree, child)

    nodeinfo = tree.info[node]
    if nodeinfo['Duplication'] != 0 and nodeinfo['taxon_name'] not in phyltree.allNames:
        child_taxa = [tree.info[ch]['taxon_name'] for ch,_ in tree.data.get(node,[])]
        #child_taxa = [t for t in child_taxa if t in phyltree.allNames]
        #assert len(set(child_taxa)) == 1, '%s:%s' % (nodeinfo['taxon_name'],
        #                                             set(child_taxa))
        #nodeinfo['taxon_name'] = child_taxa[0]
        logger.info('Reset taxon: %s -> %s (%s)',
                    nodeinfo['taxon_name'],
                    phyltree.lastCommonAncestor(child_taxa),
                    child_taxa)
        nodeinfo['taxon_name'] = phyltree.lastCommonAncestor(child_taxa)



def main(phyltreefile, forestfile=None):
    #with open(badspecieslistfile) as f:
    #    badspecies = [line.rstrip() for line in f if not line.startswith('#')]
    phyltree = myPhylTree.PhylogeneticTree(phyltreefile)

    if forestfile is None:
        forestfile = stdin
    for tree in myProteinTree.loadTree(forestfile):
        keptleaves = set((leaf for leaf in set(tree.info).difference(tree.data)
                          if tree.info[leaf]['taxon_name'] in phyltree.allNames))
        newroot, _ = thin_prottree(tree, tree.root, 0, keptleaves)
        #print('DEBUG: newroot =', newroot)
        #print('DEBUG: newdata =', tree.data)
        #print('DEBUG: newinfo =', ' '.join(str(x) for x in tree.info.keys()))
        if newroot is not None:
            fix_thinned_dups(phyltree, tree, newroot)
            tree.printTree(stdout, newroot)
        else:
            logger.warning('Discard tree %d', tree.root)


if __name__ == '__main__':
    setrecursionlimit(10000)
    logging.basicConfig()
    logger.setLevel(logging.INFO)
    parser = ap.ArgumentParser(description=__doc__)
    parser.add_argument('forestfile', nargs='?')
    parser.add_argument('phyltreefile')
    
    args = parser.parse_args()
    main(**vars(args))

