#!/usr/bin/env python3
# -*- coding: utf-8 -*-


"""Dollo's Law of Irreversibility: an organism never returns exactly to a
former state [...].

Implementation of the Maddison method to detect correlation between 2 binary traits,
where one is supposed to cause the other one. In this method, the explained trait can 
not be lost."""


import numpy as np
from dendron.climber import rev_dfw_descendants
import logging
logger = logging.getLogger(__name__)
logging.basicConfig()


get_phylchildren = lambda tree,node: [x for x,_ in tree.items.get(node,[])]


#@myTools.memoize
def place_single_events(n, phyltree, get_phylchildren=get_phylchildren):
    """Dynamic approache to place n events in a Maddison manner on a
    phylogenetic tree"""

    node_counts = {}
    #TODO: if n is None:
    # There can't be more single events than there are leaves.
    maxn = min(n, len(phyltree.listSpecies))

    init_count = np.array([1] + [0]*n)

    for parent, children in rev_dfw_descendants(phyltree, get_phylchildren,
                                                include_leaves=True):
        if not children:
            # There is exactly one way to place 0 event on a single branch,
            # and one way to place one single event on a single branch.
            node_counts[parent] = init_count.copy()
            logger.debug('* Leaf %r: counts %s', parent, node_counts[parent])
            continue

        logger.debug('* %r -> %s', parent, children)
        assert len(children) == 2, "Not implemented for non dichotomic trees."\
                " (node %s)" % parent
        
        c1 = node_counts.pop(children[0])
        c2 = node_counts.pop(children[1])
        logger.debug('children counts: %s: %s; %s: %s', children[0], c1,
                     children[1], c2)
                
        pcount = init_count.copy()

        # Iterating in decreasing order, because I need to call pcount[k+1], k+2
        # Exclude 0. The number of ways of placing 0 events is defined as 1.
        for k in range(maxn, -1, -1):
            logger.debug('  k=%d', k)
            # Zero event on the two children branches
            # All possible ways of sharing k events between 2 children
            if k>0:
                # The following matrix product is equivalent to:
                #pcount[k] = 0
                #for i in range(0, k+1):
                #    pcount[k] += c1[i] * c2[k-i]

                pcount[k] = c1[:(k+1)].dot(c2[:(k+1)][::-1])

            logger.debug('  %d ways of having: %d event(s) in the subtrees', pcount[k], k)
            #if k <= maxn-1:
            #    # One extra event occured on one of the children branches
            #    pcount[k+1] += c1[k] + c2[k]
        pcount[1:(maxn+1)] += c1[:maxn] + c2[:maxn]
        logger.debug('  +%s ways of having: %s event(s) in the subtrees +1 event',
                     ','.join((c1[:maxn] + c2[:maxn]).astype(str)),
                     ','.join(str(x) for x in range(maxn)))

        if 2 <= maxn:
            # One way to place 2 events on this fork (and no more events leafward)
            logger.debug('  +1 way of having: 2 event(s) in this fork.')
            pcount[2] += 1

        logger.debug("Counts at %r: %s", parent, pcount)

        node_counts[parent] = pcount

    assert parent == phyltree.root
    return node_counts[parent]

