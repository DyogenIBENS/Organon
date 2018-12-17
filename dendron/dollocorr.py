#!/usr/bin/env python3
# -*- coding: utf-8 -*-


"""Dollo's Law of Irreversibility: an organism never returns exactly to a
former state [...].

Implementation of the Maddison method to detect correlation between 2 binary traits,
where one is supposed to cause the other one. In this method, the explained trait can 
not be lost."""


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
    # n = len(leaves)

    for parent, children in rev_dfw_descendants(phyltree, get_phylchildren,
                                                include_leaves=True):
        if not children:
            # There is exactly one way to place 0 event on a single branch,
            # and one way to place one single event on a single branch.
            #node_counts[parent] = {0:1, 1:1}
            node_counts[parent] = {0:1}  # node without descendant branches.
            #node_counts[parent] = np.array([1] + [0]*(n-1))
            logger.debug('* Leaf %r: counts %s', parent, node_counts[parent])
            continue

        logger.debug('* %r -> %s', parent, children)
        assert len(children) == 2, "Not implemented for non dichotomic trees."\
                " (node %s)" % parent
        
        c1 = node_counts[children[0]]
        c2 = node_counts[children[1]]
        logger.debug('children counts: %s: %s; %s: %s', children[0], c1,
                     children[1], c2)
                
        pcount = {0: 1}

        # Iterating in decreasing order, because I need to call pcount[k+1], k+2
        # Exclude 0. The number of ways of placing 0 events is defined as 1.
        for k in range(n, -1, -1):
            logger.debug('  k=%d', k)
            # Zero event on the two children branches
            # All possible ways of sharing k events between 2 children
            if k>0: #if k not in pcount:
                pcount[k] = 0
                for i in range(0, k+1):
                    logger.debug('    i=%d', i)
                    
                    #if c1.get(i) and c2.get(k-i):  # Unnecessary step if properly initialized.
                    #    pcount[k] += c1[i] * c1[k-i]
                    pcount[k] += c1.get(i,0) * c2.get(k-i,0)

            logger.debug('  %d ways of having: %d event(s) in the subtrees', pcount[k], k)
            if k <= n-1:
                # One extra event occured on one of the children branches
                pcount[k+1] += c1.get(k,0) + c2.get(k,0)
                logger.debug('  +%d ways of having: %d event(s) in the subtrees +1 event',
                             c1.get(k,0) + c2.get(k, 0), k)
        if 2 <= n:
            # One way to place 2 events on this fork (and no more events leafward)
            logger.debug('  +1 way of having: 2 event(s) in this fork.')
            pcount[2] += 1

        logger.debug("Counts at %r: %s", parent, pcount)

        node_counts[parent] = pcount

    assert parent == phyltree.root
    return node_counts[parent]






