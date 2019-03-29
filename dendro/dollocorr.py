#!/usr/bin/env python3
# -*- coding: utf-8 -*-


"""Dollo's Law of Irreversibility: an organism never returns exactly to a
former state [...].

Implementation of the Maddison method to detect correlation between 2 binary traits,
where one is supposed to cause the other one. In this method, the explained trait can 
not be lost."""


import itertools as it
import numpy as np
from dendro.bates import rev_dfw_descendants
import logging
logger = logging.getLogger(__name__)
logging.basicConfig()


def integer_n_partition(total, nterms):
    """Iterator over all possible ways of summing to `total` in nterms,
    ordered, without zeros."""
    # Select nterms-1 split locations among [1, 2, ..., total-1]
    for isplits in it.combinations(range(1, total), nterms-1):
        splits = np.array( (0,) + isplits + (total,) )
        terms = splits[1:] - splits[:-1]
        yield terms
    ###FIXME: pb: I get all possible orderings, e.g (1,2) AND (2,1).
    ###       I only want each set once. (Because in the calling function, I
    ###       have already made all permutations).

def integer_partition(total, max_nterms=None):
    if max_nterms is None:
        max_nterms = total
    for nterms in range(1, max_nterms+1):
        yield from integer_n_partition(total, nterms)


def ordered_integer_n_partition(total, nterms, minval=1):
    """Iterate over all possibilities.
    Each terms of the summation is ordered."""
    assert nterms>0

    if nterms == 1:
        yield (total,)
    else:
        for i in range(minval, total//nterms + 1):
            for next_terms in ordered_integer_n_partition(total-i, nterms-1, i):
                yield (i,) + next_terms


def combinations_of_children_summing_to_k(chcounts, k):
    """Given children counts, compute the number of possibilities to observe
    k events in the tree rooted at the current node.
    
    param: chcounts := precomputed number of possibilities of any number of
                       events in the children pending subtrees. (children
                       of the current node)
                          
    Count all possible ways of observing k events at the current node.
    To do so, we need to decompose between:
    1. Observing a sum of k events in the pending subtrees, and no event in the
       branches leading to the children;
    2. observing i (1<=i<=k) events in the children branches, then k-i events
       in the pending subtrees.

    For each selected `i`, select all possible arrangements of
    (child branch, child nb of events) such that they sum to k-i.
    For a given particular arrangement, the number of possibilities is the
    product for each child of the number of possibilities of getting
    `child nb of events`.
    """
    # Generalization to any number of children.
    pcountk = 0

    # We have nch branches (one per children)
    # * We select i branches with one event each, then we get the k - i remaining events
    #   from the b-i remaining branches.
    nch = len(chcounts)
    branches = set(range(nch))
    for i in range( min(nch,k)+1 ):  # i of the k events happened on this fork.
        # For each selection of i branches among b:
        for event_branches in it.combinations(branches, i):
            noevent_branches = branches - set(event_branches)
            for permuted_noevent_br in it.permutations(noevent_branches):

                # How many ways of sharing k-i events between b-i branches?
                for nterms in range(1, min(k-i, b-i)):
                    for terms in ordered_integer_n_partition(k-i, nterms):
                        # Product of the number of ways of the descendants:
                        ki_noevent_count = 1
                        for neb, n_events in zip(permuted_noevent_br, terms):
                            ki_noevent_count *= chcounts[neb][n_events]
                        #for s0 in range(s, len(noevent_branches)):
                        #    ki_noevent_count *= 1
                        pcountk += ki_noevents_count
    return pcountk


def combinations_of_2_children_summing_to_k(chcounts, k):
    pass


get_phylchildren = lambda tree,node: [x for x,_ in tree.items.get(node,[])]


#@myTools.memoize
def place_single_events(n, phyltree, get_phylchildren=get_phylchildren):
    """Dynamic programming approach to place n events in a Maddison manner on a
    phylogenetic tree.

    Briefly, those events are considered unique/irreversible, so if it happened
    *once* in the ancestry, it can't happen again in any descendant.

    The tree must be a strictly dichotomic tree."""

    # Store intermediate results for each node.
    node_counts = {}

    #TODO: if n is None:
    # There can't be more single events than there are leaves.
    maxn = min(n, len(phyltree.listSpecies))

    # There is exactly 1 way to place 0 event on a single branch,
    # and we initialize 0 ways for all other numbers of events.
    init_count = np.array([1] + [0]*n, dtype=int)

    # Iterate from the leaves to the root.
    for parent, children in rev_dfw_descendants(phyltree, get_phylchildren,
                                                include_leaves=True):
        if not children:
            node_counts[parent] = init_count.copy()
            logger.debug('* Leaf %r: counts %s', parent, node_counts[parent])
            continue

        logger.debug('* %r -> %s', parent, children)
        assert len(children) == 2, "Not implemented for non dichotomic trees."\
                " (node %s)" % parent
        
        #nch = len(children)
        #chcounts = np.array([node_counts.pop(ch) for ch in children])
        try:
            c1, c2 = [node_counts.pop(ch) for ch in children]
        except ValueError as err:
            err.args = ("Not implemented for non dichotomic trees."\
                        " (node %s)" % parent,)
            raise


        logger.debug('children counts: %s: %s; %s: %s', children[0], c1,
                     children[1], c2)
                
        pcount = init_count.copy()

        # Iterating in decreasing order, because I need to call pcount[k+1], k+2
        # Exclude 0. The number of ways of placing 0 events is defined as 1.
        for k in range(maxn, 0, -1):
            logger.debug('  k=%d', k)
            # Zero event on the two children branches
            # All possible ways of sharing k events between 2 children

            # The following matrix product is equivalent to:
            #pcount[k] = 0
            #for i in range(0, k+1):
            #    pcount[k] += c1[i] * c2[k-i]

            pcount[k] = c1[:(k+1)].dot(c2[:(k+1)][::-1])
            #pcount[k] = combinations_of_2_children_summing_to_k(c1, c2, k)

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

