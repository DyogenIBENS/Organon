#!/usr/bin/env python3
# -*- coding: utf-8 -*-


"""Dollo's Law of Irreversibility: an organism never returns exactly to a
former state [...].

Implementation of the Maddison method to detect correlation between 2 binary traits,
where one is supposed to cause the other one. In this method, the explained trait can 
not be lost."""


import itertools as it
from copy import deepcopy
import numpy as np
from dendro.bates import rev_dfw_descendants
from LibsDyogen import myPhylTree

import logging
logger = logging.getLogger(__name__)
logging.basicConfig(format='%(levelname)s:%(module)s:%(funcName)s:%(message)s')


def nonnull_integer_n_partition(total: int, nterms: int):
    """Iterator over all possible ways of summing to `total` in nterms,
    ordered, **without zeros**.
    
    N.B.: yield all possible orderings, e.g (1,2) AND (2,1).
    """
    if total <= 0 or nterms > total or nterms < 0:
        raise ValueError("Partitioning %d in %d terms is impossible." % (total, nterms))

    # Select nterms-1 split locations among [1, 2, ..., total-1]
    for isplits in it.combinations(range(1, total), nterms-1):
        splits = np.array( (0,) + isplits + (total,) )
        terms = splits[1:] - splits[:-1]  # Convert split to block size.
        yield terms


def nonnull_integer_partition(total:int, max_nterms=None):
    if max_nterms is None:
        max_nterms = total
    for nterms in range(1, max_nterms+1):
        yield from nonnull_integer_n_partition(total, nterms)


def ordered_nonnull_integer_n_partition(total:int, nterms:int , minval=1):
    """Iterate over all possibilities.
    Each terms of the summation is ordered."""
    if nterms <= 0:
    #    yield ()
        raise StopIteration

    if nterms == 1:
        yield (total,)
    else:
        for i in range(minval, total//nterms + 1):
            for next_terms in ordered_nonnull_integer_n_partition(total-i, nterms-1, i):
                yield (i,) + next_terms


def integer_n_partition(total: int, nterms: int):
    """Iterator over all possible ways of summing to `total` in nterms,
    ordered, **without zeros**."""
    if nterms <= 0:
        if total == 0:
            yield np.array([], dtype=int)
            raise StopIteration
        else:
            raise StopIteration
            #raise ValueError("Partitioning %d in 0 terms is impossible." % total)

    # Select nterms-1 split locations among [1, 2, ..., total-1]
    for isplits in it.combinations_with_replacement(range(total+1), nterms-1):
        splits = np.array( (0,) + isplits + (total,) )
        terms = splits[1:] - splits[:-1]  # Convert split to block size.
        yield terms


def combinations_of_children_summing_to_k(chcounts, k:int):
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
    # Generalization to any number of children: "possibility count for k events"
    pcountk = 0

    # 1. We have nch branches (one per children).
    # 2. We select i branches with one event each,
    # 3. then we get the k - i remaining events from the b-i remaining branches.
    nch = len(chcounts)
    branches = set(range(nch))

    # The chcounts arrays must be extended to reach length k+1 (max index k)
    if min(len(c) for c in chcounts) < k+1:
        logger.warning("k is higher than the provided children k's, filling with 0.")
        chcounts = [np.append(c, np.zeros(k+1-len(c), dtype=int)) for c in chcounts]
    
    chcounts = np.array(chcounts)

    logger.debug('%d events among %d branches (maximum %d here).', k, nch, min(nch, k))

    # chcounts must have one more way to have 1 event (along the leading branch)
    for terms in integer_n_partition(k, nch):
        # Product of the number of ways of getting this partition
        # among the selected descendants:

        prodterms = chcounts[range(nch), terms]

        logger.debug('      + branches %s: %s', terms, '×'.join(str(x) for x in prodterms))
        pcountk += prodterms.prod()

    return pcountk


def combinations_of_2_children_summing_to_k(chcounts, k:int):
    try:
        c1, c2 = chcounts
    except ValueError:
        raise ValueError("This function expects exactly 2 children.")

    if k < 0:
        raise ValueError("Expecting k >= 0")
    #if k == 0:
    #    return 1  # Not if chcounts[:][0] == 0

    # The c1, c2 arrays must be extended to reach length k+1 (max index k)
    if min(len(c1), len(c2)) < k+1:
        logger.warning("k is higher than the provided children k's, filling with 0.")
        c1 = np.append(c1, np.zeros(k+1-len(c1), dtype=int))
        c2 = np.append(c2, np.zeros(k+1-len(c2), dtype=int))

    # The following matrix product is equivalent to:
    #pcount[k] = 0
    #for i in range(0, k+1):
    #    pcount[k] += c1[i] * c2[k-i]

    pcountk = np.dot(c1[:(k+1)], c2[k::-1])

    return pcountk


def update_all_combinations_dichotomic(chcounts, maxn, init_count):
    try:
        c1, c2 = chcounts
    except ValueError:
        raise ValueError("This function expects exactly 2 children.")

    pcount = init_count.copy()  # Faster than calling np.array() ?
    
    ## Iterating in decreasing order, because I need to call pcount[k+1], k+2
    ## Exclude 0. The number of ways of placing 0 events is defined as 1.
    #for k in range(maxn, 0, -1):
    #    logger.debug('  k=%d', k)
    #    # Zero event on the two children branches
    #    # All possible ways of sharing k events between 2 children

    #    # The following matrix product is equivalent to:
    #    #pcount[k] = 0
    #    #for i in range(0, k+1):
    #    #    pcount[k] += c1[i] * c2[k-i]

    #    pcount[k] = c1[:(k+1)].dot(c2[:(k+1)][::-1])
    #    #pcount[k] = combinations_of_2_children_summing_to_k(c1, c2, k)

    #    logger.debug('  %d ways of having: %d event(s) in the subtrees', pcount[k], k)
    #    #if k <= maxn-1:
    #    #    # One extra event occured on one of the children branches
    #    #    pcount[k+1] += c1[k] + c2[k]
    #pcount[1:(maxn+1)] += c1[:maxn] + c2[:maxn]
    #logger.debug('  +%s ways of having: %s event(s) in the subtrees +1 event',
    #             ','.join((c1[:maxn] + c2[:maxn]).astype(str)),
    #             ','.join(str(x) for x in range(maxn)))

    #if 2 <= maxn:
    #    # One way to place 2 events on this fork (and no more events leafward)
    #    logger.debug('  +1 way of having: 2 event(s) in this fork.')
    #    pcount[2] += 1
    
    for k in range(1, maxn+1):
        #combinations_of_2_children_summing_to_k(chcounts, k)
        pcount[k] = np.dot(c1[:(k+1)], c2[k::-1])

    return pcount


def update_all_combinations(chcounts, maxn, init_count):
    """Do not use"""

    pcount = init_count.copy()

    # Exclude 0. The number of ways of placing 0 events is defined as 1.
    for k in range(1, maxn+1):
        logger.debug('  k=%d', k)
        pcount[k] = combinations_of_children_summing_to_k(chcounts, k)

    return pcount


get_phylchildren = lambda tree,node: [x for x,_ in tree.items.get(node,[])]
#from dendro.any.phyltree import get_children as get_phylchildren

#@myTools.memoize
def place_single_events(n:int, phyltree, get_phylchildren=get_phylchildren,
                        root=None, leaf_counts=None):
    """Dynamic programming approach to place n events in a Maddison manner on a
    phylogenetic tree.

    Briefly, those events are considered unique/irreversible, so if it happened
    *once* in the ancestry, it can't happen again in any descendant.

    The tree must be a strictly dichotomic tree."""
    if root is None:
        root = phyltree.root

    # Store intermediate results for each node.
    node_counts = {}
    if leaf_counts is not None:
        node_counts.update(leaf_counts)

    #seen_nodes = {k: 1 for k in leaf_counts}

    #TODO: if n is None:
    # There can't be more single events than there are leaves.
    maxn = min(n, len(phyltree.species[root]))

    # There is exactly:
    # - 1 way to place 0 event on a single branch,
    # - 1 way to place 1 event on this single branch (on the branch itself)
    # and we initialize 0 ways for all other numbers of events (because there is no subtree).
    if maxn == 0:
        return np.ones(1, dtype=int)

    init_count = np.array([1, 1] + [0]*(n-1), dtype=int)

    #node_counts[root] = init_count.copy()  # In case the tree has 0 branches...
    # Iterate from the leaves to the root.
    for parent, children in rev_dfw_descendants(phyltree, get_phylchildren,
                                                include_leaves=False,
                                                queue=[root]):
        ## Check that this node name is unique, or make it so.
        ## Reminder: with myPhylTree, it is necessary unique.
        #try:
        #    # Assuming it is a string
        #    parent += '.%d' % seen_nodes[parent]
        #    seen_nodes[parent] += 1
        #except KeyError:
        #    seen_nodes[parent] = 1

        logger.debug('* %r -> %s', parent, children)
        #assert len(children) == 2, "Not implemented for non dichotomic trees."\
        #        " (node %s)" % parent

        #nch = len(children)
        chcounts = np.array([init_count]*len(children))
        for i, ch in enumerate(children):
            try:
                # The data is not a leaf, or was initialized via `leaf_counts`.
                chcounts[i] = node_counts.pop(ch)
            except KeyError:
                pass

        logger.debug('children counts: %s',
                     '; '.join('%s:%s' % t for t in zip(children, chcounts)))

        pcount = init_count.copy()

        try:
            for k in range(maxn+1):
                pcount[k] = combinations_of_2_children_summing_to_k(chcounts, k)
        except ValueError:
            #err.args = ("Not implemented for non dichotomic trees."\
            #            " (node %s)" % parent,)
            logger.info("Non dichotomic fork.")
            for k in range(maxn+1):
                pcount[k] = combinations_of_children_summing_to_k(chcounts, k)

        logger.debug("Counts at %r: %s", parent, pcount)
        # Proba that 1 event occurs on the *leading* branch:
        # Allowed only if zero events are allowed below!!
        if chcounts[:,0].all():
            pcount[1] += 1

        node_counts[parent] = pcount

    root_counts = node_counts.setdefault(root, init_count)
    root_counts[1] -= 1  # Consider the root has no leading branch.
    return root_counts


def detachAfter(phyltree, nodes):
    """Transform the given nodes into leaves (remove their subtrees)."""
    items = deepcopy(phyltree.items)
    officialname = deepcopy(phyltree.officialName)

    for node in nodes:
        for descendant in phyltree.allDescendants[node]:
            if descendant != node:
                officialname.pop(descendant)
            items.pop(descendant, None)  # Needed here because of an assertion test in
                                   # `reinitTree` that doesn't like disconnected nodes.

    newtree = myPhylTree.PhylogeneticTree((items, phyltree.root, officialname))
    newtree.reinitTree()
    return newtree


def update_subtree_combinations(n, tot, roots_1_counts, tree_0):
    p0 = 0
    p1_n = 0  # possibilities of a sum of n events in the subtrees at roots_1
    # All possible unordered ways to sum to `n` using n1 terms.
    for terms in integer_n_partition(n, n1):
        # Select the n1 terms. ('term' column for each row).
        p1_n += np.prod(roots_1_counts[range(n1),terms])

    if p1_n > 0:
        p0 = place_single_events(tot, tree_0, leaf_counts=detached_leaf_counts)

    return p0 * p1_n


def maddison_test(observed_1, tot, tree, roots_1, root=None, roots_0=None,
                  alternative='>=', check=False):
    """
    :param: `tree`:       The complete tree. Needed to compute possibilities
                          regardless of the branch background states.
    :param: `roots_1`:    list of node names whose leading branch+subtree have
                          the background state of 1
    :param: `observed_0`: Observed state changes in branches with background state 0
    :param: `observed_1`: Observed state changes in branches with background state 1.

    Return the probability of at least observed_1 changes in branches_1.
    """

    if root is None:
        root = tree.root
    if roots_0 is None:
        roots_0 = []

    observed_0 = tot - observed_1
    n1 = len(roots_1)
    n0 = len(roots_0)
    assert set(roots_0) <= set().union(*(tree.allDescendants[r] for r in roots_1))
    
    # Tree whose leaves end after the subtrees_1 (before the given roots_0)
    tree_1 = detachAfter(tree, roots_0) if roots_0 else tree
    tree_0 = detachAfter(tree_1, roots_1)

    # 1. Count all possibilities to get more than `observed_1` events below roots_1.

    # 1 row: counts of possibilities
    if roots_0:
        roots_0_counts = np.array([place_single_events(tot, tree, root=r)
                                   for r in roots_0])  #ndmin=2
        roots_0_counts[:,1] += 1  # One possibility of 1 event on the leading branches.
        init_leaf0_counts = np.zeros((n0, tot+1), dtype=int)
    else:
        roots_1_counts = np.array([place_single_events(tot, tree, root=r)
                                   for r in roots_1])
        roots_1_counts[:,1] += 1  # One possibility of 1 event on the leading branches.
    init_leaf1_counts = np.zeros((n1, tot+1), dtype=int)

    #if tree.root not in roots_1:
    #    p_leading_branch_events[1] += 1

    p1_more_than_obs = 0

    alternatives = {'>=': range(observed_1, tot+1),
                    '>': range(observed_1+1, tot+1),
                    '<=': range(observed_1+1),
                    '<': range(observed_1),
                    '=': (observed_1,)}
    # Here the number of computations would be optimized by choosing the shortest range,
    # and if needed subtracting it to the total.
    complementary = {'>=': '<', '>': '<=', '<=': '>', '<': '>='}

    for n in alternatives[alternative]:
        logger.debug('Proba of %d events in subtrees_1 & %d total events', n, tot)
        #p1_more_than_obs += update_subtree_combinations(n, tot, roots_1_counts, tree_0)
        
        # If there are reversions to background state 0 (roots_0):
        for terms0 in integer_n_partition(tot-n, n0+1):
            # Recompute the roots_1_counts based on those starting values
            logger.debug('# terms0: %s', terms0)
            if roots_0:
                selected_roots_0_counts = roots_0_counts[range(n0),terms0[:-1]]
                detached_leaf0_counts = init_leaf0_counts.copy()
                # Select the n1 terms. ('term' column for each row).
                detached_leaf0_counts[range(n0), terms0[:-1]] = selected_roots_0_counts
                detached_leaf0_counts = dict(zip(roots_0, detached_leaf0_counts))
                roots_1_counts = np.array([place_single_events(tot, tree_1, root=r,
                                                   leaf_counts=detached_leaf0_counts)
                                           for r in roots_1])
                roots_1_counts[:,1] += 1  # One possibility of 1 event on the leading branches.

            # All possible unordered ways to sum to `n` using n1 terms.
            for terms in integer_n_partition(n+terms0[:-1].sum(), n1):
                selected_roots_1_counts = roots_1_counts[range(n1),terms]
                logger.debug('⋅ selected roots_1 nb of events %s with possibilities: %s',
                             terms, '×'.join(str(x) for x in selected_roots_1_counts))
                if selected_roots_1_counts.all():
                    # Limit possibilities to the selected roots_1
                    detached_leaf1_counts = init_leaf1_counts.copy()
                    # Select the n1 terms. ('term' column for each row).
                    detached_leaf1_counts[range(n1), terms] = selected_roots_1_counts
                    detached_leaf1_counts = dict(zip(roots_1, detached_leaf1_counts))

                    # Possibilities of `tot` events constrained on `n` events after roots_1
                    ptot_n = place_single_events(tot, tree_0, root=root,
                                                 leaf_counts=detached_leaf1_counts)
                    logger.debug('+ %d', ptot_n[tot])

                    p1_more_than_obs += ptot_n[tot]

    ptot = place_single_events(tot, tree, root=root)[tot]

    return p1_more_than_obs, ptot


def maplosses(phyltree, leaf_states, root=None):
    """Simple parsimony algorithm where the (binary) character can only be
    lost (not regained)"""

    node_states = {}
    node_states.update(leaf_states)

    if root is None:
        root = phyltree.root

    loss_branches = set()

    for parent, children in rev_dfw_descendants(phyltree, get_phylchildren,
                                                queue=[root]):
        # The parent state is 0 if and only if both children are 0.
        node_states[parent] = any(node_states[ch] for ch in children)

        for ch in children:
            if node_states[parent] ^ node_states[ch]:
                loss_branches.add(ch)

    return loss_branches, node_states

