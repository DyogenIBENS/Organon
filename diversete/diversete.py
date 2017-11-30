#!/usr/bin/env python3
# -*- coding: utf-8 -*-


"""Utility functions for ete3 objects, to perform diversification analyses:
    
    - diversification net rate
    - constant rate BD fit
    - gamma statistic"""

#from sys import stderr
import numpy as np


def tot_branch_len(tree):
    tot_len = 0
    for node in tree.traverse():
        tot_len += node.dist
    return tot_len

def recurse_birth_events(node, root_dist=0, leaf_as_death=True):
    """Given a tree node, recursively yield (age of birth, number of newborns)
    for all descendants."""
    # not (not leaf_as_death and not children):
    if leaf_as_death or node.children:
        yield [root_dist, len(node.children) - 1]

        for child in node.children:
            for output in recurse_birth_events(child, root_dist+child.dist,
                                               leaf_as_death=leaf_as_death):
                yield output
    #else:
    #    # Stop at first leaf encountered
    #    yield (root_dist, 0)
    #    raise StopIteration


def get_ordered_birth_events(tree, compact=False, leaf_as_death=True,
                             stem=False, inplace=False):
    tree2 = tree #if inplace else tree.copy()
    # save the (birth_time, n_birth) list
    #birth_times = np.zeros((2, len(tree2))) # not filled when multifurc.
    root_dist = tree2.dist if stem else 0
    #tree2.add_feature('root_dist', root_dist)
    #birth_times[:, 0] = (root_dist, len(tree2.children) - 1)
    birth_times = [[0, 1]] # tree stem

    #root_dist
    #for node in tree2.traverse('levelorder'):
    #    try:
    #        parent_root_dist = node.up.root_dist
    #    except AttributeError:
    #        parent_root_dist = root_dist

    #    node.add_feature('root_dist', node.dist + parent_root_dist)
        #print(node.name, node.root_dist, len(node.children) - 1)
        #if node.children:
            #birth_times[:, i] = (node.root_dist, len(node.children) - 1)
    #    birth_times.append([node.root_dist, len(node.children) - 1])
            #i += 1
    birth_times += list(recurse_birth_events(tree, leaf_as_death=leaf_as_death))

    #birth_times = birth_times[:,:i]

    #print(birth_times)
    # Sort by birth_times.
    #birth_times[:, birth_times[0].argsort()]
    birth_times.sort()
    # Compact:
    #birth_times = [(age, b) for i, (age,b) in enumerate(birth_times]
    if compact:
        i = 1
        while i < len(birth_times):
            if birth_times[i][0] == birth_times[i-1][0]: # time interval zero
                birth_times[i-1][1] += birth_times.pop(i)[1]
            else:
                i += 1

    #births = []
    #for (a1,b1),(a2,b2) in zip(birth_times[:-1], birth_times[1:]):
    #    births.append()
    return birth_times


def get_LTT(tree, compact=False):
    """Lineages Through Time
    row 0: birth times;
    row 1: number of lineages after latest birth."""
    birth_events = np.array(get_ordered_birth_events(tree, compact)).T
    birth_events[1] = birth_events[1].cumsum() # count lineages
    return birth_events


def get_inter_birth_dist(tree):
    ltt = get_LTT(tree, compact=False)
    ltt_steps = np.vstack((ltt[0,1:] - ltt[0, :-1], ltt[1,:-1]))
    return ltt_steps


def get_cum_tot_branch_len(tree):
    ltt_steps = get_inter_birth_dist(tree)
    return ltt_steps.prod(axis=0).cumsum()


def get_root_dist(node, root_dist=0):
    if not node.children:
        yield root_dist

    else:
        for child in node.children:
            for descendant_dist in get_root_dist(child, root_dist+child.dist):
                yield descendant_dist


def is_ultrametric(tree, ndigits=2):
    #return len(set(round(tree.get_distance(leaf), ndigits) for leaf in tree)) == 1
    return len(set(round(d, ndigits) for d in get_root_dist(tree))) == 1
    #for node in tree.traverse('postorder'):


def div_gamma(tree):
    """Compute the gamma statistic of a phylogenetic tree:
    departure from a constant-rate birth-death process (Pybus & Harvey 2000)"""
    # check that branch lengths are meaningful
    # TODO: check that tree is bifurcating?
    if not is_ultrametric(tree):
        raise(ValueError("Tree: %r is not ultrametric" % tree.name))

    n = len(tree)
    T = get_cum_tot_branch_len(tree)[1:n]
    #print(T)
    try:
        Ttot = T[-1]
    except IndexError:
        #print("div_gamma WARNING: tree has only one node", file=stderr)
        return np.NaN

    try:
        gamma_denom = Ttot * np.sqrt(1. / (12 * (n-2)))
    except ZeroDivisionError:
        # Tree must have more than 2 leaves
        return np.NaN

    gamma_num =  T[:-1].mean() - Ttot/2.
    
    g = gamma_num / gamma_denom
    return g

    








