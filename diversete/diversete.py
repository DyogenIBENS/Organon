#!/usr/bin/env python3
# -*- coding: utf-8 -*-


"""Utility functions for ete3 objects, to perform diversification analyses:
    
    - diversification net rate
    - constant rate BD fit
    - gamma statistic"""


import numpy as np
import logging
logger = logging.getLogger(__name__)
#logging.basicConfig("%(levelname)s:%(funcName)s:%(message)s")


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
        logger.warning("Tree has only one node")
        return np.NaN

    try:
        gamma_denom = Ttot * np.sqrt(1. / (12 * (n-2)))
    except ZeroDivisionError:
        # Tree must have more than 2 leaves
        return np.NaN

    gamma_num =  T[:-1].mean() - Ttot/2.
    
    g = gamma_num / gamma_denom
    return g


#from dendron.climber import dfw_pairs
logging.basicConfig()
logger.setLevel(logging.DEBUG)

def speciation_numbers_on_selected_lineages(fulltree, reftree, ref2full=None):
    reftree = reftree.copy()

    if ref2full is None:
        ref2full = lambda x: x.replace(' ', '_')

    for node in reftree.get_descendants('postorder'):
        nodename = ref2full(node.name)
        if node.is_leaf():
            if not 'age' in node.features:
                node.add_feature('age', 0)

            matched = fulltree.get_leaves_by_name(nodename)
            if len(matched) == 0:
                logger.warning('Leaf %r not found (detach)', nodename)
                node.detach()
                break
            elif len(matched) > 1:
                logger.warning('Multiple matches for %r (take first)', nodename)
            matched = matched[0]
            age = node.age
        else:
            matched = node.matched
            age = matched.age

        logger.debug('* Ref node %s (match: %s)', nodename,matched.name)
        
        parent = node.up
        pname = ref2full(parent.name)
        logger.debug('  parent:  %s', pname)

        branch_speciations = 0
        branch_sistersizes = []
        branching_times = []

        while matched.up and matched.up.name != pname:

            sisters = matched.get_sisters()
            logger.debug(' - speciation: %s -> %s', matched.up.name,
                         [s.name for s in sisters])

            branch_speciations += len(sisters)
            branch_sistersizes.extend(len(sis) for sis in sisters)
            branching_times.extend([age - matched.dist]*len(sisters))

            age -= matched.dist
            matched = matched.up

        node.add_features(speciations=branch_speciations,
                          sistersizes=branch_sistersizes,
                          branching_times=branching_times)
        
        if not matched and not node.up.is_root():
            logger.warning('Reached the root in fulltree before reftree (%s)',
                           node.name)
            break
        
        age -= matched.dist
        matched = matched.up
        logger.debug('  parent match: %s', matched.name)

        matched.add_feature('age', age)
        parent.add_features(matched=matched,
                           age=(node.age - node.dist))

    return reftree


# Example execution
import os.path as op
from LibsDyogen import myPhylTree

def example_speciation_numbers(
        phyltreefile="~/ws2/DUPLI_data93/PhylTree.TimeTree201901.Ensembl93-like.goodQual.nwk",
        fulltreefile="~/ws2/databases/timetree/Opisthokonta_species.taxa.nwk",
        root='Primates'):

    phyltree = myPhylTree.PhylogeneticTree(op.expanduser(phyltreefile))
    ptree = phyltree.to_ete3()
    rtree = ptree&root

    fulltree = ete3.Tree(fulltreefile, format=1)

    fixmatches = {'Papionini': '58293', #'58249',
                  'HomoPan': '58326'}#,

    ref2full = lambda x: fixmatches.get(x, '_'.join(x.split()[:2]))

    nrtree = speciation_numbers_on_selected_lineages(fulltree, rtree, ref2full)

    def mylayout(node):
        sistersizes = getattr(node, 'sistersizes', ())
        ete3.add_face_to_node(
                ete3.TextFace(' '.join(str(s) for s in reversed(sistersizes))),
                node, 0, position='branch-top')
        #spe_number_face

