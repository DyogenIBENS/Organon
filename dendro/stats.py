#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Patristic matrix distances and the likes"""

import numpy as np
import itertools as it
from dendro.bates import rev_dfw_descendants
from dendro.any import myPhylTree as phyltree_methods, ete3 as ete3_methods

get_phylitems = phyltree_methods.get_items

def patristic_phyltree(phyltree, root=None, rootlength=0, point_loc=0):
    """Distances between every pair of branch.
    
    point_loc: fraction of branch length to slide the reference point:
        - 0 creates the patristic matrix between nodes;
        - 0.5 creates the patristic matrix between branch centers.
    """
    if root is None:
        root = phyltree.root
    if rootlength is None:
        rootlength = phyltree.rootlength if root==phyltree.root else phyltree.parent[root].distance

    branch_names =  [ch for n,ch in dfw_pairs_generalized(phyltree, get_phylchildren,
                                    queue=[(None, root)], include_root=True)]
    branch_indices = dict((br,i) for i,br in enumerate(branch_names))
    patristic_dists = np.zeros((len(branch_names), len(branch_names)))
    for (parent, dist), items in rev_dfw_descendants(phyltree, get_phylitems,
                                    queue=[(root, rootlength)]):
        #for ch, d in items:
        #    patristic_dists.loc[parent, ch] = dist/2. + d/2.
        p_i = branch_indices[parent]
        for (ch1,d1), (ch2,d2) in it.combinations(items, 2):
            ch1_i = branch_indices[ch1]
            ch2_i = branch_indices[ch2]
            if ch1_i < p_i: print('ch1 %s %d < parent %s %d' % (ch1, ch1_i, parent, p_i))
            for descendant1 in phyltree.allDescendants[ch1]:
                desc1_i = branch_indices[descendant1]
                patristic_dists[p_i, desc1_i] = (dist*(point_loc) + d1*(1-point_loc) + patristic_dists[ch1_i, desc1_i])
                for descendant2 in phyltree.allDescendants[ch2]:
                    desc2_i = branch_indices[descendant2]
                    patristic_dists[p_i, desc2_i] = (dist*point_loc + d2*(1-point_loc) + patristic_dists[ch2_i, desc2_i])
                    if desc2_i < desc1_i: print('desc2 %s %d < desc1 %s %d' % (descendant2, desc2_i, descendant1, desc1_i))
                    patristic_dists[desc1_i, desc2_i] = (d1*(1-point_loc) + d2*(1-point_loc)
                            + patristic_dists[ch1_i, desc1_i]
                            + patristic_dists[ch2_i, desc2_i])
                # Looks like this step could be done in one matrix product.
    #We filled only one half (older in rows, younger in columns).
    # Because rev_dfw iterates from young to old after visiting all descendants,
    # (young indices(column) > old indices(row)) it is the upper triangle.
    i_low, j_low = np.triu_indices(len(branch_names), k=1)
    patristic_dists[j_low, i_low] = patristic_dists[i_low, j_low]

    max_dist = patristic_dists.max()
    patristic_dists /= max_dist

    return patristic_dists, branch_names
