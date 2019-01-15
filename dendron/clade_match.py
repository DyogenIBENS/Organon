#!/usr/bin/env python3
# -*- coding: utf-8 -*-


"""Compare `tree1` to a reference `tree2`:
- for each clade in `tree1`, find the corresponding clade in `tree2`;
- if not found, find the smallest number of monophyletic clades containing all the listed
  species, and only them.
"""

import sys
import argparse
import ete3
from LibsDyogen import myPhylTree
from collections import defaultdict


def get_clade2sp(tree):
    clade2sp = []
    clade2sp_dict = {}
    for node in tree.traverse('postorder'):
        clade = set()
        if node.is_leaf():
            clade.add(node.name)
        else:
            for child in node.children:
                try:
                    clade.update(clade2sp_dict[child.name])
                except KeyError as err:
                    print('KeyError: Current clade2sp_dict %s' % clade2sp_dict,
                          file=sys.stderr)
                    err.args += ('at node %r' % node.name,)
        clade2sp.append((node.name, clade))
        clade2sp_dict[node.name] = clade
    return clade2sp


def match_clades(tree1, tree2):
    """NOT WORKING YET"""
    clades1 = get_clade2sp(tree1)
    clades1_dict = dict(clades1)

    matching_clades = []
    unmatched_sp2 = set()

    for leaf2 in tree2.iter_leaves():
        if leaf2.name in clades1_dict:  # If it is an existing leaf or clade.
            matching_clades.append((leaf2.name, leaf2.name))
        else:
            # Find a 'close-enough' clade name.
            lclade2 = leaf2.name.lower()
            for clade1, spset1 in clades1:
                lclade1 = clade1.lower()
                if lclade1 in lclade2 or lclade2 in lclade1:
                    matching_clades.append((clade1, leaf2.name))
                    leaf2.name = clade1  # All clades2 will use this species name.
                    break
            else:
                matching_clades.append(('', leaf2.name))
                unmatched_sp2.add(leaf2.name)
                #clades1.insert(0, (None, set((leaf2.name,))))

    species2 = set(tree2.get_leaf_names())
    clades2 = get_clade2sp(tree2)
    clades2_dict = dict(clades2)

    # Must be done from leaves to root, and memorize already matched descendants
    for clade1, spset1 in clades1:
        if len(spset1) == 1:
            ###TODO: keep clades with a single child/species.
            ###      merge this step in the more general step above. ! error with Cebidae
            sp1, = spset1
            if sp1 not in clades2_dict:  # values
                matching_clades.append((clade1, ''))

        else:
            # Exclude species not at all present in tree2
            spset1 &= species2

            matching_cl = []
            matched_sp = set()
            unmatched_sp = set(spset1)
            
            while matched_sp != spset1:
                includedness = [(clade2, unmatched_sp & spset2.difference(unmatched_sp2))
                                for clade2, spset2 in clades2
                                if unmatched_sp >= spset2.difference(unmatched_sp2)]
                #TODO: if several equivalent matches, select most basal.
                try:
                    best = max(includedness, key=lambda x: len(x[1]))
                except ValueError as err:
                    # Empty sequence
                    # There may be an intersection but no includedness?
                    #best = max(clades2, key=lambda x: len(unmatched_sp & x[1]))
                    #print('WARNING: the best match is %r but it contains extra species:\n%s' \
                    #        % (best[0], best[1].difference(unmatched_sp)),
                    #       file=sys.stderr)
                    break

                    #continue
                matching_cl.append(best[0])
                unmatched_sp.difference_update(best[1])
                matched_sp |= best[1]
            
            if len(matching_cl) >= 1:
                # add this (possibly polyphyletic) clade to the reference tree.
                clades2_dict['+'.join(matching_cl)] = matched_sp
                clades2.append(('+'.join(matching_cl), matched_sp))

            matching_clades.append((clade1, '+'.join(matching_cl)))
    return sorted(matching_clades)


def main_ete3(treefile1, treefile2, format=1):
    tree1 = ete3.Tree(treefile1, format=format)
    tree2 = ete3.Tree(treefile2, format=format)
    pass


def main(treefile1, treefile2):
    tree1 = myPhylTree.PhylogeneticTree(treefile1).to_ete3()
    tree2 = myPhylTree.PhylogeneticTree(treefile2).to_ete3()
    print(treefile1 + '\t' + treefile2)
    for matching in match_clades(tree1, tree2):
        print('\t'.join(matching))



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('treefile1')
    parser.add_argument('treefile2')
    
    args = parser.parse_args()
    main(**vars(args))
