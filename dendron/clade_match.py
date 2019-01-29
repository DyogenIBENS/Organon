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

import signal, traceback


interrupted = 0
def debug_variables_at_sigint(sig, frame):
    global interrupted
    if interrupted == 0:
        interrupted = 1
    else:
        # Quit for real.
        traceback.print_stack(frame)
        sys.exit(1)

signal.signal(signal.SIGINT, debug_variables_at_sigint)


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
                    raise
        clade2sp.append((node.name, clade))
        clade2sp_dict[node.name] = clade
    return clade2sp


def match_clades(tree1, tree2):
    """"""
    clades1 = get_clade2sp(tree1)
    clades1_dict = dict(clades1)

    matching_clades = []
    unmatched_sp2 = set()
    inner1_matching_sp2 = set()  # inner node in tree1 that matches leaves in tree2
    clades2_names = set(n.name for n in tree2.traverse() if not n.is_leaf())

    # First match the leaves
    #for leaf2, spset2 in clades2:
    #    if len(spset2) > 1:
    #        # Not a leaf, skip.
    #        continue

    for leaf2 in tree2.iter_leaves():  # Not `iter_leaves` because I change them.
        leaf2name = leaf2.name
        if leaf2name in clades1_dict:  # If it is an existing leaf or clade.
            matching_clades.append((leaf2name, leaf2name))
        else:
            # Find a 'close-enough' clade name.
            lclade2 = leaf2name.lower()
            lclade2words = set(lclade2.split())
            for clade1, spset1 in clades1:
                lclade1 = clade1.lower()
                lclade1words = set(lclade1.split())
                #if lclade1 in lclade2 or lclade2 in lclade1:
                if lclade1words <= lclade2words or lclade1words >= lclade2words \
                        or ((lclade2 in lclade1 or lclade1 in lclade2)
                            and len(lclade2)>5 and len(lclade1)>5):
                    # Ignore this if clade1 already exists as another node
                    # in tree2:
                    if clade1 in clades2_names or lclade1 in clades2_names:
                        # Edit the species set from tree2:
                        clade2match = clade1 if clade1 in clades2_names else lclade1

                        # Update tree2 data if it's the first time we see this name.
                        if clade2match not in inner1_matching_sp2:
                            inner1_matching_sp2.add(clade2match)
                            node2match = tree2&clade2match
                            for descendant in node2match.iter_descendants():
                                matching_clades.append(('', descendant.name))
                                if descendant.is_leaf():
                                    unmatched_sp2.add(descendant.name)
                            node2match.children = []
                            # Also replace the value in the list of tuples
                            matching_clades.append((clade1, clade2match))
                    else:
                        matching_clades.append((clade1, leaf2name))
                        leaf2.name = clade1  # All clades2 will use this species name.
                    break

            else:
                matching_clades.append(('', leaf2name))
                unmatched_sp2.add(leaf2name)
                #clades1.insert(0, (None, set((leaf2name,))))

    species2 = set(tree2.get_leaf_names())
    clades2 = get_clade2sp(tree2)
    clades2_dict = dict(clades2)

    # Must be done from leaves to root, and memorize already matched descendants
    for clade1, spset1 in clades1:
        #if clade1.startswith('Heterocephalus glaber'):
        #    import ipdb; ipdb.set_trace()
        if len(spset1) == 1:
            ###TODO: keep clades with a single child/species.
            ###      merge this step in the more general step above. ! error with Cebidae
            sp1, = spset1
            if sp1 not in clades2_dict:  # values
                matching_clades.append((clade1, ''))

            # Skip because was already matched as a leaf.

        else:
            # Exclude species not at all present in tree2
            spset1 &= species2

            matching_cl = []
            matched_sp = set()
            unmatched_sp = set(spset1)
            
            while matched_sp != spset1:
                if interrupted:
                    raise KeyboardInterrupt("In while loop at:\n"
                            + "clade1: %s\n" % clade1
                            + "spset1: %s\n" % spset1
                            + "matched_sp: %s\n" % matched_sp
                            + "unmatched_sp: %s\n" % unmatched_sp
                            + "matching_cl: %s\n" % (matching_cl if len(matching_cl) < 100 else (matching_cl[:10], 'len=%d' % len(matching_cl)),)
                            + "best inclusion: %s, size=%d" % (best,len(best[1])))

                inclusions = [(clade2,
                               unmatched_sp & spset2.difference(unmatched_sp2))
                              for clade2, spset2 in clades2
                              if unmatched_sp >= spset2.difference(unmatched_sp2)]
                #TODO: if several equivalent matches, select most basal.
                try:
                    best = max(inclusions, key=lambda x: len(x[1]))
                except ValueError as err:
                    # Empty sequence
                    # There may be an intersection but no inclusions?
                    #best = max(clades2, key=lambda x: len(unmatched_sp & x[1]))
                    #print('WARNING: the best match is %r but it contains extra species:\n%s' \
                    #        % (best[0], best[1].difference(unmatched_sp)),
                    #       file=sys.stderr)
                    break

                    #continue
                if not best[1]:
                    # Largest inclusion is the empty set
                    break

                matching_cl.append('(%s)' % best[0] if '+' in best[0] else best[0])
                unmatched_sp.difference_update(best[1])
                matched_sp |= best[1]

                # Avoid infinite loops
                if not unmatched_sp:
                    break
            
            matching_cl_str = '+'.join(matching_cl)
            if len(matching_cl) >= 1:
                # add this (possibly polyphyletic) clade to the reference tree.
                clades2_dict[matching_cl_str] = matched_sp
                clades2.append((matching_cl_str, matched_sp))

            matching_clades.append((clade1, matching_cl_str))
    return matching_clades


def main_ete3(treefile1, treefile2, format=1):
    tree1 = ete3.Tree(treefile1, format=format)
    tree2 = ete3.Tree(treefile2, format=format)
    pass


def main(treefile1, treefile2, sort=False):
    tree1 = myPhylTree.PhylogeneticTree(treefile1).to_ete3()
    tree2 = myPhylTree.PhylogeneticTree(treefile2).to_ete3()
    print(treefile1 + '\t' + treefile2)
    matched_clades = match_clades(tree1, tree2)
    if sort:
        matched_clades.sort()
    for matching in matched_clades:
        print('\t'.join(matching))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('treefile1')
    parser.add_argument('treefile2')
    parser.add_argument('-s', '--sort', action='store_true')
    
    args = parser.parse_args()
    main(**vars(args))
