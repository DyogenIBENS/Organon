#!/usr/bin/env python3
# -*- coding: utf-8 -*-


"""Compare `tree1` to a reference `tree2`:
- for each clade in `tree1`, find the corresponding clade in `tree2`;
- if not found, find the smallest number of monophyletic clades containing all the listed
  species, and only them.
"""

import sys
import argparse
#from dendron.climber import dfw_descendants_generalized
from itertools import zip_longest, product
import logging
logger = logging.getLogger(__name__)

import dendron.parsers as treeparsers
import dendron.converters as treeconverters

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
        nodename = node.name
        # Ensure this name is not a duplicate
        dup = 0
        while nodename in clade2sp_dict:
            dup += 1
            nodename = '%s_%d' % (node.name, dup)
        node.name = nodename  # For later collecting as child.

        clade = set()
        if node.is_leaf():
            clade.add(nodename)
        else:
            for child in node.children:
                try:
                    clade.update(clade2sp_dict[child.name])
                except KeyError as err:
                    print('KeyError: Current clade2sp_dict %s' % clade2sp_dict,
                          file=sys.stderr)
                    err.args += ('at node %r' % nodename,)
                    raise
        clade2sp.append((nodename, clade))
        clade2sp_dict[nodename] = clade
    return clade2sp


def name_similarity(name1, name2):
    lclade1 = str(name1).lower()
    lclade1words = set(lclade1.split())
    lclade2 = str(name2).lower()
    lclade2words = set(lclade2.split())
    l1 = len(lclade1)
    l2 = len(lclade2)
    minlen = min(l1, l2)
    if lclade1words <= lclade2words or lclade1words >= lclade2words:
        return len(lclade1words & lclade2words)
    elif lclade1words & lclade2words:
        score = len(lclade1words & lclade2words)
        words1 = lclade1words - lclade2words
        words2 = lclade2words - lclade1words
        for w1, w2 in product(words1, words2):
            if w1 in w2 or w2 in w1:
                score += float(min(len(w1), len(w2))) / max(len(w1), len(w2))
        return score
    elif (lclade2 in lclade1 or lclade1 in lclade2) and minlen>5:
        return float(minlen) / max(l1, l2)
    else:
        return 0


def match_clades(tree1, tree2, exact=False):
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
        elif not exact:
            # Find a 'close-enough' clade name.
            #if 'Heterocephalus glaber' in leaf2name: import ipdb; ipdb.set_trace()
            similars = sorted([(clade1, name_similarity(leaf2name, clade1))
                               for clade1,_ in clades1],
                              key=lambda x: x[1],
                              reverse=True)
            most_sim, sim = similars[0]
            if sim and len(similars)>1 and sim > similars[1][1]:  # Unambiguous name match
                logger.info('Fuzzy matching %r ~ %r (score=%g)',
                            leaf2name, most_sim, sim)
                lmost_sim = str(most_sim).lower()

                # Ignore this if clade1 already exists as another node
                # in tree2:
                if most_sim in clades2_names or lmost_sim in clades2_names:
                    # Edit the species set from tree2:
                    clade2match = most_sim if most_sim in clades2_names else lmost_sim

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
                        matching_clades.append((most_sim, clade2match))
                else:
                    matching_clades.append((most_sim, leaf2name))
                    leaf2.name = most_sim  # All clades2 will use this species name.
            else:
                if sim:
                    logger.warning('Ambiguous fuzzy matching %r ~ %s (score=%g)',
                                   leaf2name,
                                   [name for name,score in similars if score>=sim],
                                   sim)
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

                matching_cl.append('(%s)' % best[0] if '+' in str(best[0]) else str(best[0]))
                unmatched_sp.difference_update(best[1])
                matched_sp |= best[1]

                # Avoid infinite loops
                if not unmatched_sp:
                    break
            
            matching_cl_str = '+'.join(sorted(matching_cl))  # sort line for reproducibility
            if len(matching_cl) >= 1:
                # add this (possibly polyphyletic) clade to the reference tree.
                clades2_dict[matching_cl_str] = matched_sp
                clades2.append((matching_cl_str, matched_sp))

            matching_clades.append((clade1, matching_cl_str))
    return matching_clades


def main(treefile1, treefile2, exact=False, sort=False, parser='PhylTree',
         output='<=>!'):
    try:
        parse_trees = treeparsers.parserchoice[parser]
        to_ete3 = treeconverters.converterchoice[parser]['ete3']

        def iter_trees(treefile, *args, **kwargs):
            for tree in parse_trees(treefile, *args, **kwargs):
                yield to_ete3(tree)

    except KeyError:
        raise ValueError('Bad parser value', parser)

    parser_kwargs = {'format': 1} if parser.lower() == 'ete3' else {}

    print(treefile1 + '\t' + treefile2)
    for tree1, tree2 in zip_longest(iter_trees(treefile1, **parser_kwargs),
                                    iter_trees(treefile2, **parser_kwargs)):
        if tree1 is None and tree2 is not None:
            logger.error("More trees2 than trees1. Ignoring.")
            break
        elif tree2 is None:
            logger.error("More trees1 than trees2. Ignoring.")
            #TODO: recycle if only one tree2
            break

        matched_clades = match_clades(tree1, tree2, exact)
        if sort:
            matched_clades.sort()
        for matching in matched_clades:
            matchtype = ('<' if not matching[1] else
                         '>' if not matching[0] else
                         '=' if matching[0] == matching[1] else
                         '!')
            if matchtype in output:
                print('%s\t%s' % matching)


if __name__ == '__main__':
    logging.basicConfig()
    logger.setLevel(logging.INFO)
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('treefile1')
    parser.add_argument('treefile2')
    parser.add_argument('-e', '--exact', action='store_true',
                        help='Do not attempt fuzzy matching of not found leaves.')
    parser.add_argument('-s', '--sort', action='store_true')
    parser.add_argument('-p', '--parser', default='PhylTree',
                        help=('Choices: "ete3", "PhylTree", "ProtTree" (case '
                              'insensitive) [%(default)s]'))
    parser.add_argument('-o', '--output', default='<=>!',  # '~' for fuzzy
                        help=('Matches to output:\n'
                              '`=` identical,\n'
                              '`!` changed,\n'
                              '`<` only in left tree,\n'
                              '`>` only in right tree.'))
    
    args = parser.parse_args()
    main(**vars(args))
