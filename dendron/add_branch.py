#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from sys import stdin
import argparse
import ete3
import logging
logger = logging.getLogger(__name__)
logging.basicConfig(format="%(levelname)s:%(message)s")


def fix_negative_lengths(node):
    # Recursively transfer the negative length to descendants, until it creates
    # non-negative length.
    count = 0
    if node.dist < 0:
        logger.debug("Fixing length of %r: %g -> 0", node.name, node.dist)
        count += 1
        if node.is_leaf():
            logger.error("Can't fix negative length in leaf %r", node.name)
            return 0
        for child in node.children:
            logger.debug("Updating length of %r: %g -> %g", child.name,
                                                            child.dist,
                                                            child.dist + node.dist)
            child.dist += node.dist
            count += fix_negative_lengths(child)
        node.dist = 0
    return count



def main(infile, subtreesfile=None, check_ultrametricity=-1,
         format=1, quoted_node_names=False):
    if check_ultrametricity < 0:
        is_ultrametric = lambda tree, thresh: True
    else:
        from dendron.climber import iter_distleaves
        get_data = lambda tree,dat: [(ch, ch.dist) for ch in dat[0].children]
        def is_ultrametric(tree, thresh=0.01):
            leaf_dists = [d for _,d in iter_distleaves(tree, get_data)]
            return max(leaf_dists) - min(leaf_dists) < thresh

    tree = ete3.Tree(infile, format=format)
    if not is_ultrametric(tree, check_ultrametricity):
        logger.warning("Requested ultrametricity check but the input tree is not ultrametric.")

    stream = stdin if subtreesfile is None else open(subtreesfile)
    subtrees = [s + ';' for s in stream.read().rstrip('; \t\n').split(';')]
    if subtreesfile is None:
        stream.close()

    for subtreetxt in subtrees:
        newsubtree = ete3.Tree(subtreetxt, format=format, quoted_node_names=quoted_node_names)
        logger.info("Inserting " + newsubtree.name)
        # First leaf of newsubtree should be an existing node in the main tree.
        anchor_node = newsubtree.get_leaves()[0]
        assert anchor_node == newsubtree.children[0], \
            "Incorrect subtree specification: anchor_node should be a leaf. %r VS\n%s" % \
                (anchor_node.name, newsubtree.get_ascii())
        try:
            orig_node = tree.search_nodes(name=anchor_node.name)[0]
        except IndexError:
            raise LookupError('Node %r not found in the source tree' % anchor_node.name)

        if check_ultrametricity > 0:
            age_anchor = orig_node.get_farthest_leaf()[1]
            age_via_anchor = age_anchor + anchor_node.dist
            inserted_leafdists = [ch.dist + d
                                  for ch in newsubtree.children[1:]
                                  for _,d in iter_distleaves(ch, get_data)]
            if max(inserted_leafdists) - min(inserted_leafdists) >= check_ultrametricity:
                logger.error("Requested ultrametricity but an input subtree is not ultrametric.")
            age_via_inserted = max(inserted_leafdists)
            assert age_via_anchor - age_via_inserted < check_ultrametricity, \
                "The combination of inserted leaf ages and anchor age cannot"\
                " be made ultrametric, at %r:\n" % newsubtree.name \
                + "%r(age = %g):dist=%g\n" % (anchor_node.name, age_anchor, anchor_node.dist) \
                + "\n".join("%r(age = %g):dist=%g" % (ch.name, age_via_inserted - ch.dist, ch.dist)
                            for ch in newsubtree.children[1:])
                        
        parent = orig_node.up
        orig_dist = orig_node.dist
        anchor_dist = anchor_node.dist
        inserted_dists_diffs = [nch.dist - anchor_dist for nch in newsubtree.children]

        orig_node.detach()
        new_dist = orig_dist - anchor_dist

        parent.add_child(child=newsubtree, dist=new_dist)
        parent.swap_children()

        #for ref_child in orig_node.children:
        #    anchor_node.add_child(child=ref_child)
        #Plus all other features

        anchor_node.detach()
        newsubtree.add_child(child=orig_node, dist=anchor_dist)
        newsubtree.swap_children()

        if new_dist < 0:
            logger.warning("New branch to %r longer than the original (%g > %g), fixing.",
                           anchor_node.name, anchor_dist, orig_dist)
            ## ALL CHILDREN of newsubtree must be shortened!
            nfixes = fix_negative_lengths(newsubtree)
            logger.info("Fixed %d branch(es).", nfixes)

        assert newsubtree.dist + orig_node.dist == orig_dist,\
                "%r:%g + %r:%g != %g" % (newsubtree.name, newsubtree.dist,
                                         orig_node.name, orig_node.dist,
                                         orig_dist)
        assert all((nch.dist - orig_node.dist - inserted_diff < 1e-12)
                   for nch, inserted_diff \
                   in zip(newsubtree.children, inserted_dists_diffs)), \
                   ("Changed children dists differences, at %s.\n" % newsubtree.name
                    + "\n".join("%s: %s -> %s" % (nch.name,
                                                  nch.dist - orig_node.dist,
                                                  inserted_diff)
                                for nch, inserted_diff in
                                zip(newsubtree.children, inserted_dists_diffs)))
        assert is_ultrametric(parent, check_ultrametricity)

    assert is_ultrametric(tree, check_ultrametricity)
    logger.debug("Ultrametric = %s:\ndist to %s = %s;\ndist to %s = %s",
                 is_ultrametric(tree, check_ultrametricity),
                 *(tree.get_closest_leaf() + tree.get_farthest_leaf()))

    print(tree.write(format=format, format_root_node=True))

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('infile')
    parser.add_argument('subtreesfile', nargs='?',
                        help='subtrees to insert. The root should be the new ' \
                             'branching node, and the top leaf should be an ' \
                             'existing node.')
    parser.add_argument('-u', '--check-ultrametricity', type=float, default=-1,
                        help='Threshold to check for ultrametricity.' \
                             'Negative indicates no check. [%(default)s]')
    parser.add_argument('-f', '--format', type=int, default=1)
    parser.add_argument('-q', '--quoted-node-names', '--quoted', action='store_true')
    parser.add_argument('-v', '--verbose', action='count')
    
    args = parser.parse_args()
    if args.verbose == 1:
        logger.setLevel(logging.INFO)
    elif args.verbose == 2:
        logger.setLevel(logging.DEBUG)

    logger.debug(str(args) + str(type(args.check_ultrametricity)))
    main(args.infile, args.subtreesfile, args.check_ultrametricity,
         args.format, args.quoted_node_names)
