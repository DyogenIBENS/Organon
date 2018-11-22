#!/usr/bin/env python3


"""Edit chosen nodes from the tree: keep only their closest leaf."""

from sys import stdin, stdout
import argparse as ap
from LibsDyogen import myProteinTree
from dendron.climber import iter_leaves, iter_leaf_paths
import logging
logger = logging.getLogger(__name__)
logger.setLevel(logger.INFO)
#logging.basicConfig(format="%(levelname)s:%(message)s", level=logger.INFO)


def get_data(tree, datum):
    return tree.data.get(datum[0], [])

def get_children(tree, child):
    return [x for x, _ in tree.data.get(child, [])]


def iter_distleaves(tree, root):
    for leafpath in iter_leaf_paths(tree, get_data, [(root,0)]):
        leaf = leafpath[-1][0]
        leafdist = sum(dist for node, dist in leafpath[1:])
        yield leaf, leafdist


def keep_closest_leaf(tree, node):
    leafdists = sorted(iter_distleaves(tree, node), key=lambda datum: datum[1])
    # This does not keep intermediate nodes.
    tree.data[node] = [leafdists[0]]


def edit(badnodes, proteintrees, yield_unchanged=False):
    """Keep only one child, the one leading to the closest leaves."""
    
    n_edits = 0
    n_included = 0  # Number of nodes not edited, because already in a larger
                    # edited subtree
    n_morethan2 = 0  # Number of nodes with >2 leaves

    for tree in proteintrees:
        tree_badnodes = badnodes.intersection(tree.data)
        if not tree_badnodes:
            if yield_unchanged:
                yield tree
            continue

        # First the bad nodes must be sorted according to their rootwardness.
        # (in order to edit only the most basal in case of inclusion, and
        # ignore the included ones)
        badnode_leaves = [(badnode,
                           set(iter_leaves(tree, get_children, [badnode])))
                          for badnode in tree_badnodes]

        edited_leafsets = []
        # Indices in the above list where the leafset size decreases
        larger_size_i = 0
        # Current size
        size = len(badnode_leaves[0][1])

        for badnode, leaves in sorted(badnode_leaves, key=lambda x: len(x[1]), reverse=True):

            # Check leaf number: then only check intersection with strictly larger sets.
            new_size = len(leaves)
            if new_size < size:
                #logger.debug("Size decrease!")
                larger_size_i = len(edited_leafsets)
                size = new_size

            # Edit only if no ancestral node has already been edited.
            if not any((leaves & edited_l) for edited_l in edited_leafsets[:larger_size_i]):
                edited_leafsets.append(leaves)
                keep_closest_leaf(tree, badnode)
                n_edits += 1
            else:
                n_included += 1
        
        assert size == 2
        n_morethan2 += larger_size_i

        badnodes.difference_update(tree_badnodes)
        yield tree

    if badnodes:
        logger.warning("%d nodes not found: %s...", len(badnodes),
                        ' '.join(str(n) for n in list(badnodes)[:5]))
    logger.info("%d edited nodes\n%d implicitely edited\n%d with >2 leaves",
                 n_edits, n_included, n_morethan2)


def main(badnodelistfile, forestfile, badnode_col=0):
    with open(badnodelistfile) as f:
        header_line = next(f)
        badnodes = set(int(line.rstrip().split()[0]) for line in f)

    logger.info('%d nodes to remove.', len(badnodes))

    for tree in edit(badnodes, myProteinTree.loadTree(forestfile)):
        tree.printTree(stdout)


if __name__ == '__main__':
    parser = ap.ArgumentParser(description=__doc__)
    parser.add_argument('badnodelistfile')
    parser.add_argument('forestfile', nargs='?', default=stdin)

    args = parser.parse_args()
    main(**vars(args))

