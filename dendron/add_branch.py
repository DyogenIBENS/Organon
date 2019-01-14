#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from sys import stdin
import argparse
import ete3
import logging
logger = logging.getLogger(__name__)
#logging.basicConfig(format="%(levelname)s:%(message)s")


def main(infile, subtreesfile=None, format=1, quoted_node_names=False):
    tree = ete3.Tree(infile, format=format)

    stream = stdin if subtreesfile is None else open(subtreesfile)
    subtrees = [s + ';' for s in stream.read().split(';')]
    if subtreesfile is None:
        stream.close()

    for subtreetxt in subtrees:
        newsubtree = ete3.Tree(subtreetxt, format=format, quoted_node_names=quoted_node_names)
        # First leaf of newsubtree should be an existing node in the main tree.
        ref_leaf = newsubtree.get_leaves()[0]
        assert ref_leaf == newsubtree.children[0], "incorrect subtree specification"
        try:
            ref_node = tree.search_nodes(name=ref_leaf.name)[0]
        except IndexError:
            raise ValueError('Node %r not found in source tree' % ref_leaf.name)

        parent = ref_node.up
        ref_node.detach()
        new_node_dist = ref_node.dist - ref_leaf.dist
        if new_node_dist < 0:
            logger.warning("New branch to %r longer than the original, reducing.",
                            ref_leaf.name)
            ref_leaf.dist -= new_node_dist
            new_node_dist = 0
        parent.add_child(child=newsubtree, dist=new_node_dist)

        #for ref_child in ref_node.children:
        #    ref_leaf.add_child(child=ref_child)
        #Plus all other features

        # Do not preserve the provided order of children of the subtree
        ref_leaf.detach()
        newsubtree.add_child(child=ref_node, dist=ref_leaf.dist)

    print(tree.write(format=format, format_root_node=True))

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('infile')
    parser.add_argument('subtreesfile', nargs='?',
                        help='subtrees to insert. The root should be the new ' \
                             'branching node, and the top leaf should be an ' \
                             'existing node.')
    parser.add_argument('-f', '--format', type=int, default=1)
    parser.add_argument('-q', '--quoted-node-names', '--quoted', action='store_true')
    
    args = parser.parse_args()
    print(args)
    main(args.infile, args.subtreesfile, args.format, args.quoted_node_names)
