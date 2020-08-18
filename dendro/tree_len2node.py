#!/usr/bin/env python3
# -*- coding: utf-8 -*-


"""Break branches by inserting single-child nodes.

Can be used to weight branches more in dollocorr.
"""

import argparse as ap
import ete3


def break_branches(tree, step=1):
    for n in tree.traverse():
        for i, ch in enumerate(n.get_children()):
            breaks = ch.dist // step
            lastdist = ch.dist % step
            ch.dist = lastdist
            for j in range(int(breaks)):
                new = ete3.TreeNode(dist=step)
                ## set this new node as the parent of the current one.
                #ch.up = new
                ## take one step back (rootward)
                #ch = new
                new.add_child(ch.detach())
                ch = new
            if breaks:
                n.add_child(ch)
                #n.children[i] = ch
                #ch.up = n
            
    # Skip transforming the root node.


def main():
    parser = ap.ArgumentParser(description=__doc__)
    parser.add_argument('intree')
    parser.add_argument('-s', '--step', type=int, default=1,
                        help="Breaking distance")
    
    args = parser.parse_args()

    tree = ete3.Tree(args.intree, format=1)
    break_branches(tree, args.step)
    print(tree.write(format=1, format_root_node=True))


if __name__ == '__main__':
    main()
