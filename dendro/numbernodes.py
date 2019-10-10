#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import ete3
import argparse
from dendro.parsers import read_multinewick


def main(inputfile, order='preorder', start=None, attr='name'):
    #exclude_leaves=True
    if inputfile == '-':
        inputlines = sys.stdin.readlines()
    else:
        with open(inputfile) as f:
            inputlines = f.readlines()

    for newick in read_multinewick(inputlines):
        tree = ete3.Tree(newick, format=1)
        i = len(tree)+1 if start is None else start
        for node in tree.traverse(order):
            if not node.is_leaf():
                node.add_feature(attr, i)
                i += 1
        out_features = None if attr=='name' else [attr]
        print(tree.write(features=out_features, format=1, format_root_node=True))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('inputfile', nargs='?', default='-',
                        help='(multi-)newick tree [stdin]')
    parser.add_argument('-o', '--order', default='preorder',
                        choices=['preorder', 'levelorder', 'postorder'],
                        help='Tree traversal order [%(default)s]')
    parser.add_argument('-s', '--start', metavar='N', type=int,
                        help='First node number, or the number of leaves if None.')
    parser.add_argument('-a', '--attr', default='name',
                        help='attribute name of the number [%(default)s]')
    
    args = parser.parse_args()
    main(**vars(args))
