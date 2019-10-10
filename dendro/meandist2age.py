#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""1. Check that the tree is ultrametric (distances rounded off at 2 decimals).
2. Then average the distance to the root to annotate the age of each node."""

import argparse as ap
import ete3
from diversete.diversete import is_ultrametric
from dendro.bates import iter_distleaves
from dendro.any import ete3 as ete3_methods

get_items = ete3_methods.get_items


def main(inputtree, format=1, to_table=False):
    tree = ete3.Tree(inputtree, format=format)
    assert is_ultrametric(tree)
    
    for node in tree.traverse('postorder'):
        if not getattr(node, 'age', None):
            if node.children:
                age = float(sum(d for l,d in iter_distleaves(node,
                                                get_items))) / len(node)

            else:
                age = 0
            if to_table:
                print('%s\t%s' % (node.name, age))
            else:
                node.add_feature('age', age)

    if not to_table:
        print(tree.write(features=['age'], format=format,
              format_root_node=True))

if __name__ == '__main__':
    parser = ap.ArgumentParser(description=__doc__)
    parser.add_argument('inputtree')
    parser.add_argument('-f', '--format', type=int, default=1, 
                        help="input and output newick format")
    parser.add_argument('-t', '--to-table', action='store_true',
                        help='Output ages in tabulated list')
    
    args = parser.parse_args()
    main(**vars(args))
