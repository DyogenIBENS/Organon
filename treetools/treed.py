#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import print_function

"""Like nw_ed from newick_utils: process each node, and execute a function
conditionally on some test"""

import sys
import re
import ete3
import argparse


var = {'d': 'node.dist',
       'c': 'len(node.children)',
       'C': 'node.children',
       'n': 'node.name',
       's': 'len(node)',
       'u': 'node.up',
       'a': 'node.add_child',
       'l': 'node.is_leaf()',
       'r': 'node.is_root()',
       'L': 'node.ladderize()',
       'f': 'node.features',
       'A': 'node.get_ancestors()'}

actions = {'w': 'print(node.write(format=outfmt, format_root_node=True)); return',
           'o': 'node.delete(prevent_nondicotomic=False, preserve_branch_length=True)',
           'p': 'print(node.name)',
           'd': 'node.detach()'}

var_pattern = r'\b(' + '|'.join(var.keys()) + r')\b'
action_pattern = r'\b(' + '|'.join(actions.keys()) + r')\b'

#print(var_pattern)
#print(action_pattern)

def main(treefile, test, action, format, outfmt, strategy, is_leaf_fn,
         output=True):

    #print(re.sub(var_pattern, '{\g<0>}', test))
    #print(re.sub(action_pattern, '{\g<0>}', action))
    
    test_str = re.sub(var_pattern, '{\g<0>}', test).format(**var)
    action_str = re.sub(action_pattern, '{\g<0>}', action).format(**actions)

    #print(test_str)
    #print(action_str)

    tree = ete3.Tree(treefile, format=format)
    for node in tree.traverse(strategy, is_leaf_fn):
        if eval(test_str):
            exec(action_str)

    if output:
        print(tree.write(format=outfmt, format_root_node=True))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('treefile')
    parser.add_argument('test')
    parser.add_argument('action')
    parser.add_argument('-f', '--format', type=int, default=1, 
                        choices=[0,1,2,3,4,5,6,7,8,9,100],
                        help='input newick format [%(default)s]')
    parser.add_argument('-o', '--outfmt', type=int, default=1, 
                        choices=[-1,0,1,2,3,4,5,6,7,8,9,100],
                        help='output newick format [%(default)s]')
    parser.add_argument('-s', '--strategy', default='levelorder', 
                        choices=['levelorder', 'preorder', 'postorder'],
                        help='[%(default)s]')
    parser.add_argument('-l', '--is-leaf-fn')
    parser.add_argument('-n', dest='output', action='store_false',
                        help='Output the processed tree')
    # TODO: arguments begin and end.
    
    args = parser.parse_args()
    main(**vars(args))