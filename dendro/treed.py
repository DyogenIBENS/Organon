#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import print_function

"""Like nw_ed from newick_utils: process each node, and execute a function
conditionally on some test"""

import sys
import re
import ete3
import argparse


VAR = {'d': 'node.dist',
       'c': 'len(node.children)',
       'C': 'node.children',
       'n': 'node.name',
       's': 'len(node)',
       'u': 'node.up',
       'l': 'node.is_leaf()',
       'r': 'node.is_root()',
       'f': 'node.features',
       'A': 'node.get_ancestors()'}
       #'_': 'node'

ACTIONS = {'w': 'print(node.write(format=outfmt, format_root_node=True))',
           'o': 'node.delete(prevent_nondicotomic=False, preserve_branch_length=True)',
           'p': 'print(node.name)',
           'a': 'node.add_child',
           'L': 'node.ladderize()',
           'd': 'node.detach()',
           'F': 'node.add_feature'}

EPILOG="""SHORTCUTS
Variables:
    """ + "\n    ".join("%s: %s" % item for item in VAR.items()) \
+ """
Actions:
    """ + "\n    ".join("%s: %s" % item for item in ACTIONS.items())


VAR_PATTERN = r'\b(' + '|'.join(VAR.keys()) + r')\b'
ACTION_PATTERN = r'\b(' + '|'.join(ACTIONS.keys()) + r')\b'

#print(VAR_PATTERN)
#print(ACTION_PATTERN)

def main(treefile, test, action, format, outfmt, strategy, is_leaf_fn,
         output=True):

    #print(re.sub(VAR_PATTERN, '{\g<0>}', test))
    #print(re.sub(ACTION_PATTERN, '{\g<0>}', action))
    
    test_str = re.sub(VAR_PATTERN, '{\g<0>}', test).format(**VAR)
    action_str = re.sub(ACTION_PATTERN, '{\g<0>}', action).format(**ACTIONS)

    #print(test_str)
    #print(action_str)

    tree = ete3.Tree(treefile, format=format)
    for node in tree.traverse(strategy, is_leaf_fn):
        if eval(test_str):
            exec(action_str)

    if output:
        print(tree.write(format=outfmt, format_root_node=True))#features=


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__, epilog=EPILOG, 
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
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
