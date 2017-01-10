#!/usr/bin/env python2.7

"""USAGE:
    ./reshape_tree.py [-h] <input_tree> <output_tree>

OPTIONS:
    -h: print this help and exit

DESCRIPTION:
    produces reshaped newick tree compatible with codeml.
    Reshaping means removing the inner node labels, and removing nodes that
    start from the root and have only one child (codeml fails on such a tree).
    Requires the python module 'ete3'.
"""

import sys
import ete3


if __name__=='__main__':
    if len(sys.argv) != 3:
        print __doc__
        if len(sys.argv) == 2 and sys.argv[1] in ('-h', '--help'):
            sys.exit(0)
        sys.exit(1)

    input_tree = sys.argv[1]
    output_tree = sys.argv[2]

    tree = ete3.Tree(input_tree, format=1)
    for node in tree.traverse():
        outtree = node
        if len(node.children) > 1:
            break
    # The format compatible with codeml is format=5
    outtree.write(outfile=output_tree, format=5)
