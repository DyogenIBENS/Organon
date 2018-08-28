#!/usr/bin/env python3


"""Convert LibsDyogen treeforest to Newick format (with special tags).

USAGE:
    
    ./genomicus2newick.py [forestfile]
"""

from __future__ import print_function


from sys import argv, exit, stderr, stdout, setrecursionlimit
import LibsDyogen.myProteinTree as ProteinTree


if len(argv) != 2:
    print(__doc__, file=stderr)
    exit(1)
elif argv[1] in ('-h', '--help'):
    print(__doc__, file=stderr)
    exit()

forestfile=argv[1]

setrecursionlimit(20000)

for tree in ProteinTree.loadTree(forestfile):
    #for node, children in tree.data.items():
    #    print(node, children)
    tree.printDyogenNewick(stdout)

