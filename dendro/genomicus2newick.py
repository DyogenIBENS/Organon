#!/usr/bin/env python3


"""Convert LibsDyogen treeforest to **Newick** format (with special tags).
"""

from __future__ import print_function


from sys import stdin, stdout, setrecursionlimit
import argparse as ap
import LibsDyogen.myProteinTree as ProteinTree


if __name__ == '__main__':
    parser = ap.ArgumentParser(description=__doc__)
    parser.add_argument('forestfile', nargs='?', default=stdin,
                        type=ap.FileType('r'))
    parser.add_argument('outfile', nargs='?', default=stdout,
                        type=ap.FileType('w'))
    args = parser.parse_args()

    setrecursionlimit(20000)

    for tree in ProteinTree.loadTree(args.forestfile):
        #for node, children in tree.data.items():
        #    print(node, children)
        tree.printDyogenNewick(args.outfile)

