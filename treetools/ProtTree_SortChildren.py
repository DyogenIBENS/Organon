#!/usr/bin/env python3
# -*- coding: utf-8 -*-


"""Rewrite a LibsDyogen ProteinTree with sorted node children (by id).

USAGE: ProtTree_SortChildren.py [infile] [outfile]

Default to stdin/stdout if no argument."""


from sys import argv, exit, stdin, stdout
from LibsDyogen import myProteinTree


def sort_children(tree):
    """Sort children based on their numerical id"""
    for node, nodedata in tree.data.items():
        nodedata.sort()


if __name__ == '__main__':
    if set(('-h', '--help')) & set(argv[1:]):
        print(__doc__)
        exit()
    elif len(argv) > 3:
        print(__doc__)
        exit(1)

    else:
        try:
            outfile = argv[2]
        except IndexError:
            outfile = stdout
        try:
            infile = argv[1]
        except IndexError:
            infile = stdin

    for tree in myProteinTree.loadTree(infile):
        sort_children(tree)
        tree.printTree(outfile)
