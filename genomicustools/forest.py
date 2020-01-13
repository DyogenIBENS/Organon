#!/usr/bin/env python3
# -*- coding: utf-8 -*-


"""CLI Wrapper for all LibsDyogen.myProteinTree methods.

USAGE:

./forest.py <valid method> <proteinTreeFile> [<extra args>]

Example:

./forest.py flattenTree <proteinTreeFile> <PhylTreeFile> <recurs>
"""


from __future__ import print_function

from sys import argv, stderr, stdout, stdin
from collections import defaultdict
from UItools.autoCLI import build_cli_processor  # TODO: drop this dependency.

try:
    from LibsDyogen import myProteinTree, myPhylTree
except ImportError:
    try: 
        from LibsDyogen.utils import myProteinTree, myPhylTree
    except ImportError:
        from utils import myProteinTree, myPhylTree


def run(process, proteinTreeFile, converted_args):
    print("Give args: %s" % converted_args, file=stderr)

    if proteinTreeFile == '-':
        proteinTreeFile = stdin

    count_outputs = defaultdict(int)

    for tree in myProteinTree.loadTree(proteinTreeFile):
        r = process(tree, *converted_args)
        try:
            count_outputs[r] += 1
        except TypeError:
            count_outputs["unhashable"] += 1
            pass

        tree.printTree(stdout)

    print("Outputs counts:", count_outputs, file=stderr)


def main():
    # This is necessary for rebuildTree:
    myProteinTree.nextNodeID = int(1e8)
    process, converted_args = build_cli_processor(myProteinTree.ProteinTree,
                                    {'phyltree': myPhylTree.PhylogeneticTree},
                                    1,
                                    *argv[1:])
    run(process, argv[2], converted_args)


if __name__=='__main__':
    main()

