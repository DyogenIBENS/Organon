#!/usr/bin/env python3
# -*- coding: utf-8 -*-


"""Rotate branches according to the descendant leaf names.

Input/output trees are in Dyogen ProteinTree format.
"""


from sys import stdin, stdout, stderr, setrecursionlimit
import argparse

from sorting import leaf_sort
from LibsDyogen import myProteinTree, myFile #, myTools


setrecursionlimit(20000)


def ProteinTree_getchildren(tree, node_and_dist):
    #return [child for (child,dist) in tree.data.get(node, [])]
    return tree.data.get(node_and_dist[0], [])


def ProteinTree_getnodeattr(tree, node_and_dist, attrname="gene_name"):
    try:
        return tree.info[node_and_dist[0]][attrname]
    except KeyError:
        print(node_and_dist, tree.info[node_and_dist[0]], file=stderr)
        raise


def ProteinTree_getId(tree, node_and_dist):
    return node_and_dist[0]


def ProteinTree_assignchildren(tree, node_and_dist, children):
    tree.data[node_and_dist[0]] = children


def ProteinTree_LeafSort(tree, ProteinTree_getattribute=ProteinTree_getId):
    leaf_sort(tree, (tree.root, 0.),
              ProteinTree_getchildren, ProteinTree_assignchildren,
              get_attribute=ProteinTree_getattribute)


def main(proteinTree, outFile, sortAttr=None):
    if sortAttr:
        get_attribute = lambda *args: ProteinTree_getnodeattr(*args, attrname=sortAttr)
    else:
        get_attribute = ProteinTree_getId

    with myFile.openFile(outFile, 'w') as out:
        for tree in myProteinTree.loadTree(proteinTree):
            ProteinTree_LeafSort(tree, get_attribute)
            tree.printTree(out)


if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("proteinTree", nargs='?', default=stdin,
                        help="Default to stdin")
    parser.add_argument("outFile", nargs='?', default=stdout,
                        help="Default to stdout")
    parser.add_argument("-a", "--sortAttr", help="Leaf attribute to select:\n"\
                        "can be 'gene_name', 'protein_name', 'taxon_name', ...")

    args = vars(parser.parse_args())

    #args = myTools.checkArgs(
    #        [("proteinTree", str), ("outFile", str)],
    #        [("sortAttr", str, None)],
    #        __doc__)
    main(**args)

