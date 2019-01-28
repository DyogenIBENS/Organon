#!/usr/bin/env python3
# -*- coding: utf-8 -*-


"""Rotate branches according to the descendant leaf names.

Input/output trees are in newick format.
"""


from sys import stdin
import argparse
import ete3

from sorter import leaf_sort


def getchildren_ete3(tree, node):
    #return [child for (child,dist) in tree.data.get(node, [])]
    return node.get_children()


def getnodeattr_ete3(tree, node, attrname="name"):
    try:
        return getattr(node, attrname)
    except AttributeError as err:
        err.args = (err.args[0] + ' %s' % node,) + err.args[1:]
        raise


def assignchildren_ete3(tree, node, children):
    node.children = children


def main(treefile, outfile, sort_attr='name', like=None, format=1):  #quoted_node_names
    get_attribute = lambda t,n: getnodeattr_ete3(t, n, attrname=sort_attr)
    if like:
        reftree = ete3.Tree(like, format=format)
        attr2pos = {get_attribute(reftree, leaf): i for i, leaf in
                    enumerate(reftree.iter_leaves())}
        N = len(attr2pos)  # Number of leaves
        orig_get_attribute = get_attribute
        get_attribute = lambda t,n: attr2pos.get(orig_get_attribute(t,n), N)

    if not treefile:
        treefile = stdin.read().rstrip()

    tree = ete3.Tree(treefile, format=format)
    leaf_sort(tree, tree,
              getchildren_ete3, assignchildren_ete3,
              get_attribute=get_attribute)
    outtxt = tree.write(outfile=outfile, format=1, format_root_node=True)
    if outtxt:
        print(outtxt)


if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("treefile", nargs='?', default=None,
                        help="Default to stdin")
    parser.add_argument("outfile", nargs='?', default=None,
                        help="Default to stdout")
    parser.add_argument("-a", "--sort-attr", default='name',
                        help="Leaf attribute to select:\n"\
                        " [%(default)s]")
    parser.add_argument("-l", "--like",
                        help="Other newick tree to take the order from.")
    parser.add_argument('-f', '--format', type=int, default=1)

    args = vars(parser.parse_args())

    #args = myTools.checkArgs(
    #        [("proteinTree", str), ("outFile", str)],
    #        [("sortAttr", str, None)],
    #        __doc__)
    main(**args)

