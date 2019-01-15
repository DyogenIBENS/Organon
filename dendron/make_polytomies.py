#!/usr/bin/env python3
# -*- coding: utf-8 -*-


from sys import stdin
import argparse as ap
import ete3


def main(conversionfile, treefile=None):
    with open(conversionfile) as stream:
        conversion = {}
        for line in stream:
            if not line.startswith('#'):
                field1, field2 = line.rstrip().split('\t')
                if field1 and field2:
                    conversion[field1] = field2

    if treefile is None:
        treefile = stdin.read()
    tree = ete3.Tree(treefile, format=1)
    
    for node in tree.traverse('postorder'):
        newname = conversion.get(node.name, node.name)
        if '+' in newname:
            # This clade is polyphyletic in the reference tree.
            # Transform the parent into a multifurcation.
            node.delete(prevent_nondicotomic=False, preserve_branch_length=True)
            #TODO: check that newname.split('+') matches the node.children.

    print(tree.write(format=1, format_root_node=True))

if __name__ == '__main__':
    parser = ap.ArgumentParser(description=__doc__)
    parser.add_argument('conversionfile')
    parser.add_argument('treefile', nargs='?')
    
    args = parser.parse_args()
    main(**vars(args))

