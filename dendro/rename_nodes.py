#!/usr/bin/env python3
# -*- coding: utf-8 -*-


from sys import stdin
import argparse as ap
import ete3
from dendro.converters import PhylTree_to_ete3

def load_conversion(conversionfile):
    with open(conversionfile) as stream:
        conversion = {}
        for line in stream:
            if not line.startswith('#'):
                field1, field2, *extra = line.rstrip('\r\n').split('\t')
                if field1 and field2:
                    conversion[field1] = field2
    return conversion


def rename(tree, conversion):
    for node in tree.traverse():
        try:
            node.name = conversion[node.name]
        except KeyError:
            pass


def main(conversionfile, treefile=None, parser='ete3'):
    conversion = load_conversion(conversionfile)

    if parser == 'ete3':
        if treefile is None:
            treefile = stdin.read()
        tree = ete3.Tree(treefile, format=1)
    elif parser == 'myPhylTree':
        from LibsDyogen import myPhylTree
        if treefile is None:
            treefile = stdin
        tree = PhylTree_to_ete3(myPhylTree.PhylogeneticTree(treefile), nosinglechild=False)
    
    rename(tree, conversion)

    print(tree.write(format=1, format_root_node=True))


if __name__ == '__main__':
    parser = ap.ArgumentParser(description=__doc__)
    parser.add_argument('conversionfile')
    parser.add_argument('treefile', nargs='?')
    parser.add_argument('-p', '--parser', default='ete3',
                        choices=['ete3', 'myPhylTree', 'myProteinTree'],
                        help='How to load the tree. `myPhylTree` allows to ' \
                             'automatically set unique node names [%(default)s].')
    
    args = parser.parse_args()
    main(**vars(args))

