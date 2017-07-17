#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Converts a treebest-produced newick tree (from `treebest best`) into a tree
with internal node names (taxon + genename). Also strip the species name from 
the leaf names."""

import ete3
import argparse


def get_ascii_suffix(integer):
    assert 0 <= integer <= 2*26
    if integer <= 26:
        return chr(97 + integer)
    else:
        return chr(65 + integer)


def treebest2genomicus(treebestnwk, genetreename, outfile):
    treebesttree = ete3.Tree(treebestnwk)

    treebesttree.name = treebesttree.S + genetreename
    treebesttree.add_feature('G', genetreename)

    for node in treebesttree.traverse():
        if not node.is_leaf():
            for i, child in enumerate(node.children):
                child.add_feature('G', node.G)
                if node.D == 'Y':
                    child.G += '.' + get_ascii_suffix(i)
                if child.is_leaf():
                    child.name = child.name.split('_')[0]
                else:
                    child.name = child.S + child.G
        
    treebesttree.write(outfile=outfile, format=1, format_root_node=True)



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('treebestnwk')
    parser.add_argument('genetreename')
    parser.add_argument('outfile')
    
    args = parser.parse_args()
    treebest2genomicus(**vars(args))
