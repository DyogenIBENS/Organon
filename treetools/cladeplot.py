#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Plot tree and draw color boxes behind monophyletic groups"""

from math import log10
import ete3
import argparse
import pandas as pd


def main(treefile, outfile, cladelistfile=None, datafile=None, prune=False,
         log_dist=False, values=None, quoted_node_names=False):

    tree = ete3.Tree(treefile, format=1, quoted_node_names=quoted_node_names)
    if cladelistfile:
        with open(cladelistfile) as f:
            cladelist = set((line.rstrip().split('\t')[0] for line in f if \
                                not line.startswith('#')))
    else:
        cladelist = set()

    if datafile:
        assert values
        data = pd.read_csv(datafile, sep='\t', index_col=0, header=0)
        assert all(val in data.columns for val in values)
        value_face = {val: ete3.AttrFace(val, formatter=' %.2f ',
                       fgcolor='grey', fsize=8) for val in values}
        max_val = data[values].max().max()
        min_val = data[values].min().min()
        average_val = data[values].mean().mean()

    alt_col = ['#7FD09C', '#4cbdad']

    ns = [ete3.NodeStyle(bgcolor=col) for col in alt_col]
    default_ns = ete3.NodeStyle(size=0)
    labelface = ete3.AttrFace('name')

    #for node in tree.traverse():
    def mylayout(node):
        if log_dist:
            node.add_feature('orig_dist', node.dist)
            node.dist = log10(node.dist) if node.dist else 0

        if node.name in cladelist and not prune:
            node.set_style(ns[0])
            ns.insert(0, ns.pop()) #Cycle through colors 
            if not node.is_leaf():
                ete3.add_face_to_node(labelface, node, column=0, position='branch-right')
        else:
            node.set_style(default_ns)

        if datafile and node.is_leaf():
            node.add_feature('profile', [data[val][node.name] for val in values])
            node.add_feature('deviation', [0.0]*len(values))
            heatface = ete3.ProfileFace(max_val, min_val, average_val,
                                        width=20*len(values), height=20,
                                        style='heatmap')
            ete3.add_face_to_node(heatface, node, column=1, aligned=True)

            for i, val in enumerate(values, start=2): 
                node.add_feature(val, data[val][node.name])
                ete3.add_face_to_node(value_face[val], node, column=i, aligned=True)

    if prune:
        tree.prune(cladelist, preserve_branch_length=True)
        tree.dist = 0

    tree.render(outfile, mylayout)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('treefile')
    parser.add_argument('outfile')
    parser.add_argument('-c', '--cladelistfile')
    parser.add_argument('-d', '--datafile')
    parser.add_argument('-v', '--values', nargs='+')

    parser.add_argument('-p', '--prune', action='store_true',
                        help='keep only listed clades')
    parser.add_argument('-l', '--log-dist', action='store_true')
    parser.add_argument('-q', '--quoted-node-names', action='store_true')
    
    args = parser.parse_args()
    main(**vars(args))
