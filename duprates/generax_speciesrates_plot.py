#!/usr/bin/env python3
# -*- coding: utf-8 -*-



import argparse as ap
from math import log10
import pandas as pd
import ete3
from LibsDyogen import myPhylTree
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from dendro.any import ete3 as ete3_methods, myPhylTree as phyltree_methods
from datasci.graphs import plottree


def get_logitems(phyltree, item, offset=10):
    return [(ch, log10(10+phyltree.ages[item[0]]) - log10(10+phyltree.ages[ch])) for ch,dist in
            phyltree.items.get(item[0], [])]

def main():
    default_figsize=(15,10)
    parser = ap.ArgumentParser(description=__doc__)
    parser.add_argument('speciestreefile')
    parser.add_argument('ratefile')
    parser.add_argument('ratecolumn', nargs='?', default='D')
    parser.add_argument('-o', '--outfile', help='Plot output file')
    args = parser.parse_args()

    #tree = ete3.Tree(args.speciestreefile, format=1)
    phyltree = myPhylTree.PhylogeneticTree(args.speciestreefile)
    rates = pd.read_csv(args.ratefile, sep='\t',
                        names=['taxon', 'length', 'D', 'L', 'T'],
                        index_col=0)

    #edge_colors = pd.Series({n: rates.D[n.name] for n in tree.traverse()})
    if args.outfile is None:
        plt.switch_backend('TkAgg')

    lines, anc_coords, _ = plottree(phyltree,
                                    get_logitems, #phyltree_methods.get_items,
                                    phyltree_methods.get_label,
                                    edge_colors=rates[args.ratecolumn], #edge_colors)
                                    edge_cmap='viridis',
                                    label_nodes=True,
                                    label_params={'alpha': 0.7},
                                    age_from_root=True)
    #lines.axes.set_xscale('log') # Fail
    plt.gcf().set_size_inches(default_figsize)

    if args.outfile:
        plt.savefig(args.outfile, bbox_inches=False)
    else:
        plt.show()


if __name__=='__main__':
    main()

