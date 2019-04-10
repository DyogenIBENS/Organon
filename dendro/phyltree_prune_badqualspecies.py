#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse as ap
from LibsDyogen import myPhylTree


def main(phyltreefile):
    phyltree = myPhylTree.PhylogeneticTree(phyltreefile)
    phyltree.pruneSpecies(phyltree.lstEsp2X | phyltree.lstEsp6X, inplace=True)
    phyltree.printNewick(commonnames=True, symbols=True)

if __name__ == '__main__':
    parser = ap.ArgumentParser(description=__doc__)
    parser.add_argument('phyltreefile')
    
    args = parser.parse_args()
    main(**vars(args))
