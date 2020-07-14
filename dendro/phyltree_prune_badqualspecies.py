#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from sys import stderr
import argparse as ap
from LibsDyogen import myPhylTree


def main(phyltreefile, specieslistfile=None, reverse=False):
    phyltree = myPhylTree.PhylogeneticTree(phyltreefile)
    if specieslistfile is None:
        badspeciesset = phyltree.lstEsp2X | phyltree.lstEsp6X
    else:
        with open(specieslistfile) as f:
            badspeciesset = set(line.rstrip() for line in f)
        if badspeciesset - phyltree.listSpecies:
            print('WARNING: unknown species:',
                  ', '.join(badspeciesset-phyltree.listSpecies),
                  file=stderr)
    if reverse:
        badspeciesset = phyltree.listSpecies - badspeciesset
    phyltree.pruneSpecies(badspeciesset, inplace=True)
    phyltree.printNewick(commonnames=True, symbols=True)


if __name__ == '__main__':
    parser = ap.ArgumentParser(description=__doc__)
    parser.add_argument('phyltreefile')
    parser.add_argument('specieslistfile', nargs='?',
                        help='Species to remove. If None, remove 2X and 6X species.')
    parser.add_argument('-v', '--reverse', action='store_true',
                        help='Only keep species given in the list.')
    
    args = parser.parse_args()
    main(**vars(args))
