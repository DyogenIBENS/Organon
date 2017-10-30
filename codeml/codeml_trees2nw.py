#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import print_function

"""Convert the .trees format of phylogenetic trees (from the PAML example data)
into newick."""


import argparse
import re

def main(treesfile, outfile):
    with open(treesfile) as trees:
        line0 = trees.readline()
        line = trees.readline()
        stats = [int(s) for s in line0.split()]
        assert len(stats) == 2 and not line.rstrip()

        treelist = ['']

        while line:
            line = line.rstrip()
            if line:
                treelist[-1] += line
                if line.endswith(';'):
                    treelist.append('')

            line = trees.readline()

    last = treelist.pop()
    assert not last and len(treelist) == stats[1]

    with open(outfile, 'w') as out:
        for tree in treelist:
            out.write(tree + '\n')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('treesfile')
    parser.add_argument('outfile')
    
    args = parser.parse_args()
    main(**vars(args))
