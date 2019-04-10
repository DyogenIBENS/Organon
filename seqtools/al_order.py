#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Reorder sequences in alignment to match the leaves of a given tree, or an
ordered list of record names."""

from sys import stdin, stdout
import argparse
from Bio.Align import MultipleSeqAlignment
from Bio.AlignIO import read as alignread, write as alignwrite
from Bio.Phylo import read as phyloread


def main(infile, outfile, treefile=None, recordsfile=None, format='fasta'):
    
    align = alignread(infile, format=format)

    recnames = [rec.name for rec in align]

    if recordsfile:
        with open(recordsfile) as recf:
            records = [recnames.index(line.rstrip()) for line in recf]
    elif treefile:
        tree = phyloread(treefile, 'newick')
        records = [recnames.index(leaf.name) for leaf in tree.get_terminals()]

    align = MultipleSeqAlignment([align[r] for r in records])

    alignwrite(align, outfile, format)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('infile', nargs='?', type=argparse.FileType('r'),
                        default=stdin)
    parser.add_argument('outfile', nargs='?', type=argparse.FileType('w'),
                        default=stdout)
    order_group = parser.add_mutually_exclusive_group(required=True)
    order_group.add_argument('-t', '--treefile')
    order_group.add_argument('-r', '--recordsfile')
    parser.add_argument('-f', '--format', default='fasta')
    
    args = parser.parse_args()
    main(**vars(args))
