#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
"""Append the species name to a gene name in a sequence alignment.
Can convert ensembl gene names to species, as well as assembly names
(ex: loxAfr3)."""

import argparse
from Bio import SeqIO
from genomicustools.identify import ultimate_seq2sp


ENSEMBL_VERSION = 87


def iter_specify(inputfile, format="fasta", dots=True, ensembl_version=ENSEMBL_VERSION):
    for record in SeqIO.parse(inputfile, format):
        sp = ultimate_seq2sp(record.id, ensembl_version)
        if dots:
            sp = sp.replace(' ', '.')
        record.id += "_" + sp
        record.description = ''
        #print(record.id)
        yield record


def specify(inputfile, outfile, format="fasta", dots=True,
            ensembl_version=ENSEMBL_VERSION):
    SeqIO.write(iter_specify(inputfile, format, dots, ensembl_version),
                outfile, format)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('inputfile')
    parser.add_argument('outfile') 
    parser.add_argument('-f', '--format', default="fasta",
                        help='[%(default)s]')
    parser.add_argument('-e', '--ensembl-version', type=int,
                        default=ENSEMBL_VERSION, help='[%(default)s]')

    args = parser.parse_args()
    specify(**vars(args))
