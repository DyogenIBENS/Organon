#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
"""Append the species name to a gene name in a sequence alignment.
Can convert ensembl gene names to species, as well as assembly names
(ex: loxAfr3)."""

import os.path as op
import argparse
from Bio import SeqIO
#from codeml.select_leaves_from_specieslist import convert_gene2species 
from genomicustools.identify import convert_prot2species


def load_conversion(filename="~glouvel/ws2/UCSC_genome_releases_full.tsv",
                    fromcol=2, tocol=0):
    conversion = {}
    with open(op.expanduser(filename)) as stream:
        header = next(stream)
        for line in stream:
            fields = line.rstrip().split('\t')
            conversion[fields[fromcol]] = fields[tocol]
    return conversion


def iter_specify(inputfile, format="fasta", dots=True):
    a2sp = load_conversion()
    for record in SeqIO.parse(inputfile, format):
        try:
            sp = convert_prot2species(record.id, ensembl_version=87)
        except KeyError:
            record.id = record.id.split('_')[0]
            sp = a2sp[record.id]
        if dots:
            sp = sp.replace(' ', '.')
        record.id += "_" + sp
        record.description = ''
        #print(record.id)
        yield record


def specify(inputfile, outfile, format="fasta", dots=True):
    SeqIO.write(iter_specify(inputfile), outfile, format)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('inputfile')
    parser.add_argument('outfile') 
    parser.add_argument('-f', '--format', default="fasta", 
                        help='[%(default)s]')
    args = parser.parse_args()
    specify(**vars(args))

