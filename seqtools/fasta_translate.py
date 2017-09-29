#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import argparse
from Bio import SeqIO


def iter_translate(inputfile, format="fasta"):
    for record in SeqIO.parse(inputfile, format):
        record.seq = record.seq.translate()
        assert len(record.seq) > 0
        yield record

def main(inputfile, outfile, format="fasta"):
    SeqIO.write(iter_translate(inputfile), outfile, format)


if __name__=='__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('inputfile')
    parser.add_argument('outfile')
    parser.add_argument('-f', '--format', default="fasta", 
                        help='[%(default)s]')
    
    args = parser.parse_args()
    main(**vars(args))


