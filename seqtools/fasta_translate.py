#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import argparse
from Bio import SeqIO


def iter_translate(inputfile, format="fasta"):
    if inputfile == '-':
        inputfile = sys.stdin
    for record in SeqIO.parse(inputfile, format):
        record.seq = record.seq.translate(gap='-')
        assert len(record.seq) > 0
        yield record

def main(inputfile, outfile, format="fasta"):
    SeqIO.write(iter_translate(inputfile), outfile, format)


if __name__=='__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('inputfile', nargs='?', default=sys.stdin)
    parser.add_argument('outfile', nargs='?', default=sys.stdout)
    parser.add_argument('-f', '--format', default="fasta", 
                        help='[%(default)s]')
    
    args = parser.parse_args()
    main(**vars(args))


