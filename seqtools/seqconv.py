#!/usr/bin/env python3
# -*- coding: utf-8 -*-


from sys import stdin, stdout
import argparse
from Bio import AlignIO


def write_phylip_sequential_relaxed(al, outfile):
    """Bio implements:
    - phylip (strict, interleaved),
    - phylip-sequential (strict, sequential),
    - phylip-relaxed (relaxed, interleaved),
      
    but not a relaxed sequential (as needed for PhyloBayes).
    
    Write sequences on a single line.
    """
    
    if outfile is not stdout:
        outfile = open(outfile)
    try:
        outfile.write('%d\t%d\n' % (len(al), al.get_alignment_length()))
        for record in al:
            outfile.write('%s  %s\n' % (record.name, record.seq))
    finally:
        if outfile is not stdout:
            outfile.close()


def convert_to_phylip_sequential_relaxed(infile, outfile, fro="fasta"):
    write_phylip_sequential_relaxed(AlignIO.read(infile, fro), outfile)


def main(infile, outfile, fro="fasta", to="phylip-relaxed"):
    if to == 'phylip-sequential-relaxed':
        convert_to_phylip_sequential_relaxed(infile, outfile, fro)
    else:
        AlignIO.convert(infile, fro, outfile, to)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('infile', nargs='?', #type=argparse.FileType('r'),
                        default=stdin)
    parser.add_argument('outfile', nargs='?', #type=argparse.FileType('w'),
                        default=stdout)
    # TODO: inplace edition (hum, bad idea)
    parser.add_argument('-f', '--from', dest='fro', default='fasta')
    parser.add_argument('-t', '--to', default='phylip-relaxed')
    
    args = parser.parse_args()
    main(**vars(args))

