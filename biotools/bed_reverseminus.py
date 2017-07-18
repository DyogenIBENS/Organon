#!/usr/bin/env python3


import sys
import argparse


def load_chromsizes(chromsizes_file):
    chromsizes = {}
    with open(chromsizes_file) as f:
        for line in f:
            fields = line.split()
            chromsizes[fields[0]] = int(fields[1])
    return chromsizes

def iter_bed(bedfile, int_columns=[1,2,4]):
    with open(bedfile) as bed:
        for line in bed:
            fields = line.split()
            for col in int_columns:
                try:
                    fields[col] = int(fields[col])
                except IndexError as err:
                    err.args += ('There is no column %d in the file.' % col)
            yield fields
            
def reverse_coords(bedfile, chromsizes_file, n=0, b=1, e=2, s=5):
    chromsizes = load_chromsizes(chromsizes_file)
    for row in iter_bed(bedfile, int_columns=[b,e]):
        if row[s] == '-':
            chrom_len = chromsizes[row[n]]
            yield [row[n], chrom_len - row[e], chrom_len - row[b], row[s]]
        else:
            yield [row[n], row[b], row[e], row[s]]


def main(bedfile, chromsizes_file, outfile='-', n=0, b=1, e=2, s=5):
    out = sys.stdout if outfile == '-' else open(outfile, 'w')
    
    for processed_row in reverse_coords(bedfile, chromsizes_file, n ,b, e, s):
        out.write("%s\t%d\t%d\t%s\n" % tuple(processed_row))

    if outfile != '-': out.close()

if __name__=='__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('bedfile')
    parser.add_argument('chromsizes_file')
    parser.add_argument('-o', '--outfile', default='-', help='[stdout]')
    parser.add_argument('-n', type=int, default=0,
                        help='column number for the chrom name [%(default)s]')
    parser.add_argument('-b', type=int, default=1,
                        help='column number for the sequence start [%(default)s]')
    parser.add_argument('-e', type=int, default=2,
                        help='column number for the sequence end [%(default)s]')
    parser.add_argument('-s', type=int, default=5,
                        help='column number for the strand [%(default)s]')
    args = parser.parse_args()
    main(**vars(args))

