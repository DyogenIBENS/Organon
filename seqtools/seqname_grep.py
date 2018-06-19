#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from sys import stdout, stdin
import re
import argparse
from Bio import SeqIO, Seq, AlignIO, Align


def seqrecords_grep(records_iterable, pattern, negate=False):
    reg = re.compile(pattern)
    
    for record in records_iterable:
        match = reg.search(record.name)
        if (match and not negate) or (not match and negate):
            yield record


def algrep(alignment, pattern, negate=False):
    out_records = seqrecords_grep(alignment, pattern, negate)
    return Align.MultipleSeqAlignment(out_records)


def main(pattern, filename, negate=False, format='fasta', outfile=stdout):
    #al = AlignIO.read(alignment_file, format=format)
    sequences = SeqIO.parse(filename, format=format)
    SeqIO.write(seqrecords_grep(sequences, pattern, negate), outfile, format)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('pattern')
    parser.add_argument('filename', nargs='?', type=argparse.FileType('r'),
                        default=stdin)
    parser.add_argument('-o', '--outfile', nargs='?', type=argparse.FileType('w'),
                        default=stdout)
    parser.add_argument('-F', '--format', default='fasta')
    parser.add_argument('-v', '--negate', action='store_true')
    
    args = parser.parse_args()
    main(**vars(args))
