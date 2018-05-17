#!/usr/bin/env python3
# -*- coding: utf-8 -*-


"""Remove **columns only containing gaps**. Runs in codon mode by default"""


from sys import stderr, stdin, stdout
import argparse as ap
from itertools import islice
from Bio import AlignIO


def main(infile, format='fasta', nucl=False, outfile=None):
    align = AlignIO.read(infile, format=format)
    N = align.get_alignment_length()
    gaps = []
    curr_gap_start = None
    #curr_gap_end = None
    for i in islice(range(N), 0, N, 3):
        if all(rec.seq == '---' for rec in align[:, i:(i+3)]):
            if curr_gap_start is None:
                curr_gap_start = i
        else:
            if curr_gap_start is not None:
                gaps.append((curr_gap_start, i))
                curr_gap_start = None
    
    if curr_gap_start is not None:
        gaps.append((curr_gap_start, N))

    tot_gap_len = sum((x2 - x1) for x1,x2 in gaps)
    print("Gap length = %d/%d (Keep %d positions)" % \
            (tot_gap_len, N, N - tot_gap_len),
          file=stderr)

    try:
        gap_start, gap_end = gaps.pop(0)
    except IndexError:
        print('No need to remove gaps', file=stderr)
        return

    ungapped_al = align[:, 0:gap_start]

    while gaps:
        next_gap_start, next_gap_end = gaps.pop(0)
        ungapped_al += align[:, gap_end:next_gap_start]
        gap_end = next_gap_end

        #print(gap_end, next_gap_start, next_gap_end, ':', len(ungapped_al), ungapped_al.get_alignment_length(), file=stderr)

    ungapped_al += align[:, next_gap_end:N]

    #print(len(ungapped_al), ungapped_al.get_alignment_length(), file=stderr)

    AlignIO.write((ungapped_al,), outfile, format=format)


if __name__ == '__main__':
    parser = ap.ArgumentParser(description=__doc__)
    parser.add_argument('infile', nargs='?', type=ap.FileType('r'), default=stdin)
    parser.add_argument('-f', '--format', default='fasta',
                        help='Input sequence format [%(default)s]')
    parser.add_argument('-n', '--nucl', action='store_true',
                        help='Read as nucleotide alignment (default: codon)')
    parser.add_argument('-o', '--outfile', type=ap.FileType('w'), default=stdout)
    
    
    args = parser.parse_args()
    main(**vars(args))
