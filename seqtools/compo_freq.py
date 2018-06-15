#!/usr/bin/env python3
# -*- coding: utf-8 -*-


"""Get compositional frequencies of an alignment: GC content, gaps, 'N's.

For the gap count to be meaningful, you should remove columns containing only gaps
(use `ungap.py`)
"""


from sys import stdin
import argparse
from Bio import AlignIO


def compo_freq(alignment, format='fasta'):
    nucl_count = 0
    GC_count   = 0
    gap_count  = 0
    N_count    = 0
    CpG_count  = 0

    length = alignment.get_alignment_length()

    for record in alignment:
        gaps = record.seq.count('-')
        N = record.seq.count('N')
        A = record.seq.count('A')
        T = record.seq.count('T')
        G = record.seq.count('G')
        C = record.seq.count('C')
        CpG = record.seq.count('CG')

        assert gaps + N + A + T + G + C == length

        nucl_count += A + T + G + C
        GC_count   += G + C
        gap_count  += gaps
        N_count    += N
        CpG_count  += CpG
    
    return (length,
            float(GC_count)      / nucl_count,
            float(N_count)       / (nucl_count+N_count),
            float(gaps)          / length,
            float(CpG_count) * 2 / nucl_count)

def main(alignment_file, format='fasta'):
    alignment = AlignIO.read(alignment_file, format=format) # Alphabet=
    freqs = compo_freq(alignment, format)
    print('\t'.join(['%7s']*len(freqs)) % ('len', 'GC', 'N', 'gaps', 'CpG'))
    print('\t'.join(['%7d'] + ['%7.4f']*(len(freqs)-1)) % freqs)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('alignment_file', nargs='?',
                        type=argparse.FileType('r'), default=stdin)
    parser.add_argument('-f', '--format', default='fasta')
    
    args = parser.parse_args()
    main(**vars(args))
