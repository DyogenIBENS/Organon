#!/usr/bin/env python3
# -*- coding: utf-8 -*-


"""Get compositional frequencies of an alignment: GC content, gaps, 'N's.

For the gap count to be meaningful, you should remove columns containing only gaps
(use `ungap.py`)
"""


from sys import stdin
import numpy as np
import argparse
from Bio import AlignIO


def compo_freq(alignment, format='fasta'):
    # Statistics for each sequence
    seq_nucl = np.zeros(len(alignment))
    seq_GC   = np.zeros(len(alignment))
    seq_gap  = np.zeros(len(alignment))
    seq_N    = np.zeros(len(alignment))
    seq_CpG  = np.zeros(len(alignment))

    length = alignment.get_alignment_length()

    for i, record in enumerate(alignment):
        gaps = record.seq.count('-')
        N = record.seq.count('N')
        A = record.seq.count('A')
        T = record.seq.count('T')
        G = record.seq.count('G')
        C = record.seq.count('C')
        CpG = record.seq.count('CG')

        assert gaps + N + A + T + G + C == length

        seq_nucl[i] = A + T + G + C
        seq_GC[i]   = G + C
        seq_gap[i]  = gaps
        seq_N[i]    = N
        seq_CpG[i]  = CpG
    
    nucl_count = seq_nucl.sum()
    GC_count   = seq_GC.sum()
    gap_count  = seq_gap.sum()
    N_count    = seq_N.sum()
    CpG_count  = seq_CpG.sum()

    global_stats = (length,
                    float(GC_count)      / nucl_count,
                    float(N_count)       / (nucl_count+N_count),
                    float(gaps)          / length,
                    float(CpG_count) * 2 / nucl_count)

    seq_lengths = seq_nucl + seq_N
    seq_GC_freq = seq_GC.astype(float) / seq_nucl
    seq_N_freq  = seq_N.astype(float) / seq_lengths
    seq_gap_freq = seq_gap.astype(float) / length
    seq_CpG_freq = seq_CpG * 2. / seq_nucl
    
    seq_stats = (seq_lengths,
                 seq_GC_freq,
                 seq_N_freq,
                 seq_gap_freq,
                 seq_CpG_freq)

    seq_means   = [s.mean()     for s in seq_stats]
    seq_medians = [np.median(s) for s in seq_stats]
    seq_stds    = [s.std()      for s in seq_stats]

    return global_stats, seq_means, seq_medians, seq_stds

def main(alignment_file, format='fasta'):
    alignment = AlignIO.read(alignment_file, format=format) # Alphabet=
    global_stats, *seq_stats = compo_freq(alignment, format)

    n_stats = len(global_stats)
    print('\t'.join(['%8s']*n_stats) % ('len', 'GC', 'N', 'gaps', 'CpG'))
    print('\t'.join(['%8d'] + ['%8.4g']*(n_stats-1)) % global_stats)
    seq_stat_template = '\t'.join(['%8.4g']*n_stats)

    for desc, seq_stat in zip(['mean', 'med', 'std'], seq_stats):
        print(desc)
        print(seq_stat_template % tuple(seq_stat))
    
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('alignment_file', nargs='?',
                        type=argparse.FileType('r'), default=stdin)
    parser.add_argument('-f', '--format', default='fasta')
    
    args = parser.parse_args()
    main(**vars(args))
