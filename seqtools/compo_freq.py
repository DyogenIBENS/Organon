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
from collections import Counter
from seqtools.IUPAC import gaps, nucleotides, unknown, ambiguous, stop_codons

import logging
logger = logging.getLogger(__name__)


def weighted_std(values, weights, axis=None):
    average = np.average(values, axis, weights=weights)
    variance = np.average((values-average)**2, axis, weights=weights)
    return np.sqrt(variance)


def get_seq_counts(alignment):
    """Counts for each sequence"""
    seq_nucl = np.zeros((4, len(alignment)))  # counts of A, C, G, T
    seq_gap  = np.zeros(len(alignment))
    seq_N    = np.zeros(len(alignment))
    seq_CpG  = np.zeros(len(alignment))
    seq_stops = np.zeros(len(alignment))

    length = alignment.get_alignment_length()

    for i, record in enumerate(alignment):

        counts = Counter(record.seq)
        CpG = record.seq.count('CG')

        seq_nucl[:, i] = [counts[n] for n in nucleotides]
        seq_gap[i]  = sum(counts[g] for g in gaps)
        seq_N[i]    = counts[unknown]
        for a, possible in ambiguous.items():
            seq_N[i] += counts[a] * (1-len(possible)/4.)
            for p in possible:
                seq_nucl[nucleotides.index(p), i] += 0.25 * counts[a]

        assert seq_nucl[:, i].sum() + seq_gap[i] + seq_N[i] == length
        seq_stops[i] = sum(str(record.seq[k:k+3]) in stop_codons for k in range(0,length,3))

        seq_CpG[i]  = CpG

    return length, seq_nucl, seq_gap, seq_N, seq_CpG, seq_stops


def get_seq_freqs(length, seq_nucl, seq_gap, seq_N, seq_CpG, seq_stops):
    """composition frequencies **per sequence**"""
    seq_nucltot = seq_nucl.sum(axis=0)
    seq_lengths = seq_nucltot + seq_N
    
    seq_nucl_freq = seq_nucl.astype(float) / seq_nucltot
    seq_GC_freq   = seq_nucl_freq[1:3].sum(axis=0)       # Only way to get the std(G+C)
    seq_N_freq    = seq_N.astype(float)    / seq_lengths
    seq_gap_freq  = seq_gap.astype(float)  / length
    seq_CpG_freq  = seq_CpG * 2. / seq_nucltot  # Is the *2 necessary?
    seq_stop_freq = seq_stops * 3. / seq_lengths
    
    return (seq_lengths,) + tuple(seq_nucl_freq) + \
           (seq_GC_freq, seq_N_freq, seq_gap_freq, seq_CpG_freq, seq_stop_freq)


def get_al_compo_summary(length, seq_counts, seq_freqs):
    """Stats over the whole alignment
    
    Return: for each summary type, a list of values for each compositional element.
    """
    #seq_nucl, seq_gap, seq_N, seq_CpG = seq_counts
    seq_nucl, seq_gap, seq_N, seq_CpG, seq_stops = seq_counts
    
    seq_lengths = seq_freqs[0]

    nucl_count = seq_nucl.sum(axis=1)  # total number of each nucleotide (summed over sequences)
    GC_count   = nucl_count[1] + nucl_count[2]
    nucl_tot   = nucl_count.sum()
    gap_count  = seq_gap.sum()
    N_count    = seq_N.sum()
    CpG_count  = seq_CpG.sum()
    stop_count = seq_stops.sum()

    global_stats = (length,) \
                   + tuple(nucl_count.astype(float) / nucl_tot) \
                   + (float(GC_count)      / (nucl_tot),
                   +  float(N_count)       / (nucl_tot+N_count),
                      float(gap_count)     / length,
                      float(CpG_count) * 2 / nucl_tot, float(stop_count)* 3 / nucl_tot)

    seq_means   = [s.mean()     for s in seq_freqs]
    seq_medians = [np.median(s) for s in seq_freqs]
    seq_stds    = [s.std()      for s in seq_freqs]

    # Same as glob in all cases but "len" and "gaps".
    # Both variables don't need to be weighted by the number of non-gaps.
    seq_w_means = [np.average(s, weights=seq_lengths)   for s in seq_freqs]

    seq_w_stds  = [weighted_std(s, weights=seq_lengths) for s in seq_freqs]

    # TODO: return a np.array
    return global_stats, seq_means, seq_medians, seq_stds, seq_w_means, seq_w_stds


def make_al_compo(alignment, byseq=False):
    length, *seq_counts = get_seq_counts(alignment)

    seq_freqs = get_seq_freqs(length, *seq_counts)
    
    if byseq:
        stats = np.array(seq_freqs).T
        stat_names = [r.name for r in alignment]  # (1 seq per row)

    else:
        stats = get_al_compo_summary(length, seq_counts, seq_freqs)
        stat_names = ['glob', 'mean', 'med', 'std', 'w_mean', 'w_std']  # summary stats
    return stat_names, stats


def main(alignment_file, format='fasta', byseq=False):
    alignment = AlignIO.read(alignment_file, format=format)  # Alphabet=
    stat_names, stats = make_al_compo(alignment, byseq)

    n_stats = len(stats[0])
    
    # column names
    print('\t'.join(['%8s']*(n_stats+1)) % (
          '', 'len', 'A', 'C', 'G', 'T', 'GC', 'N', 'gaps', 'CpG', 'stops'))
    
    row_template = '\t'.join(['%8s'] + ['%8g']*n_stats)

    for desc, stat in zip(stat_names, stats):
        print(row_template % tuple([desc] + list(stat)))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('alignment_file', nargs='?',
                        type=argparse.FileType('r'), default=stdin)
    parser.add_argument('-f', '--format', default='fasta', help='[%(default)s]')
    parser.add_argument('-s', '--byseq', action='store_true')
    
    args = parser.parse_args()
    main(**vars(args))
