#!/usr/bin/env python3
# -*- coding: utf-8 -*-


"""Merge selected disjoint sequences in an alignment"""


from sys import stdin, stdout, stderr
import argparse

import numpy as np
from functools import reduce

from Bio import Align, AlignIO
from Bio.Seq import Seq


def get_chunks(a):
    """Compute the boundaries of chunks (intervals) of True values in an array"""
    a = a.astype(bool).astype(int)
    # Add outer boundaries:
    a = np.concatenate(([0], a, [0]))
    gradient = a[1:] - a[:-1] # Now 1 means start, -1 means end (excluded).
    starts, = np.where(gradient == 1)
    ends,   = np.where(gradient == -1)
    assert len(starts) == len(ends)
    #chunks = np.array(list(zip(starts, ends)))
    return np.stack((starts, ends), axis=-1)
    
def intersection_type(target_chunks, test_chunk):
    """- target_chunks: 2D array: chunks in rows (column 0: start, column 1: end)"""
    start_in = (test_chunk[0] >= target_chunks[:,0]) & (test_chunk[0] <  target_chunks[:,1])
    end_in   = (test_chunk[1] >  target_chunks[:,0]) & (test_chunk[1] <= target_chunks[:,1])
    return start_in, end_in

def which_chunk_intersect(target_chunks, test_chunk):
    """Return indices of the target chunks that are intersected by the test_chunk."""
    start_in, end_in = intersection_type(target_chunks, test_chunk)
    return (start_in | end_in).nonzero()[0]

def chunk_contained_in(target_chunks, test_chunk):
    start_in, end_in = intersection_type(target_chunks, test_chunk)
    return (start_in & end_in)


def al_merge_splits(align, split_seqs):
    split_seqs.sort()
    seq_index = {record.id: i for i, record in enumerate(align)}
    merged_records =     [align[seq_index[s]] for s in split_seqs]
    merged_seqs =        [np.array(rec.seq) for rec in merged_records]
    
    # Check conflicting chunks
    merged_seqs_nogaps = [s != '-' for s in merged_seqs]
    merged_nogap_intersection = reduce(np.bitwise_and, merged_seqs_nogaps)

    if merged_nogap_intersection.any():
        print("WARNING: Sequences coordinates intersect (%d residues, identical: %s)" \
                         % (merged_nogap_intersection.sum(),
                            np.array_equal(merged_seqs[0][merged_nogap_intersection],
                                  merged_seqs[1][merged_nogap_intersection])),
              file=stderr)
        print("conflicting slices:",
              " ".join("%d-%d" % tuple(ch) for ch in get_chunks(merged_nogap_intersection)),
              file=stderr)

    # Start merging
    new_rec = merged_records[0]
    new_seq = merged_seqs[0]

    for i, merged_seq in enumerate(merged_seqs[1:], start=1):
        merged_positions = (merged_seq != '-')
        newseq_notgap = (new_seq != '-')
        
        if any(newseq_notgap[merged_positions]):
            # Conflicting chunk
            newseq_chunks    = get_chunks(newseq_notgap)
            mergedseq_chunks = get_chunks(merged_positions)
            conflict_chunks  = get_chunks(newseq_notgap & merged_positions)

            for conflict_chunk in conflict_chunks:
                newseq_conflict_src = newseq_chunks[which_chunk_intersect(
                                                        newseq_chunks,
                                                        conflict_chunk),
                                                    ]
                mergedseq_conflict_src = mergedseq_chunks[which_chunk_intersect(
                                                            mergedseq_chunks,
                                                            conflict_chunk),
                                                    ]
                
                assert len(newseq_conflict_src) == 1 and len(mergedseq_conflict_src) == 1

                if np.array_equal(conflict_chunk, newseq_conflict_src[0]):
                    if np.array_equal(conflict_chunk, mergedseq_conflict_src):
                        raise ValueError("Can't decide")
                    else:
                        # The smallest fragment comes from 'newseq', so overwrite it.
                        print('Discarding smaller conflicting fragment %s from %s'\
                                %(conflict_chunk, split_seqs[0]), file=stderr)

                elif np.array_equal(conflict_chunk, mergedseq_conflict_src[0]):
                    # The merged seq has the smallest fragment, discard it.
                    merged_positions[conflict_chunk[0]:conflict_chunk[1]] = False
                    print('Discarding smaller conflicting fragment %s from %s'\
                            %(conflict_chunk, split_seqs[i]), file=stderr)
                else:
                    assert chunk_contained_in(newseq_conflict_src, conflict_chunk).all() and \
                           chunk_contained_in(mergedseq_conflict_src, conflict_chunk).all(), \
                           "Coordinate error: %s from %s or %s" % (conflict_chunk,
                                    newseq_conflict_src,
                                    mergedseq_conflict_src)

                    src_chunk_lengths = [(ch[1] - ch[0]) for ch in \
                                         (newseq_conflict_src[0],
                                          mergedseq_conflict_src[0])]
                    if src_chunk_lengths[0] > src_chunk_lengths[1]:
                        # newseq has the longest source. keep.
                        print('Discarding smaller conflicting fragment %s from %s'\
                                %(conflict_chunk, split_seqs[i]), file=stderr)
                    else:
                        print('Discarding smaller conflicting fragment %s from %s'\
                                %(conflict_chunk, split_seqs[0]), file=stderr)

        new_seq[merged_positions] = merged_seq[merged_positions]

    new_rec.seq = Seq(''.join(new_seq), alphabet=new_rec.seq.alphabet)
    new_rec.description = 'Split gene merged from: ' + ','.join(split_seqs)
    
    new_al = Align.MultipleSeqAlignment([rec for rec in align \
                                         if rec.id not in split_seqs[1:]])
    return new_al


def main_fromlist(input_al, split_seqs_list, output_al, format='fasta'):
    align = AlignIO.read(input_al, format=format)
    for split_seqs in split_seqs_list:
        align = al_merge_splits(align, split_seqs)

    AlignIO.write(align, output_al, format)


def main(input_al, split_seqs_file, output_al, format='fasta'):
    #with open(split_seqs_file) as f:
    split_seqs_list = [line.split() for line in split_seqs_file]

    main_fromlist(input_al, split_seqs_list, output_al, format)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input_al')
    parser.add_argument('split_seqs_file', nargs='?', type=argparse.FileType('r'),
                        default=stdin)
    parser.add_argument('-o', '--output-al', type=argparse.FileType('w'),
                        default=stdout)
    
    args = parser.parse_args()
    main(**vars(args))


