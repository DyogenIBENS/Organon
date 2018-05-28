#!/usr/bin/env python3
# -*- coding: utf-8 -*-


"""Merge selected disjoint sequences in an alignment"""


from sys import stdin, stdout
import argparse

import numpy as np
from functools import reduce

from Bio import Align, AlignIO
from Bio.Seq import Seq


def al_merge_splits(align, split_seqs):
    split_seqs.sort()
    seq_index = {record.id: i for i, record in enumerate(align)}
    merged_records =     [align[seq_index[s]] for s in split_seqs]
    merged_seqs =        [np.array(rec.seq) for rec in merged_records]
    merged_seqs_nogaps = [s != '-' for s in merged_seqs]
    
    merged_nogap_intersection = reduce(np.bitwise_and, merged_seqs_nogaps)

    assert not merged_nogap_intersection.any()

    new_rec = merged_records[0]
    new_seq = merged_seqs[0]

    for merged_seq in merged_seqs[1:]:
        merged_positions = (merged_seq != '-')
        new_seq[merged_positions] = merged_seq[merged_positions]

    new_rec.seq = Seq(''.join(new_seq), alphabet=new_rec.seq.alphabet)
    
    new_al = Align.MultipleSeqAlignment([rec for rec in align \
                                         if rec.id not in split_seqs[1:]])
    return new_al


def main(input_al, split_seqs_filehandle, output_al, format='fasta'):
    #with open(split_seqs_file) as f:
    split_seqs_list = [line.split() for line in split_seqs_filehandle]

    align = AlignIO.read(infile, format=format)
    for split_seqs in split_seqs_list:
        align = al_merge_splits(align, split_seqs)

    AlignIO.write(align, output_al, format)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input_al')
    parser.add_argument('split_seqs_file', nargs='?', type=argparse.FileType('r'),
                        default=stdin)
    parser.add_argument('-o', '--output-al', type=argparse.FileType('w'),
                        default=stdout)
    
    args = parser.parse_args()
    main(**vars(args))










