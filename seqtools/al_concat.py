#!/usr/bin/env python3
# -*- coding: utf-8 -*-


from sys import stdin, stdout
from io import StringIO
from Bio.Alphabet import generic_dna
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from Bio import AlignIO
import argparse as ap


def parse_many_alignments(infiles, fmt='fasta'):
    """
    Iterate either over multiple files, or over one single file split by '//'
    """
    current = ''
    if not infiles or (len(infiles)==1 and infiles[0] in ('-', '/dev/stdin')):
        for line in stdin:
            if line.startswith('//'):
                yield AlignIO.read(StringIO(current), fmt)
                current = ''
            else:
                current += line
    else:
        for infile in infiles:
            yield from AlignIO.parse(infile, fmt)


def al_concat(infiles, fmt='fasta', sort=True):
    concat = {}  # seq id: seqrecord
    length = 0
    seq_set = set()
    for align in parse_many_alignments(infiles, fmt):
        align_len = align.get_alignment_length()
        this_align_ids = set()
        for record in align:
            if record.id in this_align_ids:
                raise ValueError('Identical ids in the same alignment (%r)' % record.id)
            this_align_ids.add(record.id)
            oldseq = concat.setdefault(record.id, '-'*length)
            concat[record.id] += str(record.seq)
        for unseen in seq_set - this_align_ids:
            concat[unseen] += '-' * align_len
        seq_set |= this_align_ids

        length += align_len
    msa = MultipleSeqAlignment([SeqRecord(Seq(seq), id=seqname, description='') for seqname, seq
                                in sorted(concat.items(), key=lambda t: t[0])])
    return msa


def main():
    parser = ap.ArgumentParser(description=__doc__)
    parser.add_argument('infiles', nargs='*',
        help='Either several files, or, if stdin, many alignments joined by "//" lines.')
    parser.add_argument('-f', '--fmt', default='fasta', help='input/output format [%(default)s]')
    args = parser.parse_args()
    msa = al_concat(args.infiles, args.fmt)
    AlignIO.write(msa, stdout, args.fmt)


if __name__ == '__main__':
    main()
