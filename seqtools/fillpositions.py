#!/usr/bin/env python3
# -*- coding: utf-8 -*-


from sys import stdin, stdout
import argparse as ap
from Bio import SeqIO
import logging
logger = logging.getLogger(__name__)


"""Expects HmmCleaner output, such as:

ENSMICG00000004206
        1918-5851
        18565-18615
ENSANAG00000022814
ENSCCAG00000037294
        7681-7775
ENSPCOG00000028606
ENSSBOG00000020815
        7673-7756
ENSMLEG00000030251
        7675-7749
ENSMMUG00000016936
        1414-5851
        18565-18615
ENSG00000105926
ENSCSAG00000014108
ENSPPYG00000017722
"""


def parse_seqranges(positionsfile, to_codons=False):  #should be FROM codons/FROM aa.
    ranges = {}
    multiply = 3 if to_codons else 1
    with open(positionsfile) as f:
        current_seq = None
        for line in f:
            if not line.startswith('\t'):
                current_seq = line.rstrip()
                ranges[current_seq] = []
            else:
                seqrange = tuple(int(x) for x in line.strip().split('-'))
                # Convert to 0-based coordinates with end excluded.
                seqrange = ((seqrange[0]-1)*multiply,
                            seqrange[1]*multiply)
                ranges[current_seq].append(seqrange)
    return ranges


def main(positionsfile, alignmentfile=stdin, to_codons=False, fillchar='-',
         format='fasta', outputfile=stdout):

    ranges = parse_seqranges(positionsfile, to_codons)

    if alignmentfile == '-':
        alignmentfile = stdin

    recordlist = []
    for record in SeqIO.parse(alignmentfile, format):
        slen = len(record.seq)
        record.seq = record.seq.tomutable()  # Edit in place.
        seqranges = ranges[record.name]
        logger.debug('Filling %d ranges in %s', len(seqranges), record.name)
        for start, end in seqranges:
            logger.debug('range %d-%d', start, end)
            length = end - start
            logger.debug('Before: %s', record.seq[start:end])
            record.seq[start:end] = fillchar * length
            logger.debug('After:  %s', record.seq[start:end])
        assert slen == len(record.seq)
        recordlist.append(record)

    SeqIO.write(recordlist, outputfile, format)


if __name__ == '__main__':
    logging.basicConfig(format=logging.BASIC_FORMAT)
    parser = ap.ArgumentParser(description=__doc__)
    parser.add_argument('positionsfile')
    parser.add_argument('alignmentfile', nargs='?', default=stdin)
    parser.add_argument('-f', '--format', default='fasta',
                        help='Input and output sequence format [%(default)s]')
    parser.add_argument('-o', '--outputfile', default=stdout)
    parser.add_argument('-c', '--to-codons', action='store_true',
                        help='Multiply ranges boundaries by 3 to fill codon '\
                        'alignments using amino-acid positions')
    parser.add_argument('-F', '--fillchar', default='-',
                        help='Filling character [%(default)s]')
    parser.add_argument('-v', '--verbose', action='count', default=0)
    args = parser.parse_args()

    if args.verbose > 1:
        logger.setLevel(logging.DEBUG)
    elif args.verbose > 0:
        logger.setLevel(logging.INFO)
    delattr(args, 'verbose')

    main(**vars(args))
