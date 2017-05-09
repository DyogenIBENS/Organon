#!/usr/bin/env python3

"""Join sequences from multifasta file, and reverse complement when on strand -.
"""

"""
USAGE:
    ./join_revcomp_fasta.py <input_fasta> <info_file> <outputfile>
"""

import sys
import re
import bz2
import argparse

HEADER_FMT = r'^(\w+):(\d+)-(\d+)$'

def myopen(filename, *args, **kwargs):
    if filename.endswith('.bz2'):
        f = bz2.open(filename, *args, **kwargs)

        orig_readline = f.readline
        def readdecodeline(*args, **kwargs):
            return orig_readline().decode()
        f.readline = readdecodeline
        return f
    else:
        return open(filename, *args, **kwargs)


COMPL = {'A': 'T',
         'T': 'A',
         'G': 'C',
         'C': 'G',
         'N': 'N',
         'a': 't',
         't': 'a',
         'g': 'c',
         'c': 'g',
         'n': 'n'}

def reverse_complement(seq):
    return ''.join(COMPL[nucl] for nucl in reversed(seq))

def load_info(info_file):
    strand_info = {}
    with open(info_file) as stream:
        for line in stream:
            seq, start, end, strand = line.split()
            info_seq = strand_info.setdefault(seq, {})
            info_seq['%s-%s' % (start, end)] = strand
    return strand_info


def iter_multifasta(input_fasta):
    """iterate over a multifasta file. Yields (seqname, seqrange, joined sequence)"""
    with myopen(input_fasta) as fin:
        line = fin.readline()
        #i = 0

        while line:
            header = line.rstrip().lstrip('>')
            #print("seq %02d: %s:%s" % (i, seqname, seqrange))
            #i += 1

            nextline = fin.readline()
            bloc = ''
            while nextline and not nextline.startswith('>'):
                bloc += nextline.rstrip()
                nextline = fin.readline()

            yield header, bloc
            
            line = nextline


def write_fastaseq(out, header, sequence, width=None):
    out_txt = '>' + header + '\n'
    if width:
        while sequence:
            out_txt += sequence[:width] + '\n'
            sequence = sequence[width:]
    else:
        out_txt += sequence + '\n'
    out.write(out_txt)


def parse_header(header, fmt=HEADER_FMT, types=(str, int, int)):
    """parse the header line of a fasta according to this format:
    'seqname:start-end'
    """
    #seqname, seqrange = header.split(':')
    #start, end = [int(coord) for coord in seqrange.split('-')]
    #return seqname, start, end
    return [typ(val) for val, typ in zip(re.search(fmt, header).groups(), types)]


def main(input_fasta, outputfile, info_file=None, header_fmt=HEADER_FMT,
         fill=79, rev_coords=False, revcomp=False):
    width = 79
    default_fill_gap = 'N' * fill
    output_seq = ''
    output_seqnames = []

    strand_info = load_info(info_file) if info_file else {}
    prev_seqname = None
    prev_start, prev_end = None, None
    prev_strand = None

    for header, seq in iter_multifasta(input_fasta):
        seqname, start, end = parse_header(header, header_fmt)
        output_seqnames.append('%s:%d-%d' % (seqname, start, end))
        strand = strand_info[seqname]['%d-%d' % (start, end)] if info_file else '+'
        fill_gap = default_fill_gap if prev_end else ''
        if strand == '-':
            #print('reverse-complementing')
            if revcomp: seq = reverse_complement(seq)
            output_seqnames[-1] += '[revcomp]'
        if seqname == prev_seqname and strand == prev_strand:
            # works when coords are given for the + strand
            if (start - prev_end - 1) >= 0:
                fill_gap = 'N' * (start - prev_end - 1)

        output_seq += fill_gap + seq

        prev_seqname = seqname
        prev_start, prev_end = start, end
        prev_strand = strand


    if outputfile == '-':
        out = sys.stdout
    else:
        out = open(outputfile, 'w')

    write_fastaseq(out, '; '.join(output_seqnames), output_seq, width)

    if outputfile != '-':
        out.close()

if __name__=='__main__':
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument('input_fasta')
    parser.add_argument('outputfile')
    parser.add_argument('-i', '--info_file', help='tabulated format with 4 columns: '\
                        'contig name, start, end, strand.')
    parser.add_argument('-f', '--fill', type=int, default=79,
                        help="number of 'N' to put between joined contigs [%(default)s].")
    parser.add_argument('--rev-coords', action='store_true', help='Not implemented')
    parser.add_argument('--revcomp', action='store_true',
                        help='reverse complement sequences on minus strand.')
    parser.add_argument('-H', '--header-fmt', default=HEADER_FMT, 
                        help='format of header [%(default)r]')
    args = parser.parse_args()

    main(**vars(args))
