#!/usr/bin/env python3

import sys
import re
import bz2
import argparse

from fasta_join_revcomp import iter_multifasta, write_fastaseq


def main(pattern, fastafile, outfile='-', dealign=False, sort_starts=False, 
         width=None):
    out_seqs = []

    if fastafile.endswith('.gz'):
        pattern = bytes(pattern, encoding='utf8')
    regex = re.compile(pattern)
    regex_nucl = re.compile(r'[^-]')
    #fasta = bz2.open(fastafile) if fastafile.endswith('.bz2') else open(fastafile)
    out = sys.stdout if outfile == '-' else open(outfile, 'w')
    for header, sequence in iter_multifasta(fastafile):
        if regex.search(header):
            #out.write('>' + header + '\n' + sequence + '\n')
            start_nucl = regex_nucl.search(sequence).start()
            end_nucl = regex_nucl.search(''.join(reversed(sequence))).start()
            end_nucl = len(sequence) - end_nucl # end nucleotide not included
            header += '[start_nucl=%d;end_nucl=%d]' % (start_nucl, end_nucl)
            if dealign:
                sequence = sequence.replace('-', '')

            out_seqs.append((start_nucl, end_nucl, header, sequence))

    if sort_starts:
        out_seqs.sort()

    for _, _, header, sequence in out_seqs:
        write_fastaseq(out, header, sequence, width)

    if outfile != '-': out.close()


class getPattern(argparse.Action):
    """Action class for the command-line parser: read the pattern from
    a file."""
    def __call__(self, parser, namespace, values, option_string=None):
        #print(namespace)
        if namespace.file:
            delattr(namespace, 'file')
            filename = values
            with open(filename) as f:
                pattern = r'|'.join(line.rstrip('\n') for line in f)
            print("Read %r:\n  %r" % (filename, pattern))
        else:
            pattern = values
            
        setattr(namespace, self.dest, pattern)


#class patternFromCLI(argparse.Action):
#    """Action class for the command-line parser: read the pattern from
#    a file."""
#    def __call__(self, parser, namespace, values, option_string=None):
#        got_pattern
#        setattr(namespace, self.dest, pattern)
#        print(namespace)


if __name__=='__main__':
    parser = argparse.ArgumentParser()
    #pattern_parser = parser.add_mutually_exclusive_group(required=True)
    parser.add_argument('pattern', nargs='?', action=getPattern)
    parser.add_argument('-f', '--file', action='store_true',
                        help='Obtain pattern from file (one per line)')
    parser.add_argument('fastafile')
    parser.add_argument('-o', '--outfile', default='-')
    parser.add_argument('-d', '--dealign', action='store_true',
                        help='Remove dashes from the sequence')
    parser.add_argument('-w', '--width', type=int, help='wrap width')
    parser.add_argument('-s', '--sort-starts', action='store_true',
                        help='sort by index of first non-dash character.')
    args = parser.parse_args()
    print(args)

    main(**vars(args))

