#!/usr/bin/env python3

"""Extract given sequences at given coordinates from an axt alignment,
reverse-complement sequences on the minus strand, and output a fasta file.
"""

import sys
import argparse
import gzip


def myopen(filename, *args, **kwargs):
    if filename.endswith('.gz'):
        return gzip.open(filename, *args, **kwargs)
    else:
        return open(filename, *args, **kwargs)


def get_1based_coords_from_bed(bedfile):
    with open(bedfile) as bed:
        for line in bed:
            fields = line.split()
            fields[1] = int(fields[1]) + 1
            fields[2] = int(fields[2])
            return fields # end at first line.


def iter_axt(axtfile):
    counter = 0
    with myopen(axtfile) as axt:
        for line in axt:
            line = line.decode().rstrip()
            if not line.startswith('#'):
                counter += 1
                fields = line.split()
                if len(fields) != 9:
                    raise RuntimeError('There should be 9 fields in the first'\
                                        ' line of each block.')
                seqs = [next(axt).decode().rstrip(),
                        next(axt).decode().rstrip()]

                if next(axt).rstrip():
                    raise RuntimeError('Each block in the axt should consist '\
                                        'of 3 lines + a blank line.')

                for numeric_field in (0,2,3,5,6,8):
                    fields[numeric_field] = int(fields[numeric_field])

                yield fields, seqs
    print("Finished iterating %d blocks." % counter)


def axt2fasta(axtfile, outfile, outinfo, seqname=None, sequence=1, begin=1, end=None):

    # WARNING: coordinates in .axt are 1-based, and the base at the end included.
    # this function converts to 0-based half-open (python-style)

    #fastaname = '>'
    #fastaseq = ''
    #lastcoords = (None, None) # seqname, end

    no_overlap = 0
    overlap = 0
    with open(outfile, 'w') as out, open(outinfo, 'w') as info:
        for fields, seqs in iter_axt(axtfile):
            _, tname, tstart, tend, qname, qstart, qend, strand, score = fields
            if (seqname and seqname != tname) \
                    or tend < begin \
                    or (end and tstart > end):
                # do not overlap
                no_overlap += 1
                continue
            else:
                #print('Overlap: %d-%d (%d, %d)' % (tstart, tend, tend - tstart,
                #                                   len(seqs[sequence])))
                #left_strip = max(0, begin - tstart)
                #right_strip = max(0, tend - end) if end else 0
                #out.write('>%s:%d-%d\n' % (qname, qstart + left_strip, qend - right_strip))
                #info.write('%s\t%d\t%d\t%s\n' % (qname, qstart + left_strip,
                #                                 qend - right_strip, strand))
                #out.write(seqs[sequence][(left_strip+1):(-right_strip)] + '\n')
                out.write('>%s:%d-%d\n' % (qname, qstart-1, qend))
                info.write('%s\t%d\t%d\t%s\n' % (qname, qstart-1, qend, strand))
                out.write(seqs[sequence].replace('-', '') + '\n')
                overlap += 1
    print("Found %d blocks overlapping, %d not overlapping." % (overlap, no_overlap))


if __name__=='__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('axtfile')
    parser.add_argument('outfile')
    parser.add_argument('outinfo')
    parser.add_argument('-s', '--sequence', type=int, default=1,
                        help='sequence to extract: default 1 (VS reference is 0)')
    parser.add_argument('-S', '--seqname', help='name of the contig/chromosome in the target')
    parser.add_argument('-b', '--begin', type=int, default=1,
                        help='beginning, in the reference coordinates (1-based)')
    parser.add_argument('-e', '--end', type=int,
                        help='end, in the reference coordinates (included)')
    parser.add_argument('-B', '--bed',
                        help='bed file with a single line (0-based coordinates)')


    dictargs = vars(parser.parse_args())
    bed = dictargs.pop('bed')
    
    if bed:
        print('Taking coordinates from bed file (overrides command-line options)',
              file=sys.stderr)
        seqname, begin, end = get_1based_coords_from_bed(bed)
        dictargs['seqname'] = seqname
        dictargs['begin'] = begin
        dictargs['end'] = end

    axt2fasta(**dictargs)
