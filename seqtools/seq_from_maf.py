#!/usr/bin/env python3


import re
import sys
import argparse


def parse_assignments(line):
    """return dictionary corresponding to a line like 'score=36 bla=bluh'"""
    d = {}
    for assignment in line.split():
        key, value = assignment.split('=')
        d[key] = value
    return d


def process_maf(maffilename, select_pattern=''):
    """Iterate over alignment blocks (paragraphs) of a MAF file.
    
    Use `select_pattern` to ignore any sequence names not matching it. Use '|'
    to separate alternative patterns, e.g 'hg19|panTro'
    """

    select_regex = re.compile(select_pattern)
    
    with open(maffilename) as maf:
        alignment = None
        sequences = None
        for lineno, line in enumerate(maf):
            line = line.rstrip()
            if line.startswith('#'):
                continue
            elif line.startswith('a '):
                alignment = parse_assignments(line[2:])
                sequences = {}
            else:
                linetype, *linecontent = line.split()
                if select_regex.match(linecontent[0]):
                    if linetype == 's':
                        seqname, start, length, strand, tot_len, seq = linecontent
                        start = int(start)
                        length = int(length)
                        #species, contigname = seqname.split('.')
                        sequences[seqname] = {'start': start,
                                              'length': length,
                                              'strand': strand,
                                              'seq': seq}
                    elif linetype == 'i':
                        seqname, leftStatus, leftCount, rightStatus, rightCount = linecontent
                        sequences[seqname].update({'leftS':  leftStatus,
                                                   'leftC':  int(leftCount),
                                                   'rightS': rightStatus,
                                                   'rightC': int(rightCount)})
                    elif linetype == 'e':
                        seqname, start, length, strand, tot_len, status = linecontent

                    yield alignment, sequences


def get_seq(maffilename, seqname_start, ignore_overlap=False):
    """Retrieve all sequences from a species (for eg), check if there are
    overlapping fragments on the reference genome, and concatenate sequences."""
    sequence, prev_seqname, prev_start, prev_length, prev_line = \
          '',         None,          0,           0,      None
    with open(maffilename) as maf:
        for lineno, line in enumerate(maf):
            if line.startswith('s ' + seqname_start):
                try:
                    _, seqname, start, length, strand, _, seq = line.split()
                except ValueError as err:
                    # Too many values to unpack, or not enough.
                    err.args += ("line %d:  %r" % (lineno, line),)
                    raise
                species, contigname = seqname.split('.')
                #sp = sequences.setdefault(species, {})
                #sp.setdefault(seqname, {})
                #if species:
                #if sp:
                #    # This assumes sequences are ordered on the reference seq.
                #    prev_start  = sp['start']
                #    prev_length = sp['length']

                start = int(start)
                length = int(length)
                gap_size =  start - (prev_start + prev_length)
                
                if seqname == prev_seqname and gap_size < 0:
                    msg = ("sequences overlap on the target (%s):\n"
                           "%d: %r\n%d: %r" % (seqname, prev_lineno, prev_line,
                                               lineno, line.rstrip()))
                    if ignore_overlap:
                        print(msg, file=sys.stderr)
                    else:
                        raise AssertionError(msg)

                #sp['seq'] += 'N'*gap_size + seq
                sequence += seq
                prev_seqname = seqname
                prev_start = start
                prev_length = length
                prev_line = line.rstrip()
                prev_lineno = lineno
    return sequence

                #else:
                #    sp['seq'] = seq
                #    sp['start'] = start
                #    # Include dash characters in the length
                #    sp['length'] = len(seq)

def get_seq2(maffilename, select_seq='', ref='hg19'):
    # for one species, contig name: (start, end, sequence)
    targets = {}

    select_pattern = select_seq + '|' + ref
    for alignment, sequences in process_maf(maffilename, select_pattern):
        refkey, = [seqname for seqname in sequences.keys() if seqname.startswith(ref)]
        refstart = sequences[refkey]['start']
        # Since there cannot be overlapping sequences in the reference seq,
        # I propose to skip this check.
        #
        refend = refstart + sequences[refkey]['length']

        for seqname, seqdata in sequences.items():
            if not seqname.startswith(ref):
                target_contigs = targets.setdefault(seqname, {})
                start, end, seq = target_contigs.get(seqname, (0, 0, ''))
                if seqdata['leftS'] == 'T':
                    print('WARNING: This region is a (tandem) duplicate. Skipping',
                          file=sys.stderr)
                elif seqdata['leftS'] == 'I':
                    gap_size = seqdata['leftC']
                    # Would be better to retrieve the actual sequence.
                    gap_filler = 'N' * gap_size
                    seq += gap_filler
                    end += gap_size
                else:
                #elif seqdata['rightS'] in ('C', 'N', 'n', 'M'):
                    pass

                seq += seqdata['seq'].replace('-', '')
                end += length



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('maffilename')
    parser.add_argument('seqname_start', 
                        help='start of the names of targets to search')
    parser.add_argument('--ignore-overlap', action='store_true', 
                        help='do not raise an error if two consecutive '\
                             'sequences overlap on the target genome.')
    
    args = parser.parse_args()
    sequence = get_seq(**vars(args))
    print('>%s from %s' % (args.seqname_start, args.maffilename))
    print(sequence.replace('-', ''))

