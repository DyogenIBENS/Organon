#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import print_function

"""Convert the .nuc format of sequences (from the PAML example data) into fasta."""


import argparse
import re


def main(nucfile, outfile):
    reg_label = re.compile(r'[0-9]+\.\S+$')
    reg_seq   = re.compile(r'[ATGC ]+[0-9]+$')
    with open(nucfile) as nuc:
        line0 = nuc.readline()
        line = nuc.readline()
        stats = line0.split()
        assert len(stats) == 3 and not line.rstrip()
        stats[0] = int(stats[0])
        stats[1] = int(stats[1])

        sequences = []
        while line:
            line = line.rstrip()
            if line:
                if reg_label.match(line):
                    sequences.append([line, ''])
                elif reg_seq.match(line):
                    sequences[-1][1] += ''.join(line.split()[:-1])
                else:
                    break
            line = nuc.readline()

    try:
        assert len(sequences) == stats[0]
    except AssertionError:
        print(stats[0], len(sequences))
        print(sequences)
        raise
    #print(stats[0], len(sequences))
    #print(sequences)
    
    assert all(len(s)==stats[1] for _, s in sequences)

    with open(outfile, 'w') as out:
        for label, seq in sequences:
            out.write('>%s\n%s\n' % (label, seq))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('nucfile')
    parser.add_argument('outfile')
    
    args = parser.parse_args()
    main(**vars(args))
