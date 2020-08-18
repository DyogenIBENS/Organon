#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
USAGE:
    ./find_non_overlapping_codeml_results.py <<<"mlcfile [mlcfile...]"

Split genes can be detected in alignments when two or more gene from the same species in the same alignment do not overlap at all.

It can be detected in the `.mlc` file from codeml in the Nei-Gojobori matrix, when the distance between two sequences is -1.0000.
"""

import sys
import os
import fileinput
import re
import os.path as op
from pamliped.codemlparser2 import parse_mlc


# space preceded by closing parenthesis
#RE_COL_SEP = re.compile('(?<=\)) ')
#COL_WIDTH = 23
RE_FIELD = re.compile('-?\d\.\d+ \(-?\d+\.\d+ -?\d\.\d+\)')
NAval = '-1.0000 (-1.0000 -1.0000)'

def parse_NG(mlc):
    NGlines = mlc['Nei & Gojobori']['matrix']
    assert len(NGlines[0]) == len(NGlines[1]), "Bad parsing of NG matrix."
    seqnames, nblines = zip(*NGlines)
    return seqnames, [RE_FIELD.findall(line) for line in nblines]


def quickcheck_NG(mlc):
    seqnames, values = parse_NG(mlc)
    for i, linevalues in enumerate(values):
        if NAval in linevalues:
            # There is no intersection between the two sequence positions.
            return (seqnames[i], seqnames[linevalues.index(NAval)])
            # Return now and omit any other problematic values

def list_nonoverlapping_NG(mlc):
    seqnames, values = parse_NG(mlc)
    nonoverlapping = []
    for i, linevalues in enumerate(values):
        for j, cellval in enumerate(linevalues):
            if cellval == NAval:
                # There is no intersection between the two sequence positions.
                nonoverlapping.append((seqnames[i], seqnames[j]))
                # Return now and omit any other problematic values
    return nonoverlapping


def main(mlcfiles):
    for mlcfile in mlcfiles:
        r = quickcheck_NG(parse_mlc(mlcfile))
        if r:
            print('%s: %s - %s' % (op.basename(mlcfile), *r))


if __name__ == '__main__':
    if len(sys.argv) == 1 and os.isatty(0):
        print(__doc__, file=sys.stderr)
        sys.exit()
    main((line.rstrip() for line in fileinput.input() if not line.startswith('#')))

