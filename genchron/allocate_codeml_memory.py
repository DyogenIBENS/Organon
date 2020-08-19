#!/usr/bin/env python3


"""
Draft script, not actually used in the pipelines.


Print the estimated memory (max) needed for running codeml with this\
 alignment.
USAGE:
    ./allocate_codeml_memory.py phylipfile
    ./allocate_codeml_memory.py -n phylipfile # like echo -n, no newline.
EXAMPLE:
    condor_descript.py condorjob.txt codeml -a gene.ctl\
 -m `./allocate_codeml_memory.py gene.phy`
 """

from __future__ import print_function

import sys


def allocate_memory(phylip_file):
    """Decide how much memory should be allocated to the codeml run on the
    cluster"""
    with open(phylip_file) as IN:
        line = IN.readline()
        fields = line.split()
        nseq = int(fields[0])
        length = int(fields[1])
    if nseq < 50:
        return "500M"
    elif nseq < 100:
        return "1G"
    elif nseq < 200:
        return "2G"
    else:
        return "10G"


if __name__=='__main__':
    endnewline = True
    if len(sys.argv) != 2:
        if len(sys.argv) == 3 and sys.argv[1] == '-n':
            sys.argv.remove('-n')
            endnewline = False
        else:
            print(__doc__, file=sys.stderr)
            sys.exit(2)

    mem = allocate_memory(sys.argv[1])
    if endnewline: mem += '\n'
    print(mem, end='')
