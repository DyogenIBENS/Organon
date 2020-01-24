#!/usr/bin/env python3
# -*- coding: utf-8 -*-


"""Read output file from ALE (.uml_rec: samples of reconciled gene trees), and 
reformat them as for TreeBest.

A new NHX tag is needed: T= for transfer.
"""

import argparse as ap
import re
import ete3



def parse_ALEoutput(aleoutputfile):
    #TODO: archiparse
    outputs = {}
    with open(aleoutputfile) as aleoutput:
        header = next(aleoutput).rstrip()
        assert header.startswith('#') and 'ALE' in header

        line = next(aleoutput).rstrip()
        while not line:
            line = next(aleoutput).rstrip()

        assert line.startswith('S:\t'), line
        # Found the species tree
        outputs['Stree'] = line

        line = next(aleoutput).rstrip()
        while not line:
            line = next(aleoutput).rstrip()

        assert line.startswith('Input ale from:'), line
        outputs['ale input'] = line.split('\t', 1)

        line = next(aleoutput).rstrip()
        assert line.startswith('>logl:'), line
        outputs['logl'] = float(line.split(':')[1])

        line = next(aleoutput).rstrip()
        assert line == 'rate of\t Duplications\tTransfers\tLosses', repr(line)

        line = next(aleoutput).rstrip()
        assert line.startswith('ML \t'), repr(line)
        outputs['rate_Duplications'], \
        outputs['rate_Transfers'], \
        outputs['rate_Losses'] = [float(r) for r in line.split('\t')[1:]]

        line = next(aleoutput).rstrip()
        try:
            Nrec = int(re.compile(r'(\d+) reconciled G-s:').search(line).group(1))
        except AttributeError as err:
            err.args += ('at line %r' % line,)
            raise
        line = next(aleoutput).rstrip()
        while not line:
            line = next(aleoutput).rstrip()

        outputs['reconciliations'] = [line]
        for i in range(Nrec-1):
            outputs['reconciliations'].append(next(aleoutput).rstrip())

        line = next(aleoutput).rstrip()
        assert line.startswith('# of\t Duplications\tTransfers\tLosses\tSpeciations'), repr(line)
        line = next(aleoutput).rstrip()
        assert line.startswith('Total \t'), repr(line)
        outputs['total_Duplications'], \
        outputs['total_Transfers'], \
        outputs['total_Losses'], \
        outputs['total_Speciations'] = [float(n) for n in line.split('\t')[1:]]

        line = next(aleoutput).rstrip()
        while not line:
            line = next(aleoutput).rstrip()

        table = {}
        assert line.startswith('# of\t Duplications\tTransfers\tLosses\tOriginations\tcopies\tsingletons\textinction_prob'), repr(line)
        table['header'] = line.split('\t')[1:]
        for line in aleoutput:
            fields = line.rstrip().split('\t')
            table[fields[1]] = [float(x) for x in fields[2:]]
        outputs['branch_stats'] = table

        return outputs

def ale_species_numbers_2_names(ale_Stree, labelled_tree):
    pass

def ale2treebest(tree):
    pass




