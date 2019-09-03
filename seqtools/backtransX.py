#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Written to circumvent `treebest backtrans` which *deletes* sequences
containing "X" without any warning.
"""


from sys import stdout
import argparse as ap
from Bio import SeqIO, SeqRecord, Seq
from collections import OrderedDict

from seqtools import IUPAC
#from seqtools.plot_al_conservation import al2array  # ~~> move

from Bio import Alphabet  # Actually all the non empty alphabets are in Alphabet.IUPAC
#from Bio.Tools import Translate  # Doesn't exist.
import logging
logger = logging.getLogger(__name__)


DNAGAP = '-'
FORMAT = 'fasta'


def locate_differences(s1, s2):
    if len(s1) != len(s2):
        return None
    return [i for i, (x1, x2) in enumerate(zip(s1, s2)) if x1 != x2]


def display_translated_diff(dna, s1, s2):
    positions = locate_differences(s1, s2)
    if positions is None:
        return 'lengths: %d != %d' % (len(s1), len(s2))
    return '\n'.join(' '.join(row)
                     for row in zip(*(('%5d' % p,
                            '%5s' % dna[p*3:(p+1)*3],
                            '%5s' % s1[p],
                            '%5s' % s2[p]) for p in positions)))


def backtrans(prots: dict, dnas: dict):
    #TODO: Extra checks:
    # - take `prots` as a MultipleAlignment object.
    # - verify that `out_cds` can be converted into a MultipleAlignment (compatible lengths)
    out_cds = []
    for name, prot in prots.items():
        dnaseq = str(dnas[name].seq).replace(DNAGAP, '')  #.replace('.', '')
        assert not len(dnaseq) % 3
        i = 0
        backseq = ''
        for aa in prot.seq:
            if aa in IUPAC.gaps:
                backseq += '---'
            else:
                backseq += dnaseq[i*3:(i+1)*3]
                i += 1
        backdna = SeqRecord.SeqRecord(Seq.Seq(backseq, Alphabet.IUPAC.ambiguous_dna),
                                      id=name, name=name, description=prot.description)
        # Check
        retranslated = str(backdna.translate(gap='-').seq)
        try:
            assert retranslated == str(prot.seq)
        except AssertionError as err:
            err.args = ('At %s:\n%s' % (name,
                            display_translated_diff(str(backdna.seq),
                                                    retranslated,
                                                    str(prot.seq))),)
            if all('X' in (retranslated[p], prot.seq[p])
                   for p in locate_differences(retranslated, str(prot.seq))):
                logger.warning('Inconsistent treatment of ambiguous codons. ' +
                               err.args[0])
            else:
                raise

        # This works even with 'NNN' codons corresponding to 'X'.
        out_cds.append(backdna)

    return out_cds


def backtransIO(inputprot, inputdna, outfile, format=FORMAT):
    prots = OrderedDict((seq.id, seq)
                        for seq in SeqIO.parse(inputprot, format,
                            Alphabet.IUPAC.extended_protein))
    dnas = {seq.id: seq
            for seq in SeqIO.parse(inputdna, format,
                                   Alphabet.IUPAC.ambiguous_dna)}

    SeqIO.write(backtrans(prots, dnas), outfile, format)


if __name__ == '__main__':
    logging.basicConfig(format=logging.BASIC_FORMAT)
    parser = ap.ArgumentParser(description=__doc__)
    parser.add_argument('inputprot', help='Typically an alignment file')
    parser.add_argument('inputdna', help='Gaps are ignored.')
    parser.add_argument('outfile', nargs='?', default=stdout,
                        type=ap.FileType('w'))
    parser.add_argument('-f', '--format', default=FORMAT, 
                        help='Input and output sequence format [%(default)s]')
    
    args = parser.parse_args()
    backtransIO(**vars(args))

