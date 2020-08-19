#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import argparse as ap
import os.path as op
from LibsDyogen import myPhylTree
import csv
import logging
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


# Examples
phyltreefile = op.expanduser('~/ws7/DUPLI_data93/PhylTree.TimeTree201901.Ensembl93-like.noQual.nwk')
refphyltreefile = op.expanduser('~/GENOMICUS93/PhylTree.Ensembl.93.conf')
genomequal = op.expanduser('~/GENOMICUS93/GenQualS.txt')
column_qual = 'Qual'


def main(phyltreefile, genomequal, refphyltreefile=None, column_qual='Qual'):
    phyltree = myPhylTree.PhylogeneticTree(phyltreefile)
    if refphyltreefile: refphyltree = myPhylTree.PhylogeneticTree(refphyltreefile)

    with open(genomequal, newline='') as gf:
        csvrd = csv.DictReader(gf, dialect='excel-tab')
        species_qual = {row['Species'].rstrip().replace('.', ' '): int(row[column_qual])
                        for row in csvrd}
        # Warning multiple subspecies per species with different qualities

    lstEsp2X = set()
    lstEsp6X = set()
    lstEspFull = set()
    commonNames = {}

    for sp in phyltree.listSpecies:
        try:
            q = species_qual[sp]
        except KeyError:
            logger.error('No quality found for %s', sp)
            q = 1
        if q <= 3:
            lstEspFull.add(sp)
        elif q <= 4:
            lstEsp6X.add(sp)
        else:
            lstEsp2X.add(sp)

        if refphyltreefile:
            try:
                names = refphyltree.commonNames[sp]
                commonNames[sp] = [n for n in names if isinstance(n, str)]
            except KeyError:
                logger.error('No %r common names found in reference tree', sp)

    setattr(phyltree, 'lstEsp2X', lstEsp2X)
    setattr(phyltree, 'lstEsp6X', lstEsp6X)
    setattr(phyltree, 'lstEspFull', lstEspFull)
    if refphyltreefile: setattr(phyltree, 'commonNames', commonNames)

    phyltree.printNewick(commonnames=True, symbols=True)



if __name__ == '__main__':
    logging.basicConfig()

    parser = ap.ArgumentParser(__doc__)
    parser.add_argument('phyltreefile')
    parser.add_argument('genomequal')
    parser.add_argument('refphyltreefile', nargs='?')
    parser.add_argument('-c', '--column-qual', default='Qual')

    main(**vars(parser.parse_args()))
