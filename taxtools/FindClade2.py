#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import sys

from itertools import chain

import argparse
from ete3 import NCBITaxa


def findclade(namelist, ranks='family|genus'):
    #rankregex = re.compile('^(%s)$' % ranks)
    ncbi = NCBITaxa()
    name2taxid = ncbi.get_name_translator(namelist)
    lineages = ncbi.get_lineage_translator([v[0] for v in name2taxid.values()])
    cladetaxids = []
    for name in namelist:
        lineage = lineages[name2taxid[name][0]]
        #print(name, name2taxid[name], lineage)
        rank2clade = {rk: taxid for taxid, rk in ncbi.get_rank(lineage).items()}
        cladetaxids.append([rank2clade.get(rank, 0) for rank in ranks.split('|')])

    #print(cladetaxids)
    taxid2clade = ncbi.get_taxid_translator(chain(*cladetaxids))

    for name, taxidlist in zip(namelist, cladetaxids):
        yield name, [taxid2clade.get(t, '') for t in taxidlist]


def main(infile, ranks='family|genus'):
    #stream = open(filename) if filename else sys.stdin

    print('\t' + ranks.replace('|', '\t'))
    for name, clades in findclade([line.rstrip() for line in infile], ranks):
        print(name + '\t' + '\t'.join(clades))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('infile', nargs='?', type=argparse.FileType('r'),
                        default=sys.stdin)
    parser.add_argument('-r', '--ranks', default='family|genus')
    
    args = parser.parse_args()
    main(**vars(args))
