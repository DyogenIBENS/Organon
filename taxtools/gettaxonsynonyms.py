#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#

"""Get all alternative names of a list of taxa."""


import argparse
import fileinput
from collections import defaultdict


def load_synonyms(namesdmp, wanted='synonym'):
    names2id = {}
    id2wanted = defaultdict(list)

    with open(namesdmp) as names:
        for line in names:
            ID, name, uniq, nameclass = line.rstrip('\t|\n').split('\t|\t')
            names2id[name] = ID
            if nameclass == wanted:
                id2wanted[ID].append(name)

    return names2id, id2wanted


def main(namesdmp, taxa, wanted='synonym', oneline=False):
    #logger.debug(len(taxa), taxa)

    if not taxa:
        # read from stdin
        #taxa = []
        taxa = [line.rstrip() for line in fileinput.input([])]
        #with fileinput.input([]) as stream:
        #    while stream.isstdin():
        #        taxa.append(next(stream).rstrip())
                #logging.debug(stream.filename(), stream.fileno(), stream.isstdin())

    names2id, id2wanted = load_synonyms(namesdmp, wanted)

    if oneline:
        def printout(taxon, synonyms):
            print(taxon + '\t' + '\t'.join(synonyms))
    else:
        def printout(taxon, synonyms):
            if synonyms:
                print('\n'.join(['%s\t%s' % (taxon, syn) for syn in synonyms]))
    
    #logging.debug(len(taxa), taxa)

    for taxon in taxa:
        try:
            printout(taxon, id2wanted[names2id[taxon]])
        except KeyError as err:
            err.args = (err.args[0] + (': %s' % taxon),) + err.args[1:]
            raise


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('namesdmp')
    parser.add_argument('taxa', nargs='*')
    parser.add_argument('-w', '--wanted', default='synonym')
    parser.add_argument('--oneline', action='store_true')
    
    args = parser.parse_args()
    main(**vars(args))
