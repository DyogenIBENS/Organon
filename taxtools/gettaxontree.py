#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#

"""Extract the taxon taxonomic tree from the ncbi taxonomy files
(nodes.dmp and names.dmp)"""

import re
from collections import defaultdict
import ete3
import argparse
import logging
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
#logging.basicConfig(format='%(levelname)s:%(funcName)s:%s(message)s',
#                    level=logging.INFO)



EXCLUDE=r'\b(sp|spp|x)\b'


def load_names2id(namesdmp, wantedname='scientific name'):
    names2id = {}
    id2names = {}
    with open(namesdmp) as names:
        for line in names:
            ID, name, uniq, nameclass = line.rstrip('\t|\n').split('\t|\t')
            names2id[name] = ID
            if nameclass == wantedname or not id2names.get(ID):
                id2names[ID] = name
    return names2id, id2names


def load_descent(nodesdmp, cut='subspecies'):
    descent = defaultdict(set)
    with open(nodesdmp) as nodes:
        for line in nodes:
            tax_id, parent_tax_id, rank, *_ = line.split('\t|\t')
            if rank != cut:
                descent[parent_tax_id].add(tax_id)
    return descent


def build_nodes_ete3(descent, id2names, node, ID, exclude=EXCLUDE):
    children = descent.get(ID, set())
    for child in children:
        childname = id2names.get(child)
        if not exclude or not re.search(exclude, childname):
            nextnode = node.add_child(name=childname)
            build_nodes_ete3(descent, id2names, nextnode, child, exclude)

def build_tree_ete3(descent, id2names, names2id, taxon, exclude=EXCLUDE):
    tree = ete3.Tree(name=taxon)
    root_id = names2id[taxon]
    build_nodes_ete3(descent, id2names, tree, root_id, exclude)
    return tree


def extractree(nodesdmp, namesdmp, taxon, cut='subspecies',
               wantedname='scientific name', exclude=EXCLUDE):
    logger.info("Loading name-to-id conversion")
    names2id, id2names = load_names2id(namesdmp, wantedname)
    logger.info("Loading node children")
    descent = load_descent(nodesdmp, cut)
    logger.info("Building tree")
    return build_tree_ete3(descent, id2names, names2id, taxon, exclude)


def main(nodesdmp, namesdmp, taxon, outfile=None, cut='subspecies',
         wantedname='scientific name', exclude=EXCLUDE):
    tree = extractree(nodesdmp, namesdmp, taxon, cut, wantedname, exclude)
    print("Writing tree of size %d" % len(tree))
    if not outfile:
        #for fmt in [3, 5, 6, 7, 8]:
        #    print('---')
        print(tree.write(format=8, format_root_node=True))
    else:
        tree.write(outfile=outfile, format=8, format_root_node=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('nodesdmp')
    parser.add_argument('namesdmp')
    parser.add_argument('taxon')
    parser.add_argument('-o', '--outfile', help='default to stdout')
    parser.add_argument('-c', '--cut', default='subspecies')
    parser.add_argument('-w', '--wantedname', default='scientific name',
                        choices=['scientific name', 'common name', 'authority',
                                 'genebank common name', 'synonym'])
    parser.add_argument('--exclude', default=EXCLUDE,
                        help='Remove nodes matching [%(default)r]')
    
    args = parser.parse_args()
    main(**vars(args))

