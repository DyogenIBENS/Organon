#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import argparse as ap
from LibsDyogen import myPhylTree
import ete3
from dendro.parsers import read_multinewick
from dendro.trimmer import thin_ete3 as thin
from genomicustools.identify import ultimate_seq2sp

ENSEMBL_VERSION = 93


def main(treefile, speciesfile, ensembl_version, from_phyltree=False, 
         outfile=None, keep_single_node_trees=False):
    if from_phyltree:
        phyltree = myPhylTree.PhylogeneticTree(speciesfile)
        specieslist = phyltree.listSpecies
    else:
        with open(speciesfile) as f:
            specieslist = [line.rstrip() for line in f if not line.startswith('#')]

    outtrees = []

    with open(treefile) as newick:
        for tree in read_multinewick(newick):
            tree = ete3.Tree(tree, format=1)
            features = set.union(*(n.features for n in tree.traverse())) \
                        - set(('name', 'dist', 'support'))

            if from_phyltree:
                keptleaves = [l for l in tree.iter_leaves()
                              if (ultimate_seq2sp(l.name, ensembl_version)
                                  in specieslist)]
            else:
                keptleaves = [l for l in tree.iter_leaves()
                              if (ultimate_seq2sp(l.name, ensembl_version)
                                  not in specieslist)]
            newtree = thin(tree, keptleaves)
            if newtree and (len(newtree)>1 or keep_single_node_trees):
                newnewick = tree.write(outfile=None, format=1, format_root_node=True)
                outtrees.append(newnewick)

    if outfile is not None and outtrees:
        outfile = open(outfile, 'w')
    if outtrees:
        print('\n'.join(outtrees), file=outfile)
    if outfile is not None: outfile.close()


if __name__ == '__main__':
    parser = ap.ArgumentParser(description=__doc__)
    parser.add_argument('treefile')
    parser.add_argument('speciesfile', help='Species to exclude (except if taken from a PhylTree)')
    parser.add_argument('-o', '--outfile')
    parser.add_argument('-e', '--ensembl-version', type=int,
                        default=ENSEMBL_VERSION)
    parser.add_argument('-p', '--from-phyltree', action='store_true')
    parser.add_argument('-1', '--keep-single-node-trees', action='store_true')
    
    args = parser.parse_args()
    main(**vars(args))
