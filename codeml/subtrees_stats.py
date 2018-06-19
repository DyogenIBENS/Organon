#!/usr/bin/env python3
# -*- coding: utf-8 -*-


"""Iterate over all subtrees and output alignment/tree/codeml statistics"""

from sys import stdin, stderr
import re
import argparse
import os.path as op
from glob import glob
from Bio import AlignIO
from codeml.select_leaves_from_specieslist import SP2GENEID
from seqtools import ungap, algrep, compo_freq
import LibsDyogen.myPhylTree as PhylTree

ENSEMBL_VERSION = 85

def get_ensembl_ids_from_anc(ancestor, phyltree, ensembl_version=ENSEMBL_VERSION):
    return [SP2GENEID[ensembl_version][sp] for sp in phyltree.species[ancestor]]


def main(genetreelistfile, ancestor, phyltreefile, rootdir='.',
         subtreesdir='subtreesCleanO2', ensembl_version=ENSEMBL_VERSION):

    phyltree = PhylTree.PhylogeneticTree(phyltreefile)
    ensembl_ids_anc = get_ensembl_ids_from_anc(ancestor, phyltree, ensembl_version)
    pattern = '^(' + '|'.join(ensembl_ids_anc) + ')'
    
    stats_header = ['subtree', 'genetree']
    stats_names  = [typ + '_' + measure \
                        for typ in ('glob', 'mean', 'med', 'std') \
                            for measure in ('len', 'GC', 'N', 'gaps', 'CpG')]
    stats_header += stats_names + ['ingroup_' + name for name in stats_names]

    print('\t'.join(stats_header))

    countlines = 0
    countalfiles = 0
    for line in genetreelistfile:
        countlines += 1
        genetree = line.rstrip()
        alfiles_pattern = op.join(rootdir, genetree, subtreesdir,
                                  ancestor + genetree + '*_genes.fa')
        alfiles_reg = re.compile(re.escape(alfiles_pattern).replace('\\*', '(.*)'))
        
        print(alfiles_pattern, alfiles_reg.pattern, file=stderr)
        for alfile in glob(alfiles_pattern):
            countalfiles += 1
            subtreesuffix = alfiles_reg.search(alfile).group(1)
            subtree = ancestor + genetree + subtreesuffix

            al = AlignIO.read(alfile, format='fasta')
            al = ungap(al)
            al_stats = compo_freq(al)
            
            ingroup_al = ungap(algrep(al, pattern))
            ingroup_al_stats = compo_freq(ingroup_al)
            
            stats_row = ['%g' % s for stat in (al_stats + ingroup_al_stats) for s in stat]

            print('\t'.join([subtree, genetree] + stats_row))
            
        #treefiles_pattern = alfiles_pattern.replace('_genes.fa', '.nwk')
    print('%d lines, %d subtrees' % (countlines, countalfiles), file=stderr)

            
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('genetreelistfile', nargs='?',
                        type=argparse.FileType('r'), default=stdin)
    parser.add_argument('ancestor')
    parser.add_argument('phyltreefile')
    parser.add_argument('-e', '--ensembl-version', type=int, default=ENSEMBL_VERSION,
                        help="[%(default)s]")
    parser.add_argument('-r', '--rootdir', default='.')
    parser.add_argument('-s', '--subtreesdir', default='subtreesCleanO2',
                        help="[%(default)s]")
    
    args = parser.parse_args()
    main(**vars(args))

