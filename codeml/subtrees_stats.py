#!/usr/bin/env python3
# -*- coding: utf-8 -*-


"""Iterate over all subtrees and output alignment/tree/codeml statistics"""

from sys import stdin, stderr
import re
import argparse
import os.path as op
from glob import glob
import numpy as np
from Bio import AlignIO
from codeml.select_leaves_from_specieslist import SP2GENEID
from seqtools import ungap, algrep, compo_freq
from codeml import codemlparser
import LibsDyogen.myPhylTree as PhylTree

ENSEMBL_VERSION = 85

def get_ensembl_ids_from_anc(ancestor, phyltree, ensembl_version=ENSEMBL_VERSION):
    return [SP2GENEID[ensembl_version][sp] for sp in phyltree.species[ancestor]]

def iter_glob_subtree_files(genetreelistfile, ancestor, filesuffix, rootdir='.',
                             subtreesdir='subtreesCleanO2'):
    countlines = 0
    countalfiles = 0
    for line in genetreelistfile:
        countlines += 1
        genetree = line.rstrip()
        alfiles_pattern = op.join(rootdir, genetree, subtreesdir,
                                  ancestor + genetree + '*' + filesuffix)
        alfiles_reg = re.compile(re.escape(alfiles_pattern).replace('\\*', '(.*)'))
        
        print(alfiles_pattern, alfiles_reg.pattern, file=stderr)
        for subtreefile in glob(alfiles_pattern):
            countalfiles += 1
            subtreesuffix = alfiles_reg.search(subtreefile).group(1)
            subtree = ancestor + genetree + subtreesuffix
            yield subtreefile, subtree, genetree

    print('%d lines, %d subtrees' % (countlines, countalfiles), file=stderr)


def make_al_stats(genetreelistfile, ancestor, phyltreefile, rootdir='.',
         subtreesdir='subtreesCleanO2', ensembl_version=ENSEMBL_VERSION):
    """Gather characteristics of the **input alignments**, and output them as
    a tsv file."""

    phyltree = PhylTree.PhylogeneticTree(phyltreefile)
    ensembl_ids_anc = get_ensembl_ids_from_anc(ancestor, phyltree, ensembl_version)
    pattern = '^(' + '|'.join(ensembl_ids_anc) + ')'
    
    stats_header = ['subtree', 'genetree']
    stats_names  = [typ + '_' + measure \
                        for typ in ('glob', 'mean', 'med', 'std') \
                            for measure in ('len', 'GC', 'N', 'gaps', 'CpG')]
    stats_header += stats_names + ['ingroup_' + name for name in stats_names]

    print('\t'.join(stats_header))

    for alfile, subtree, genetree in iter_glob_genetree_files(genetreelistfile,
                                                              ancestor,
                                                              '_genes.fa',
                                                              root_dir,
                                                              subtrees_dir):
        al = AlignIO.read(alfile, format='fasta')
        al = ungap(al)
        al_stats = compo_freq(al)
        
        ingroup_al = ungap(algrep(al, pattern))
        ingroup_al_stats = compo_freq(ingroup_al)
        
        stats_row = ['%g' % s for stat in (al_stats + ingroup_al_stats) for s in stat]

        print('\t'.join([subtree, genetree] + stats_row))
            
        #treefiles_pattern = alfiles_pattern.replace('_genes.fa', '.nwk')


def make_codeml_stats(genetreelistfile, ancestor, phyltreefile, rootdir='.',
         subtreesdir='subtreesCleanO2', ensembl_version=ENSEMBL_VERSION):
    """Gather characteristics of the **codeml results**, and output them as
    a tsv file."""

    print('\t'.join(stats_header))

    stats_header = ['subtree', 'genetree']
    stats_name   = ['ls', 'ns', 'Nbranches', 'treelen', 'brlen_mean', 'brlen_std', 'kappa',
                    'dN_treelen', 'dS_treelen', 'brdS_mean', 'brdS_std', 'time used']
    
    print('\t'.join(stats_header + stats_name))

    for mlcfile, subtree, genetree in iter_glob_genetree_files(genetreelistfile,
                                                              ancestor,
                                                              '_m1w04.mlc',
                                                              root_dir,
                                                              subtrees_dir):
        mlc = codemlparser.parse_mlc(mlcfile)
        stats_row = [mlc['nsls']['ls'],
                     mlc['nsls']['ns'],
                     len(mlc['output']['lnL']['branches']),
                     mlc['output']['lnL']['tree length'],
                     np.mean(mlc['output']['lnL']['branch_lengths']),
                     np.std(mlc['output']['lnL']['branch_lengths']),
                     mlc['output']['kappa'],
                     mlc['output']['tree length for dN'],
                     mlc['output']['tree length for dS'],
                     np.mean(mlc['output'][]),
                     np.std(
                     ]
        print('\t'.join([subtree, genetree] + stats_row))



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
    make_al_stats(**vars(args))

