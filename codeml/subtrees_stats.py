#!/usr/bin/env python3
# -*- coding: utf-8 -*-


"""Iterate over all subtrees and output alignment/tree/codeml statistics"""

import argparse
import os.path as op
from glob import glob
from Bio import AlignIO
from codeml.select_leaves_from_specieslist import SP2GENEID
from seqtools import ungap, seqrecords_grep, compo_freq
import LibsDyogen.myPhylTree as PhylTree


def get_ensembl_ids_from_anc(ancestor, phyltree):
    return [SP2GENEID[sp] for sp in phyltree.species(ancestor)]


def main(genetreelistfile, ancestor, rootdir='.', subtreesdir='subtreesCleanO2',
         phyltreefile):

    phyltree = PhylTree.PhylogeneticTree(phyltreefile)
    ensembl_ids_anc = get_ensembl_ids_from_anc(ancestor)
    pattern = '^(' + '|'.join(ensembl_ids_anc) + ')'
    
    stats_header = ['subtree', 'genetree']
    stats_names  = [typ + '_' + measure \
                        for typ in ('glob', 'mean', 'med', 'std') \
                            for measure in ('len', 'GC', 'N', 'gaps', 'CpG')]
    stats_header += stats_names + ['ingroup_' + name for name in stats_names]

    print('\t'.join(stats_header))

    with open(genetreelistfile) as genetreelines:
        for line in genetreelines:
            genetree = line.rstrip()
            alfiles_pattern = op.join(rootdir, subtreesdir,
                                      ancestor + genetree + '*_genes.fa')
            alfiles_reg = re.compile(re.escape(alfiles_pattern).replace(r'\\*', '(.*)'))
            
            for alfile in glob(alfiles_pattern):
                subtreesuffix = alfiles_reg.search(alfile).group(1)
                subtree = ancestor + genetree + subtreesuffix

                al = AlignIO.read(alfile, format='fasta')
                al = ungap(al)
                al_stats = compo_freq(al)
                
                ingroup_al = ungap(seqrecords_grep(al, pattern))
                ingroup_al_stats = compo_freq(ingroup_al)
                
                stats_row = [s for stat in al_stats for s in stat] + \
                            [s for stat in ingroup_al_stats for s in stat]

                print('\t'.join([subtree, genetree] + stats_row))
                

            
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('genetreelistfile')
    parser.add_argument('ancestor')
    parser.add_argument('phyltreefile')
    parser.add_argument('-r', '--rootdir', default='.')
    parser.add_argument('-s', 'subtreesdir', default='subtreesCleanO2')
    
    args = parser.parse_args()
    main(**vars(args))

