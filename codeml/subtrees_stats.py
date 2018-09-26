#!/usr/bin/env python3
# -*- coding: utf-8 -*-


"""Iterate over all subtrees and output alignment/tree/codeml statistics"""

from sys import stdin, stderr
import re
import argparse
import os.path as op
from glob import glob
import numpy as np
from scipy.stats import skew

from Bio import AlignIO
import ete3

from genomicustools.identify import SP2GENEID, convert_gene2species
from seqtools import ungap, algrep, make_al_stats
from codeml.codemlparser2 import parse_mlc
from codeml.prune2family import split_species_gene
import LibsDyogen.myPhylTree as PhylTree

ENSEMBL_VERSION = 85


def get_ensembl_ids_from_anc(ancestor, phyltree, ensembl_version=ENSEMBL_VERSION):
    return [SP2GENEID[ensembl_version][sp] for sp in phyltree.species[ancestor]]


def iter_glob_subtree_files(genetreelistfile, ancestor, filesuffix, rootdir='.',
                            subtreesdir='subtreesCleanO2', exclude=None):
    countlines = 0
    countalfiles = 0
    for line in genetreelistfile:
        countlines += 1
        genetree = line.rstrip()
        alfiles_pattern = op.join(rootdir, genetree, subtreesdir,
                                  ancestor + genetree + '*' + filesuffix)
        alfiles_reg = re.compile(re.escape(alfiles_pattern).replace('\\*', '(.*)'))
        exclude_reg = re.compile(exclude)

        print(alfiles_pattern, alfiles_reg.pattern, file=stderr)
        for subtreefile in glob(alfiles_pattern):
            if exclude is None or not exclude_reg.search(subtreefile):
                countalfiles += 1
                subtreesuffix = alfiles_reg.search(subtreefile).group(1)
                subtree = ancestor + genetree + subtreesuffix
                yield subtreefile, subtree, genetree

    print('%d lines, %d subtrees' % (countlines, countalfiles), file=stderr)


def get_al_stats(genetreelistfile, ancestor, phyltreefile, rootdir='.',
                 subtreesdir='subtreesCleanO2', ensembl_version=ENSEMBL_VERSION):
    """Gather characteristics of the **input alignments**, and output them as
    a tsv file."""

    phyltree = PhylTree.PhylogeneticTree(phyltreefile)
    ensembl_ids_anc = get_ensembl_ids_from_anc(ancestor, phyltree, ensembl_version)
    pattern = '^(' + '|'.join(ensembl_ids_anc) + ')'
    
    stats_header = ['subtree', 'genetree']
    stats_names  = [typ + '_' + measure
                    for typ in ('glob', 'mean', 'med', 'std', 'w_mean', 'w_std')
                    for measure in ('len', 'A','C','G','T', 'GC', 'N',
                                    'gaps', 'CpG')]
    stats_header += stats_names + ['ingroup_' + name for name in stats_names]

    print('\t'.join(stats_header))

    for alfile, subtree, genetree in iter_glob_subtree_files(genetreelistfile,
                                                             ancestor,
                                                             '_genes.fa',
                                                             rootdir,
                                                             subtreesdir):
        al = AlignIO.read(alfile, format='fasta')
        al = ungap(al)
        _, al_stats = make_al_stats(al)
        
        ingroup_al = ungap(algrep(al, pattern))
        _, ingroup_al_stats = make_al_stats(ingroup_al)
        
        stats_row = ['%g' % s for stat in (al_stats + ingroup_al_stats) for s in stat]

        print('\t'.join([subtree, genetree] + stats_row))
            
        #treefiles_pattern = alfiles_pattern.replace('_genes.fa', '.nwk')


def get_tree_stats(genetreelistfile, ancestor, phyltreefile, rootdir='.',
                   subtreesdir='subtreesCleanO2',
                   ensembl_version=ENSEMBL_VERSION):
    """Currently only determining the robustness of the tree.

    To find the robust trees from the given ancestor only, (excluding the
    outgroup) use `subtreesdir="subtreesClean"`."""

    phyltree = PhylTree.PhylogeneticTree(phyltreefile)
    #ensembl_ids_anc = get_ensembl_ids_from_anc(ancestor, phyltree, ensembl_version)
    print('subtree\tgenetree\troot_location\trobust_leaves\tsingle_gene_nodes')

    ancgene2sp = re.compile(r'('
                            + r'|'.join(list(phyltree.species[ancestor]) + 
                                        sorted(phyltree.getTargetsAnc(ancestor),
                                               key=lambda a:len(a),
                                               reverse=True)).replace(' ','\.')
                            + r')(.*)$')
    
    for subtreefile, subtree, genetree in iter_glob_subtree_files(genetreelistfile,
                                                              ancestor,
                                                              '.nwk',
                                                              rootdir,
                                                              subtreesdir,
                                                              exclude='_codeml\.nwk$'):
        tree = ete3.Tree(subtreefile, format=1)
        root_taxon, _ = split_species_gene(tree.name, ancgene2sp)
        
        root_location = 'O' if root_taxon is None \
                        else 'I' if root_taxon != ancestor \
                        else '='
        expected_species = phyltree.species[root_taxon or ancestor]

        if root_location == 'O':
            for node in tree.traverse('levelorder'):
                if split_species_gene(node.name, ancgene2sp)[0] is not None:
                    tree = node
                    break
        
        species_counts = {sp: 0 for sp in expected_species}
        for leaf in tree.iter_leaf_names():
            try:
                species_counts[convert_gene2species(leaf, ensembl_version)] += 1
            except KeyError as err:
                err.args += (subtreefile,)
                raise

        robust_leaves = all(c == 1 for c in species_counts.values())

        single_child_nodes = any(len(n.children) == 1 for n in tree.traverse())
        
        print('\t'.join((subtree, genetree, root_location,
                         str(int(robust_leaves)),
                         str(int(single_child_nodes))
                         )))


def get_codeml_stats(genetreelistfile, ancestor, rootdir='.',
                     subtreesdir='subtreesCleanO2'):
    """Gather characteristics of the **codeml results**, and output them as
    a tsv file."""

    stats_header = ['subtree', 'genetree']
    stats_name   = ['ls', 'ns', 'Nbranches',
                    'NnonsynSites', 'NsynSites', 'kappa',
                    'treelen', 'dN_treelen', 'dS_treelen',
                    'brlen_mean', 'brlen_std', 'brlen_med', 'brlen_skew',
                    'brOmega_mean', 'brOmega_std', 'brOmega_med', 'brOmega_skew',
                    'brdS_mean', 'brdS_std', 'brdS_med', 'brdS_skew',
                    'brdN_mean', 'brdN_std', 'brdN_med', 'brdN_skew',
                    'lnL', 'Niter', 'time used']
    
    # TODO: number of dN or dS values of zero
    br_len_reg = re.compile(r': ([0-9]+\.[0-9]+)[,)]')

    print('\t'.join(stats_header + stats_name))

    for mlcfile, subtree, genetree in iter_glob_subtree_files(genetreelistfile,
                                                              ancestor,
                                                              '_m1w04.mlc',
                                                              rootdir,
                                                              subtreesdir):
        try:
            mlc = parse_mlc(mlcfile)
            
            Nbr = len(mlc['output']['branches'])
            br_lengths = mlc['output']['branch lengths + parameters'][:Nbr]
            br_omegas = mlc['output']['omega']
            dNdS_rows = list(mlc['output']['dNdS']['rows'].values())
            dS_len = [float(x) for x in br_len_reg.findall(mlc['output']['dS tree'])]
            dN_len = [float(x) for x in br_len_reg.findall(mlc['output']['dN tree'])]
            assert len(dS_len) == Nbr
            
            #TODO: exclude the branches below the root from some values
            #      (e.g. brlen_mean)
            stats_row = [mlc['nsls']['ls'],
                         mlc['nsls']['ns'],
                         Nbr,
                         \
                         dNdS_rows[0][1],  # NnonsynSites
                         dNdS_rows[0][2],  # NsynSites
                         mlc['output']['kappa'],
                         \
                         mlc['output']['tree length'],
                         mlc['output']['tree length for dN'],
                         mlc['output']['tree length for dS'],
                         \
                         np.mean(br_lengths),    # brlen_mean
                         np.std(br_lengths),     # brlen_std
                         np.median(br_lengths),  # brlen_med
                         skew(br_lengths),       # brlen_skew
                         np.mean(br_omegas),
                         np.std(br_omegas),
                         np.median(br_omegas),
                         skew(br_omegas),
                         np.mean(dS_len),
                         np.std(dS_len),
                         np.median(dS_len),
                         skew(dS_len),
                         np.mean(dN_len),
                         np.std(dN_len),
                         np.median(dN_len),
                         skew(dN_len),
                         \
                         mlc['output']['lnL']['loglik'],
                         mlc['output']['lnL']['ntime']
                         ]

            print('\t'.join([subtree, genetree] + ['%g' % s for s in stats_row] + [mlc['Time used']]))
        except BaseException as err:
            err.args = (err.args[0] + ". File %s" % mlcfile,) + err.args[1:]
            raise


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    
    def make_subparser_func(func):
        """Transform a normal function so that it takes arguments from the Argparse args."""
        def subp_func(args):
            dictargs = vars(args)
            dictargs.pop('commands')
            dictargs.pop('func')
            return func(**dictargs)
        return subp_func

    parent_parser = argparse.ArgumentParser(add_help=False)
    parent_parser.add_argument('genetreelistfile', nargs='?',
                               type=argparse.FileType('r'), default=stdin)
    parent_parser.add_argument('ancestor')
    parent_parser.add_argument('-r', '--rootdir', default='.', help='[%(default)s]')
    parent_parser.add_argument('-s', '--subtreesdir', default='subtreesCleanO2',
                               help="[%(default)s]")
    
    subp = parser.add_subparsers(dest='commands', help='type of statistics to compile')
    
    codemlstats_parser = subp.add_parser('codeml', parents=[parent_parser])
    codemlstats_parser.set_defaults(func=make_subparser_func(get_codeml_stats))
    
    alstats_parser = subp.add_parser('alignment', parents=[parent_parser], aliases=['al'])
    alstats_parser.add_argument('phyltreefile')
    alstats_parser.add_argument('-e', '--ensembl-version', type=int,
                                default=ENSEMBL_VERSION, help="[%(default)s]")
    alstats_parser.set_defaults(func=make_subparser_func(get_al_stats))
    
    treestats_parser = subp.add_parser('tree', parents=[parent_parser], aliases=['tr'])
    treestats_parser.add_argument('phyltreefile')
    treestats_parser.add_argument('-e', '--ensembl-version', type=int,
                                  default=ENSEMBL_VERSION, help="[%(default)s]")
    treestats_parser.set_defaults(func=make_subparser_func(get_tree_stats))

    args = parser.parse_args()
    args.func(args)

