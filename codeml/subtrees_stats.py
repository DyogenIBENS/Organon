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

import LibsDyogen.myPhylTree as PhylTree

from genomicustools.identify import SP2GENEID, \
                                    convert_gene2species
from dendron.reconciled import get_taxon, \
                               get_taxon_treebest, \
                               infer_gene_event_taxa
from seqtools import ungap, \
                     algrep, \
                     make_al_stats
from codeml.codemlparser2 import parse_mlc
from codeml.prune2family import split_species_gene


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
        exclude_reg = exclude if exclude is None else re.compile(exclude)

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

        
def simple_robustness_test(tree, expected_species, ensembl_version=ENSEMBL_VERSION):
    """Determine robustness after the expected set of species at the tips."""
    
    species_counts = {sp: 0 for sp in expected_species}
    for leaf in tree.iter_leaf_names():
        species_counts[convert_gene2species(leaf, ensembl_version)] += 1

    leaves_robust = all(c == 1 for c in species_counts.values())

    single_child_nodes = any(len(n.children) == 1 for n in tree.traverse())
    
    return (leaves_robust, single_child_nodes)


def any_get_taxon(node, ancgene2sp, ensembl_version=ENSEMBL_VERSION):
    try:
        return get_taxon(node, ancgene2sp, ensembl_version=ENSEMBL_VERSION)
    except ValueError:
        return get_taxon_treebest(node)


def per_node_events(tree, phyltree, aberrant_dist=10000,
                    get_taxon=any_get_taxon, *args):
    """Check each node and its type (duplication, deletion, speciation).
    """
    aberrant_dists = 0
    wrong_duptypes = 0  # NotImplemented (where duplication != 0 or 2 or 3)

    leaves = 0

    events = ('dup', 'spe', 'speloss', 'duploss')
    # TreeBest and taxa agree.
    sure_events = {evt: 0 for evt in events}

    only_treebest_events = {evt: 0 for evt in events}

    for node in tree.traverse('preorder'):
        if node.dist >= aberrant_dist:
            aberrant_dists += 1

        if node.is_leaf():
            leaves += 1
            continue

        treebest_isdup = getattr(node, 'D', None)

        taxon = get_taxon(node, *args)
        expected_speciated_taxa = set(t[0] for t in phyltree.items.get(taxon, []))
        children_taxa = set(get_taxon(ch, *args) for ch in node.children)

        event = infer_gene_event_taxa(node, taxon, children_taxa)

        assert event != 'leaf'

        counter = sure_events
        if event == 'dup' and treebest_isdup == 'N':
            counter = only_treebest_events
            event = 'spe'
        elif event == 'spe' and treebest_isdup == 'Y':
            counter = only_treebest_events
            event = 'spe'

        if (event == 'dup' and counter is sure_events and
               len(node.children) < 2) or \
           (event == 'spe' and len(expected_speciated_taxa)>0 and \
               len(node.children) < len(expected_speciated_taxa)):
            event += 'loss'

        counter[event] += 1

    #return tuple(sure_events[e] for e in events) + \
    #       tuple(only_treebest_events[e] for e in events)
    return sure_events, only_treebest_events, aberrant_dists


def make_ancgene2sp(ancestor, phyltree):
    return re.compile(r'('
                      + r'|'.join(list(phyltree.species[ancestor]) +
                                  sorted(phyltree.getTargetsAnc(ancestor),
                                         key=lambda a:len(a),
                                         reverse=True)).replace(' ','\.')
                      + r')(.*)$')


def get_tree_stats(genetreelistfile, ancestor, phyltreefile, rootdir='.',
                   subtreesdir='subtreesCleanO2',
                   ensembl_version=ENSEMBL_VERSION,
                   ignore_outgroups=False, extended=False):
    """Currently only determining the robustness of the tree.

    To find the robust trees from the given ancestor only, (excluding the
    outgroup) use `subtreesdir="subtreesClean"`."""

    phyltree = PhylTree.PhylogeneticTree(phyltreefile)
    #ensembl_ids_anc = get_ensembl_ids_from_anc(ancestor, phyltree, ensembl_version)
    print('subtree\tgenetree\troot_location\tleaves_robust\tsingle_gene_nodes'
            + ('\tnodes_robust\tonly_treebest_spe\taberrant_dists' if extended else ''))

    ancgene2sp = make_ancgene2sp(ancestor, phyltree)
    all_ancgene2sp = make_ancgene2sp(phyltree.root, phyltree)
    
    for subtreefile, subtree, genetree in iter_glob_subtree_files(genetreelistfile,
                                                              ancestor,
                                                              '.nwk',
                                                              rootdir,
                                                              subtreesdir,
                                                              exclude='_codeml\.nwk$'):
        tree = ete3.Tree(subtreefile, format=1)
        root_taxon, _ = split_species_gene(tree.name, ancgene2sp)
        
        # Determine if root_taxon is inside (I) or outside (O) of the clade.
        root_location = 'O' if root_taxon is None \
                        else 'I' if root_taxon != ancestor \
                        else '='

        if root_location == 'O':
            if ignore_outgroups:
                for node in tree.traverse('levelorder'):
                    if split_species_gene(node.name, ancgene2sp)[0] is not None:
                        tree = node
                        break
            else:
                root_taxon, _ = split_species_gene(tree.name, all_ancgene2sp)

        expected_species = phyltree.species[root_taxon or ancestor]
        try:
            leaves_robust, single_child_nodes = simple_robustness_test(tree,
                                                            expected_species)
        except KeyError as err:  # Error while converting genename to species
            err.args += (subtreefile,)
            raise
        output = (int(leaves_robust), int(single_child_nodes))

        if extended:
            sure_events, only_treebest_events, aberrant_dists = \
                    per_node_events(tree,
                                    phyltree,
                                    10000,
                                    any_get_taxon,
                                    all_ancgene2sp,
                                    ensembl_version)
            nodes_robust = not any((sure_events['dup'],
                                    sure_events['speloss'],
                                    sure_events['duploss'],
                                    only_treebest_events['dup'],
                                    only_treebest_events['speloss'],
                                    only_treebest_events['duploss']))
            output += (int(nodes_robust), only_treebest_events['spe'], aberrant_dists)


        print('\t'.join((subtree, genetree, root_location) +
                         tuple(str(x) for x in output)))


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
    treestats_parser.add_argument('-E', '--extended', action='store_true',
                                  help='Perform robustness test on nodes (instead of leaves VS root)')
    treestats_parser.set_defaults(func=make_subparser_func(get_tree_stats))

    args = parser.parse_args()
    args.func(args)

