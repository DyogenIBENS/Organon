#!/usr/bin/env python3
# -*- coding: utf-8 -*-


"""Iterate over all subtrees and output alignment/tree/codeml statistics"""

from sys import stdin
import re
import argparse
import os.path as op
from glob import glob
import numpy as np
from scipy.stats import skew

from Bio import AlignIO
import ete3

import LibsDyogen.myPhylTree as PhylTree

from UItools.autoCLI import make_subparser_func
from genomicustools.identify import SP2GENEID, \
                                    convert_gene2species
from dendro.reconciled import get_taxon, \
                              get_taxon_treebest, \
                              infer_gene_event_taxa
from dendro.bates import iter_distleaves
from dendro.trimmer import fuse_single_child_nodes_ete3
from seqtools import ungap, \
                     algrep, \
                     make_al_compo
from seqtools.plot_al_conservation import reorder_al, get_position_stats, parsimony_score
from codeml.codemlparser2 import parse_mlc
from codeml.find_non_overlapping_codeml_results import list_nonoverlapping_NG
from codeml.prune2family import split_species_gene
from seqtools.compo_freq import get_seq_counts
from seqtools.gblocks_parser import parse_gb_html
from seqtools.fillpositions import parse_seqranges

import logging
logger = logging.getLogger(__name__)

ENSEMBL_VERSION = 85


# Avoid the numpy warning when the array is empty
def mean(a):
    return np.mean(a) if len(a) else np.NaN

def std(a):
    return np.std(a) if len(a) else np.NaN

def median(a):
    return np.median(a) if len(a) else np.NaN


def get_ensembl_ids_from_anc(ancestor, phyltree, ensembl_version=ENSEMBL_VERSION):
    return [SP2GENEID[ensembl_version][sp] for sp in phyltree.species[ancestor]]


def iter_glob_subtree_files(genetreelistfile, ancestor, filesuffix, rootdir='.',
                            subtreesdir='subtreesCleanO2', exclude=None):
    countlines = 0
    countfiles = 0
    for line in genetreelistfile:
        countlines += 1
        genetree = line.rstrip()
        files_pattern = op.join(rootdir, genetree, subtreesdir,
                                  ancestor + genetree + '*' + filesuffix)
        files_reg = re.compile(re.escape(files_pattern).replace('\\*', '(.*)'))
        exclude_reg = exclude if exclude is None else re.compile(exclude)

        logger.info("pattern: '%s' '%s'", files_pattern, files_reg.pattern)
        for subtreefile in glob(files_pattern):
            if exclude is None or not exclude_reg.search(subtreefile):
                countfiles += 1
                subtreesuffix = files_reg.search(subtreefile).group(1)
                subtree = ancestor + genetree + subtreesuffix
                logger.debug('next file: %s', subtreefile)
                yield subtreefile, subtree, genetree

    logger.info('%d lines, %d subtrees' % (countlines, countfiles))


def get_children(tree, node):
    return node.children


def get_al_stats(genetreelistfile, ancestor, phyltreefile, rootdir='.',
                 subtreesdir='subtreesCleanO2', filesuffix='_genes.fa',
                 ensembl_version=ENSEMBL_VERSION, ignore_outgroups=False,
                 ignore_error=True):
    """Gather characteristics of the **input alignments**, and output them as
    a tsv file."""

    #phyltree = PhylTree.PhylogeneticTree(phyltreefile)
    #ensembl_ids_anc = get_ensembl_ids_from_anc(ancestor, phyltree, ensembl_version)
    #pattern = '^(' + '|'.join(ensembl_ids_anc) + ')'
    
    stats_names  = [typ + '_' + measure
                    for typ in ('glob', 'mean', 'med', 'std', 'w_mean', 'w_std')
                    for measure in ('len', 'A','C','G','T', 'GC', 'N',
                                    'gaps', 'CpG')]
    stats_names += ['glob_stops']
    # /!\ WARNING for future self: the below summary stats are **column-wise**!!!
    # (VS sequence-wise above)
    stats_names += ['%s_%s_%s' %(seqtype, measure, typ)
                    for seqtype in ('nucl', 'codon')
                    for measure in ('entropy', 'parsimony')
                    for typ in ('mean', 'median', 'std')]
    if ignore_outgroups:
        ##TO REMOVE
        stats_names = ['ingroup_'+s for s in stats_names] + ['outgroups']

    print('\t'.join(['subtree', 'genetree'] + stats_names))

    for alfile, subtree, genetree in iter_glob_subtree_files(genetreelistfile,
                                                             ancestor,
                                                             filesuffix,
                                                             rootdir,
                                                             subtreesdir):
        try:
            al = AlignIO.read(alfile, format='fasta')
            al = ungap(al)
            subtreefile = alfile.replace(filesuffix, '.nwk')
            tree = ete3.Tree(subtreefile, format=1)

            if ignore_outgroups:
                tree, outgroups = find_ingroup_marked(tree)
                pattern = r'^(' + '|'.join(re.escape(s)
                                           for s in tree.get_leaf_names()) + ')$'

                orig_Nseq = len(al)
                al = ungap(algrep(al, pattern))
                outgroupsize = orig_Nseq - len(al)
                if outgroupsize != 2:
                    logger.error("Removed outgroup of size %d ≠ 2 in %s",
                                 outgroupsize, subtree)

            # Compositional stats
            _, compo_stats = make_al_compo(al)

            # Coding sequence stats
            #seq_stops = []
            #for seq in al:
            #    stops = sum(seq[i:(i+3)] in stop_codons
            #                for i in range(0, al.get_alignment_length(), 3))

            # ~Evolutionary stats (conservation): arrays of column-wise values
            seqlabels = tree.get_leaf_names()
            al = reorder_al(al, seqlabels)

            evo_stats = []
            tree = fuse_single_child_nodes_ete3(tree, copy=False)

            ## By nucleotide column, then by codon.
            for nucl, minlength in [(True,6), (False,66)]:
                _, entropy, alint = get_position_stats(al, nucl=nucl, allow_N=True)

                pars_score = parsimony_score(alint, tree, seqlabels,
                                             minlength=minlength,
                                             get_children=get_children)
                evo_stats.extend((entropy.mean(),
                                  median(entropy),
                                  entropy.std(),
                                  pars_score.mean(),
                                  median(pars_score),
                                  pars_score.std()))

            al_stats = ['%g' % s for stat in compo_stats for s in stat[:-1]]
            al_stats += ['%g' % compo_stats[0][-1]] + ['%g' % s for s in evo_stats]
            if ignore_outgroups:
                al_stats += [','.join(l.name for out in outgroups
                             for l in out.iter_leaves())]

            print('\t'.join([subtree, genetree] + al_stats))
                
            #treefiles_pattern = alfiles_pattern.replace(filesuffix, '.nwk')
        except BaseException as err:
            if ignore_error and not isinstance(err, KeyboardInterrupt):
                logger.exception('At file %s', alfile)
            else:
                err.args = (str(err.args[0]) + '. At file %s' % alfile,) + err.args[1:]
                raise


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
        return get_taxon(node, ancgene2sp, ensembl_version=ensembl_version)
    except ValueError:
        return get_taxon_treebest(node)


# Move to dendro.reconciled
def per_node_events(tree, phyltree, aberrant_dist=10000,
                    get_taxon=any_get_taxon, *args, **kwargs):
    """Check each node and its type (duplication, deletion, speciation).
    """
    aberrant_dists = 0
    rebuilt_topo   = 0  # If the tree was created with myProteinTree.ProteinTree.flattenTree
    # and .rebuildTree, with the option "indicator=True".
    wrong_duptypes = 0  # NotImplemented (where duplication != 0 or 2 or 3)

    leaves = 0

    events = ('dup', 'spe', 'speloss', 'duploss')
    # TreeBest and taxa agree.
    sure_events = {evt: 0 for evt in events}

    only_treebest_events = {evt: 0 for evt in events}

    for node in tree.traverse('preorder'):
        if node.dist >= aberrant_dist:
            aberrant_dists += 1
        if getattr(node, '_r', None):
            rebuilt_topo += 1

        if node.is_leaf():
            leaves += 1
            continue

        treebest_isdup = getattr(node, 'D', None)

        taxon = get_taxon(node, *args, **kwargs)
        expected_speciated_taxa = set(t[0] for t in phyltree.items.get(taxon, []))
        children_taxa = set(get_taxon(ch, *args, **kwargs) for ch in node.children)

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
    return sure_events, only_treebest_events, aberrant_dists, rebuilt_topo


def make_ancgene2sp(ancestor, phyltree):
    return re.compile(r'('
                      + r'|'.join(re.escape(s) for s in
                                  list(phyltree.species[ancestor]) +
                                  sorted(phyltree.getTargetsAnc(ancestor),
                                         key=lambda a:len(a),
                                         reverse=True)).replace(' ', '.')
                      + r')([^a-z].*|)$')


def find_ingroup(tree, ancestor, phyltree, ensembl_version, outgroupsize=2,
                 get_taxon=any_get_taxon, ancgene2sp=None):
    # Minimum number of sequences: 2.
    # Can include ingroup species if the sequences are paralogs.
    leaf2sp = {leafname: convert_gene2species(leafname, ensembl_version)
               for leafname in tree.get_leaf_names()}
    #out_sp = [sp for sp in leaf2sp.values() if sp not in phyltree.species[ancestor]]

    ingroup = tree
    found_outgroups = set()
    outgroup_numbers = 0

    def ingroup_like(child):
        sp_set = set(leaf2sp[l] for l in child.get_leaf_names())
        not_all_out = bool(sp_set & phyltree.species[ancestor])
        not_any_out = not sp_set - phyltree.species[ancestor]
        size_ok = len(child) > outgroupsize-outgroup_numbers
        size_delta = abs(len(child) - outgroupsize-outgroup_numbers)
        # If all groups are smaller, the most likely outgroup is the one whose size matches.

        if ancgene2sp is not None:
            child_sp = get_taxon(child, ancgene2sp, ensembl_version)
            parent_sp = get_taxon(child.up, ancgene2sp, ensembl_version)
            no_missing_ancestor = len(phyltree.dicLinks[parent_sp][child_sp]) <= 2
        else:
            no_missing_ancestor = None
        return (not_all_out,
                size_ok,
                not_any_out,
                size_delta,
                no_missing_ancestor)

    while outgroup_numbers < outgroupsize:
        nodechildren = ingroup.get_children()  # Copy the list.
        #out_children = []
        #for child in nodechildren:
        #    # ANY of the leaves in outgroup species.
        #    if set(leaf2sp[l] for l in child.get_leaf_names()) - phyltree.species[ancestor]:
        #        out_children.append(child)
        #if len(out_children) == len(node.children):
        #    out_children = [child for child in out_children
        #                    if not set() - phyltree.species[ancestor]]
        #if not out_children or (len(node.children) - len(out_children) > 1):
        #    if out_children:
        #        logger.error("Unexpected case: >2 children (multifurcation), "
        #                     "including paralogy+orthology. Trying to get the "
        #                     "next outgroups.")
        #        for child in out_children:
        #            nodechildren.remove(child)
        #    # We must find other leaves that are outgroups (paralogs)
        #    # Let's take the smallest subtree.
        #    sorted_children = sorted(nodechildren, key=len)
        #    if all(len(ch) <= outgroupsize for ch in sorted_children):
        #        logger.warning("Can't select paralog outgroup of an ingroup with size < %d (%s)",
        #                     outgroupsize, list(leaf2sp)[0])
        #    elif all(len(ch) > outgroupsize for ch in sorted_children):
        #        logger.warning("Can't select paralog outgroup of size > %d (%s)",
        #                       outgroupsize, list(leaf2sp)[0])
        #    ingroup = sorted_children.pop()
        #    out_children.extend(sorted_children)

        nodechildren.sort(key=ingroup_like)
        # Pop the ingroup:
        try:
            ingroup = nodechildren.pop()
        except IndexError:
            logger.error('Reached node without children: %s '
                           '(outgroup_numbers = %d).',
                           ingroup.name, outgroup_numbers)
            break
        if all(ingroup_like(out) == ingroup_like(ingroup) for out in nodechildren):
            logger.error("%r->%r: Unable to choose outgroup: %s.",
                         ingroup.up, ingroup, ingroup_like(ingroup))
            # Chooses the first one then
            nodechildren.append(ingroup)
            ingroup = nodechildren.pop(0)

        found_outgroups.update(nodechildren)
        outgroup_numbers = sum(len(outg) for outg in found_outgroups)

    assert outgroup_numbers <= outgroupsize, '%s: size of outgroup is %d > %d'\
                               % (ingroup.name, outgroup_numbers, outgroupsize)

    return ingroup, found_outgroups


def find_ingroup_marked(tree, outgroupsize=2):
    """Find the ingroup as the furthest node with sisters having the `is_outgroup` tag."""
    ingroup = tree
    found_outgroups = set()
    while any(getattr(child, 'is_outgroup', 0) for child in ingroup.children):
        try:
            ingroup, = [child for child in ingroup.children
                        if not getattr(child, 'is_outgroup', 0)]
        except ValueError:
            # Unpacking error: this is either a leaf or all children are unmarked.
            break
        found_outgroups.update(ingroup.get_sisters())

    outgroup_numbers = sum(len(outg) for outg in found_outgroups)
    assert outgroup_numbers <= outgroupsize, '%s: size of outgroup is %d > %d'\
                               % (ingroup.name, outgroup_numbers, outgroupsize)

    return ingroup, found_outgroups


def get_childdist_ete3(tree, nodedist):
    return [(ch, ch.dist) for ch in nodedist[0].children]


def count_zero_combinations(tree, exclusive=False):
    """Count the proportions of branch of length 0:
    - triplet;
    - sister (if exclusive is True: excluding triplet);
    - consecutive (if exclusive is True: excluding triplet);
    """
    consecutive_zeros = 0
    n_branches = 0

    sister_zeros = 0
    triplet_zeros = 0
    n_nodes = 0  # number of strictly internal nodes

    for node in tree.traverse():
        all_child_zeros = False
        if not node.is_leaf():
            all_child_zeros = all(ch.dist==0 for ch in node.children)
            sister_zeros += all_child_zeros

        if not node.is_root():  # If using the ingroup, no node is the root...
            # Consecutive
            n_branches += 1
            if node.dist == 0 and node.up.dist == 0:
                consecutive_zeros += 1

            # Triplet
            if not node.is_leaf():
                n_nodes += 1
                if node.dist == 0 and all_child_zeros:
                    triplet_zeros += 1
                    if exclusive:
                        consecutive_zeros -= 1
                        sister_zeros -= 1
        # Exclusive:
        #is_triplet0 = False
        #if not node.is_leaf():
        #    all_child_zeros = all(ch.dist==0 for ch in node.children)
        #    if not node.is_root():
        #        # Triplet
        #        n_nodes += 1
        #        is_triplet0 = (node.dist == 0 and all_child_zeros):
        #        triplet_zeros += int(is_triplet0)

        #    if not is_triplet_zeros or not exclusive:
        #        sister_zeros += all_child_zeros

        #if not node.is_root():  # If using the ingroup, no node is the root...
        #    # Consecutive
        #    n_branches += 1
        #    if (not is_triplet_zeros or not exclusive):
        #        consecutive_zeros += (node.dist == 0 and node.up.dist == 0)


    return (float(consecutive_zeros) / n_branches,
            float(sister_zeros) / (n_nodes+1),
            float(triplet_zeros) / n_nodes)



def get_tree_stats(genetreelistfile, ancestor, phyltreefile, rootdir='.',
                   subtreesdir='subtreesCleanO2',
                   ensembl_version=ENSEMBL_VERSION,
                   ignore_outgroups=False, extended=False, ignore_error=True):
    """Determine the robustness of the tree, and its clock-likeliness.

    To find the robust trees from the given ancestor only, (excluding the
    outgroup) use `subtreesdir="subtreesClean"`."""

    phyltree = PhylTree.PhylogeneticTree(phyltreefile)
    #ensembl_ids_anc = get_ensembl_ids_from_anc(ancestor, phyltree, ensembl_version)
    print('subtree\tgenetree\troot_location\tleaves_robust\tsingle_child_nodes'
          '\troot2tip_mean\troot2tip_sd'
          + ('\tnodes_robust\tonly_treebest_spe\taberrant_dists\trebuilt_topo\t'
             'bootstrap_min\tbootstrap_mean\tconsecutive_zeros\tsister_zeros\t'
             'triplet_zeros'
             if extended else '')
          + ('\toutgroups' if ignore_outgroups else ''))

    ancgene2sp = make_ancgene2sp(ancestor, phyltree)
    all_ancgene2sp = make_ancgene2sp(phyltree.root, phyltree)
    
    for subtreefile, subtree, genetree in iter_glob_subtree_files(genetreelistfile,
                                                              ancestor,
                                                              '.nwk',
                                                              rootdir,
                                                              subtreesdir,
                                                              exclude='_codeml\.nwk$'):
        try:
            tree = ete3.Tree(subtreefile, format=1)
            root_taxon, _ = split_species_gene(tree.name, ancgene2sp)
            
            # Determine if root_taxon is inside (I) or outside (O) of the clade.
            root_location = 'O' if root_taxon is None \
                            else 'I' if root_taxon != ancestor \
                            else '='

            # What about paralog outgroups?? Should NOT happen if prune2family WITHOUT `--latest`.
            #if root_location == 'O':
            if ignore_outgroups:
                # Use the `is_outgroup` mark. When available, this is the safest.
                ingroupmarked, outgroups = find_ingroup_marked(tree, 2)
                if ingroupmarked == tree:
                    logger.error('find_ingroup_marked:No outgroup found! %s',
                                 tree.name)

                # Double-check:
                ingroup, outgroups = find_ingroup(tree, ancestor, phyltree,
                                        ensembl_version,
                                        2,
                                        any_get_taxon,
                                        ancgene2sp=all_ancgene2sp)
                if ingroupmarked != ingroup:
                    while len(ingroupmarked.children)==1:
                        ingroupmarked, = ingroupmarked.children
                        if ingroupmarked == ingroup:
                            break
                    else:
                        logger.warning(
                                '%s: Found 2 ≠ ingroups with 2 methods: '
                                'find_ingroup -> %r ≠ find_ingroup_mark: %r',
                                subtree, ingroup.name, ingroupmarked.name)
                        # This is ignored.

                tree = ingroupmarked
                
            else:
                root_taxon, _ = split_species_gene(tree.name, all_ancgene2sp)

            #logger.debug('Considered genetree root: %s; root taxon: %s; original: %s'
            #      tree.name, root_taxon or ancestor, ancestor)
            expected_species = phyltree.species[root_taxon or ancestor]
            try:
                leaves_robust, single_child_nodes = simple_robustness_test(tree,
                                                                expected_species,
                                                                ensembl_version)
            except KeyError as err:  # Error while converting genename to species
                err.args += (subtreefile,)
                raise

            root_to_tips = np.array([leafdist for _, leafdist in
                                     iter_distleaves(tree,
                                                     get_childdist_ete3,
                                                     tree.get_tree_root())])

            output = (int(leaves_robust), int(single_child_nodes),
                      mean(root_to_tips), std(root_to_tips))

            if extended:
                sure_events, only_treebest_events, aberrant_dists, rebuilt_topo = \
                        per_node_events(tree,
                                        phyltree,
                                        10000,
                                        any_get_taxon,
                                        ancgene2sp=all_ancgene2sp,
                                        ensembl_version=ensembl_version)
                nodes_robust = not any((sure_events['dup'],
                                        sure_events['speloss'],
                                        sure_events['duploss'],
                                        only_treebest_events['dup'],
                                        only_treebest_events['speloss'],
                                        only_treebest_events['duploss']))
                output += (int(nodes_robust), only_treebest_events['spe'], aberrant_dists, rebuilt_topo)
                # Bootstrap values
                B_values = []
                for node in tree.traverse():
                    if not node.is_leaf() and hasattr(node, 'B'):
                        B_values.append(int(node.B))

                B_values = np.array(B_values)
                output += (B_values.min(), B_values.mean())
                output += count_zero_combinations(tree, exclusive=True)
            if ignore_outgroups:
                output += (','.join(l.name for out in outgroups for l in out.iter_leaves()),)

            print('\t'.join((subtree, genetree, root_location) +
                             tuple(str(x) for x in output)))
        except BaseException as err:
            if err.args:
                err.args = (str(err.args[0]) + ". At file %s" % mlcfile,) + err.args[1:]
            else:
                err.args = ("At file %s" % mlcfile,)
            if ignore_error and not isinstance(err, KeyboardInterrupt):
                logger.exception('Unknown exception')
            else:
                raise


#def tree_autocorr(tree):
#    """Compute the standard deviation of rates per branch."""
#    fork_std = []
#    for parent in tree.traverse():
#        if not parent.is_leaf():
#            # Can't have the rate value here.
#

def get_codeml_stats(genetreelistfile, ancestor, phyltreefile, rootdir='.',
                     subtreesdir='subtreesCleanO2', filesuffix='_m1w04.mlc',
                     ensembl_version=ENSEMBL_VERSION,
                     ignore_outgroups=False, ignore_error=True):
    """Gather characteristics of the **codeml results**, and output them as
    a tsv file."""

    phyltree = PhylTree.PhylogeneticTree(phyltreefile)

    stats_header = ['subtree', 'genetree']
    stats_name   = ['ls', 'ns', 'Nbranches',
                    'NnonsynSites', 'NsynSites', 'kappa',
                    'prop_splitseq', 'codemlMP', 'convergence_warning',
                    'treelen', 'dS_treelen', 'dN_treelen',
                    'brlen_mean', 'brlen_std', 'brlen_med', 'brlen_skew',
                    'brOmega_mean', 'brOmega_std', 'brOmega_med', 'brOmega_skew',
                    'brdS_mean', 'brdS_std', 'brdS_med', 'brdS_skew',
                    'brdN_mean', 'brdN_std', 'brdN_med', 'brdN_skew',
                    'r2t_t_mean', 'r2t_t_std',
                    'r2t_dS_mean', 'r2t_dS_std',
                    'r2t_dN_mean', 'r2t_dN_std',
                    'consecutive_zeros_t',
                    'sister_zeros_t',
                    'triplet_zeros_t',
                    'consecutive_zeros_dS',
                    'sister_zeros_dS',
                    'triplet_zeros_dS',
                    'triplet_zeros_dN',
                    'sister_zeros_dN',
                    'consecutive_zeros_dN',
                    'lnL', 'Niter', 'time used', 'outgroups']

    # TODO: number of dN or dS values of zero
    br_len_reg = re.compile(r': ([0-9]+\.[0-9]+)[,)]')

    print('\t'.join(stats_header + stats_name))

    for mlcfile, subtree, genetree in iter_glob_subtree_files(genetreelistfile,
                                                              ancestor,
                                                              filesuffix,
                                                              rootdir,
                                                              subtreesdir):
        try:
            mlc = parse_mlc(mlcfile)

            Nbr = len(mlc['output']['branches'])
            br_lengths = mlc['output']['branch lengths + parameters'][:Nbr]
            br_omegas = mlc['output']['omega']
            dNdS_rows = list(mlc['output']['dNdS']['rows'].values())
            dS_lengths = [float(x) for x in br_len_reg.findall(mlc['output']['dS tree'])]
            dN_lengths = [float(x) for x in br_len_reg.findall(mlc['output']['dN tree'])]

            # Also load the trees and compute the root2tip mean/std
            t_tree = ete3.Tree(mlc['output']['labelled tree'], name='t_' + subtree)
            dS_tree = ete3.Tree(mlc['output']['dS tree'], name='dS_' + subtree)
            dN_tree = ete3.Tree(mlc['output']['dN tree'], name='dN_' + subtree)

            if ignore_outgroups:
                t_tree, outgroups = find_ingroup(t_tree, ancestor, phyltree, ensembl_version, 2)
                dS_tree, _ = find_ingroup(dS_tree, ancestor, phyltree, ensembl_version, 2)
                dN_tree, _ = find_ingroup(dN_tree, ancestor, phyltree, ensembl_version, 2)

                Nbr = len(t_tree.get_descendants())
                br_lengths = [n.dist for n in t_tree.iter_descendants()]

                w_tree = ete3.Tree(mlc['output']['w tree'], format=1)
                # Reformat.
                for wnode in w_tree.iter_descendants():
                    name, omega = wnode.name.split('#')
                    wnode.name = name.rstrip()
                    wnode.add_feature('omega', float(omega))
                w_tree, _ = find_ingroup(w_tree, ancestor, phyltree, ensembl_version, 2)
                br_omegas = [n.omega for n in w_tree.iter_descendants()]

                dS_lengths = [n.dist for n in dS_tree.iter_descendants()]
                dN_lengths = [n.dist for n in dN_tree.iter_descendants()]
            else:
                outgroups = []

            assert len(dS_lengths) == Nbr

            # Not so relevant if keeping the outgroups.
            t_root_to_tips = np.array([leafdist for _, leafdist in
                                 iter_distleaves(t_tree, get_childdist_ete3)])
            dS_root_to_tips = np.array([leafdist for _, leafdist in
                                 iter_distleaves(dS_tree, get_childdist_ete3)])
            dN_root_to_tips = np.array([leafdist for _, leafdist in
                                 iter_distleaves(dN_tree, get_childdist_ete3)])

            # TO ADD:
            # - Number of NAvalues in Nei & Gojobori table.
            ns = mlc['nsls']['ns']
            prop_nonoverlap = len(list_nonoverlapping_NG(mlc)) * 2./(ns*(ns-1))

            stats_row = [mlc['nsls']['ls'],
                         mlc['nsls']['ns'],
                         Nbr,
                         \
                         dNdS_rows[0][1],  # NnonsynSites
                         dNdS_rows[0][2],  # NsynSites
                         mlc['output']['kappa'],
                         prop_nonoverlap,
                         mlc['output']['numbered topology']['MP score'],
                         int(mlc['output']['warnings']['check convergence']),
                         \
                         sum(br_lengths),  # tree length
                         sum(dS_lengths),      # tree length for dS
                         sum(dN_lengths),      # tree length for dN,
                         \
                         mean(br_lengths),    # brlen_mean
                         std(br_lengths),     # brlen_std
                         median(br_lengths),  # brlen_med
                         skew(br_lengths),       # brlen_skew
                         mean(br_omegas),
                         std(br_omegas),
                         median(br_omegas),
                         skew(br_omegas),
                         mean(dS_lengths),
                         std(dS_lengths),
                         median(dS_lengths),
                         skew(dS_lengths),
                         mean(dN_lengths),
                         std(dN_lengths),
                         median(dN_lengths),
                         skew(dN_lengths),
                         \
                         mean(t_root_to_tips),
                         std(t_root_to_tips),
                         mean(dS_root_to_tips),
                         std(dS_root_to_tips),
                         mean(dN_root_to_tips),
                         std(dN_root_to_tips),
                         \
                         *count_zero_combinations(t_tree),
                         *count_zero_combinations(dS_tree),
                         *count_zero_combinations(dN_tree),
                         \
                         mlc['output']['lnL']['loglik'],
                         mlc['output']['lnL']['ntime']
                         ]

            print('\t'.join([subtree, genetree]
                            + ['%g' % s for s in stats_row]
                            + [mlc['Time used'],
                               ','.join(l.name for out in outgroups
                                               for l in out.iter_leaves())]))
        except BaseException as err:
            if err.args:
                err.args = (str(err.args[0]) + ". At file %s" % mlcfile,) + err.args[1:]
            else:
                err.args = ("At file %s" % mlcfile,)
            if ignore_error and not isinstance(err, KeyboardInterrupt):
                logger.exception('Unknown exception')
            else:
                raise


def get_cleaning_stats(genetreelistfile, ancestor, 
                       rootdir='.', subtreesdir='subtreesCleanO2',
                       filesuffix='_genes.fa', ignore_error=True, **kwargs):
    """Data removed from alignment with Gblocks/Hmmcleaner"""
    if kwargs:
        logger.warning('Ignored kwargs: %s', kwargs)

    stats_header = ['subtree', 'genetree', 'gb_Nblocks', 'gb_percent',
                    'hmmc_nseqs', 'hmmc_propseqs', 'hmmc_max',
                    'hmmc_mean_onlycleaned', 'hmmc_mean']

    print('\t'.join(stats_header))
    regex = re.compile(r'\.'+ filesuffix.split('.')[-1] + r'$')
    hmmc_replacement = '_prot_hmm.log'
    # Ad Hoc FIX:
    if filesuffix == '_fsa.fa':
        regex = re.compile(r'/' + re.escape(subtreesdir) + '/'
                           + re.escape(ancestor)
                           + r'(ENSGT[0-9]+[A-Za-z.])'
                           + re.escape(filesuffix) + r'$')
        hmmc_replacement = r'/{}/realign/{}\1_protfsa_hmm.log'.format(
                                subtreesdir, ancestor)
    for alfile, subtree, genetree in iter_glob_subtree_files(genetreelistfile,
                                                             ancestor,
                                                             filesuffix,
                                                             rootdir,
                                                             subtreesdir):
        try:
            # parse Gblocks output
            gb_logfile = alfile + '-gb.htm'
            if op.exists(gb_logfile):
                gb = parse_gb_html(gb_logfile)
                #TODO: fix those percentages relatively to the ungapped al
                output = [gb['Nblocks'], gb['positions']['percent']]
            else:
                output = [None, None]
            
            # parse hmmc output
            hmmc_logfile = regex.sub(hmmc_replacement, alfile, count=1)
            if op.exists(hmmc_logfile):
                hmmc_ranges = parse_seqranges(hmmc_logfile)

                al = AlignIO.read(alfile, format='fasta')
                seqlabels = [record.name for record in al]
                length, seq_nucls, seq_gaps, seq_N, _ = get_seq_counts(al)

                ## stats by sequence:
                ## - prop cleaned / alignment length
                ## - prop cleaned / nongaps
                ## - prop cleaned / known nucl
                seq_stats = np.zeros((len(seqlabels), 3))

                for i, label in enumerate(seqlabels):
                    Ncleaned = float(sum(end-start
                                         for (start, end) in hmmc_ranges[label]))
                    seq_stats[i, :] = [Ncleaned / length,
                                       Ncleaned / (length - seq_gaps[i]),
                                       Ncleaned / (length - seq_gaps[i] - seq_N[i])]

                cleaned_props = seq_stats[:,0]  # FIXME: seq_stats[:,1]
                cleaned_seqs = (cleaned_props > 0).sum()
                
                output += [cleaned_seqs,
                           float(cleaned_seqs)/len(seqlabels),
                           cleaned_props.max(),
                           cleaned_props[cleaned_props > 0].mean(),
                           cleaned_props.mean()]
            else:
                output += [None]*5
            # was the outgroup cleaned?

            print('\t'.join([subtree, genetree]
                            + ['' if x is None else ('%g' % x)
                                for x in output]))
                
            #treefiles_pattern = alfiles_pattern.replace(filesuffix, '.nwk')
        except BaseException as err:
            if ignore_error and not isinstance(err, KeyboardInterrupt):
                logger.exception('At file %s', alfile)
            else:
                err.args = (str(err.args[0]) + '. At file %s' % alfile,) + err.args[1:]
                raise


def read_beastsummary(filename, convert=None):
    # TODO: convert the infinity symbol to 'Inf'
    summary = {}
    with open(filename) as f:
        colnames = next(f).split()[1:]
        if convert is None:
            for line in f:
                fields = line.split()
                summary[fields[0]] = dict(zip(colnames, fields[1:]))
        else:
            for line in f:
                fields = line.split()
                summary[fields[0]] = {k: convert(x) for k,x in zip(colnames, fields[1:])}
    return summary


def get_beast_stats(genetreelistfile, ancestor, 
                    rootdir='.', subtreesdir='subtreesCleanO2',
                    filesuffix='_beastS-summary.txt', ignore_error=True, **kwargs):
    """Data removed from alignment with Gblocks/Hmmcleaner"""
    if kwargs:
        logger.warning('Ignored kwargs: %s', kwargs)
    stats_header = ['%s_%s' % (stat, stype)
                    for stat in ('posterior', 'likelihood', 'treeL_12', 'treeL_3',
                                 'TreeHeight', 'gammaShape', 'rateAG',
                                 'ucldMean_12', 'ucldMean_3', 'ucldStdev_12', 'ucldStdev_3',
                                 'rate_12_mean', 'rate_12_var', 'rate_3_mean', 'rate_3_var',
                                 'birthRateY')  # 'mrca_age_primates', 'mrca_age_simii'
                    for stype in ('mean', 'stdev', 'med')]  # Mean,std over the iterations.

    select_vars = ('posterior', 'likelihood', 'treeLikelihood.myalignment_1,2',
                   'treeLikelihood.myalignment_3',
                   'TreeHeight', 'gammaShape', 'rateAG',
                   'ucldMean.1,2', 'ucldMean.3', 'ucldStdev.1,2', 'ucldStdev.3',
                   'rate.1,2.mean', 'rate.1,2.variance', 'rate.3.mean',
                   'rate.3.variance', 'birthRateY')
    print('\t'.join(['subtree', 'genetree'] + stats_header))
    for beastsummary, subtree, genetree in iter_glob_subtree_files(genetreelistfile,
                                                             ancestor,
                                                             filesuffix,
                                                             rootdir,
                                                             subtreesdir):
        try:
            #if not op.exists(beastsummary):
            #    output = [None]*len(stats_header)
            #else:
            summary = read_beastsummary(beastsummary)
            output = [summary[var][stype] for var in select_vars
                      for stype in ('mean', 'stddev', 'median')]

            print('\t'.join([subtree, genetree] + output))

        except BaseException as err:
            info = 'At file %s' % beastsummary
            if ignore_error and not isinstance(err, KeyboardInterrupt):
                logger.exception(info)
            else:
                if err.args:
                    err.args = (str(err.args[0]) + '. ' + info,) + err.args[1:]
                else:
                    err.args = (info,)
                raise


if __name__ == '__main__':
    logger.setLevel(logging.INFO)
    try:
        from UItools import colorlog
        colorlog.install()
    except ImportError:
        logging.basicConfig(format=logging.BASIC_FORMAT)

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-d', '--debug', action='store_true',
                        help='Set logging to DEBUG level.')
    
    parent_parser = argparse.ArgumentParser(add_help=False)
    parent_parser.add_argument('genetreelistfile', nargs='?',
                               type=argparse.FileType('r'), default=stdin)
    parent_parser.add_argument('ancestor')
    parent_parser.add_argument('phyltreefile')
    parent_parser.add_argument('-r', '--rootdir', default='.', help='[%(default)s]')
    parent_parser.add_argument('-s', '--subtreesdir', default='subtreesCleanO2',
                               help="[%(default)s]")
    parent_parser.add_argument('-e', '--ensembl-version', type=int,
                               default=ENSEMBL_VERSION, help="[%(default)s]")
    parent_parser.add_argument('-i', '--ignore-outgroups', action='store_true',
                               help='Do not take outgroup species into account'\
                                    ' to compute stats.')
    parent_parser.add_argument('-I', '--no-ignore-error', dest='ignore_error',
                               action='store_false')

    subp = parser.add_subparsers(dest='commands', help='type of statistics to compile')

    codemlstats_parser = subp.add_parser('codeml', parents=[parent_parser], aliases=['co'])
    codemlstats_parser.add_argument('-S', '--filesuffix', default='_m1w04.mlc',
                                    help='file suffix of the globbing pattern [%(default)s]')
    codemlstats_parser.set_defaults(func=make_subparser_func(get_codeml_stats))

    alstats_parser = subp.add_parser('alignment', parents=[parent_parser], aliases=['al'])
    alstats_parser.add_argument('-S', '--filesuffix', default='_genes.fa',
                                help='file suffix of the globbing pattern [%(default)s]')
    alstats_parser.set_defaults(func=make_subparser_func(get_al_stats))

    treestats_parser = subp.add_parser('tree', parents=[parent_parser], aliases=['tr'])
    treestats_parser.add_argument('-E', '--extended', action='store_true',
                                  help='Perform robustness test on nodes (instead of leaves VS root)')
    treestats_parser.set_defaults(func=make_subparser_func(get_tree_stats))
    cleaning_parser = subp.add_parser('cleaning', parents=[parent_parser],
                                      aliases=['cl'])
    cleaning_parser.add_argument('-S', '--filesuffix', default='_genes.fa',
                                 help='file suffix of the globbing pattern [%(default)s]')
    cleaning_parser.set_defaults(func=make_subparser_func(get_cleaning_stats))
    beaststats_parser = subp.add_parser('beast', parents=[parent_parser], aliases=['be'])
    beaststats_parser.add_argument('-S', '--filesuffix', default='_beastS-summary.txt',
                             help='file suffix of the globbing pattern [%(default)s]')
    beaststats_parser.set_defaults(func=make_subparser_func(get_beast_stats))

    args = parser.parse_args()
    if args.debug:
        logger.setLevel(logging.DEBUG)
    delattr(args, 'debug')

    args.func(args)

