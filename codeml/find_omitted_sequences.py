#!/usr/bin/env python3
# -*- coding: utf-8 -*-


from sys import stdin
from Bio import AlignIO
import ete3
import argparse
from UItools.autoCLI import make_subparser_func
from codeml.subtrees_stats import iter_glob_subtree_files
from codeml.codemlparser2 import parse_mlc
from dendro.trimmer import fuse_single_child_nodes_ete3
import logging


def find_omitted_seqs_in_al(genetreelist, ancestor, rootdir='.',
                            subtreesdir='subtreesGoodQualO2'):

    count_ok = 0
    count_bad = 0
    for alfile, subtree, genetree in iter_glob_subtree_files(genetreelist,
                                                             ancestor,
                                                             '_genes.fa',
                                                             rootdir,
                                                             subtreesdir):
        al = AlignIO.read(alfile, format='fasta')
        #al = ungap(al)
        subtreefile = alfile.replace('_genes.fa', '.nwk')
        tree = ete3.Tree(subtreefile, format=1)

        al_seqnames = set([record.name for record in al])
        tree_leaves = set(tree.get_leaf_names())

        al_lacking = al_seqnames - tree_leaves
        tree_lacking = tree_leaves - al_seqnames
        if al_lacking:
            count_bad += 1
            print('%s:AL lacks: %s' % (subtree, ' '.join(al_lacking)))
        elif tree_lacking:
            count_bad += 1
            print('%s:TREE lacks: %s' % (subtree, ' '.join(al_lacking)))
        else:
            count_ok += 1

        if al_lacking and tree_lacking:
            print('AL & TREE mismatched!!!')

    print('%s ok; %s bad' % (count_ok, count_bad))


def find_mismatching_codemltree(genetreelist, ancestor, rootdir='.',
                                subtreesdir='subtreesGoodQualO2'):
    count_ok = 0
    count_bad = 0
    count_missing = 0

    for subtreefile, subtree, genetree in iter_glob_subtree_files(
                                                genetreelist,
                                                ancestor,
                                                '.nwk',
                                                rootdir,
                                                subtreesdir,
                                                exclude='_codeml\.nwk$'):
        original_tree = ete3.Tree(subtreefile, format=1)
        codeml_tree = ete3.Tree(subtreefile.replace('.nwk', '_codeml.nwk'),
                                format=1)
        mlcfile = subtreefile.replace('.nwk', '_m1w04.mlc')
        try:
            mlc = parse_mlc(mlcfile)
        except FileNotFoundError:
            print('Missing file: %s', mlcfile)
            continue
        mlc_tree = ete3.Tree(mlc['output']['labelled tree'])
        
        # Probably unnecessary
        normed_orig_tree = fuse_single_child_nodes_ete3(original_tree)
        
        compare_codeml = normed_orig_tree.compare(codeml_tree)
        compare_mlc = normed_orig_tree.compare(mlc_tree)
        
        #if compare_codeml['ref_edges_in_source'] != 1 \
        #        or compare_codeml['source_edges_in_ref'] != 1:
        #    count_bad += 1
        #    print('%s: ≠ CODEML tree (rf = %g)' % (subtree,
        #                                           compare_codeml['rf']))

        orig_leaves = set(original_tree.get_leaf_names())

        bad = False
        if orig_leaves != set(codeml_tree.get_leaf_names()):
            bad |= True
            print('%s: ≠ CODEML tree (sequences)' % subtree)
        if orig_leaves != set(mlc_tree.get_leaf_names()):
            bad |= True
            print('%s: ≠ MLC tree (sequences)' % subtree)
        if compare_codeml['rf'] and compare_codeml['effective_tree_size']>0:
            bad |= True
            print('%s: ≠ CODEML tree (rf = %g)' % (subtree,
                                                   compare_codeml['rf']))
        if compare_mlc['rf'] and compare_mlc['effective_tree_size']>0:
            bad |= True
            print('%s: ≠ MLC tree (rf = %g)' % (subtree,
                                                compare_mlc['rf']))
        if bad:
            count_bad += 1
        else:
            count_ok += 1

    print('%s ok; %s bad' % (count_ok, count_bad))


if __name__ == '__main__':
    logging.basicConfig()

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-q', '--quiet', action='store_true')

    common_parser = argparse.ArgumentParser(add_help=False)
    common_parser.add_argument('genetreelist', nargs='?',
                        type=argparse.FileType('r'), default=stdin)
    common_parser.add_argument('ancestor')
    common_parser.add_argument('-r', '--rootdir', default='.', help='[%(default)s]')
    common_parser.add_argument('-s', '--subtreesdir', default='subtreesGoodQualO2',
                        help='[%(default)s]')

    subp = parser.add_subparsers(dest='commands', help='Type of check')

    altree_p = subp.add_parser('altree', parents=[common_parser])
    altree_p.set_defaults(func=make_subparser_func(find_omitted_seqs_in_al))

    codemltree_p = subp.add_parser('codemltree', parents=[common_parser])
    codemltree_p.set_defaults(func=make_subparser_func(find_mismatching_codemltree))

    args = parser.parse_args()
    if not args.quiet:
        logging.getLogger('codeml.subtrees_stats').setLevel(logging.INFO)
    delattr(args, 'quiet')

    print(args)
    args.func(args)

