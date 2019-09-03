#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
"""Append the species name to a gene name in a sequence alignment.
Can convert ensembl gene names to species, as well as assembly names
(ex: loxAfr3)."""

import argparse
from Bio import SeqIO
from genomicustools.identify import ultimate_seq2sp


ENSEMBL_VERSION = 87


# ~~> genomicustools.identify
def iter_specify(inputfile, labelfmt='{gene}_{sp}', format="fasta",
                 dots=True, ensembl_version=ENSEMBL_VERSION):
#def iter_specify(iterator, name_attr='.id', **)
    #yield element, (newname)

    seen_labels = {}
    for record in SeqIO.parse(inputfile, format):
        sp = ultimate_seq2sp(record.id, ensembl_version)
        if dots:
            sp = sp.replace(' ', '.')
        record.id = labelfmt.format(gene=record.id, sp=sp)

        if record.id in seen_labels:
            record.id += '.%d' % seen_labels[record.id]
            seen_labels[record.id] += 1
        else:
            seen_labels[record.id] = 1

        record.description = ''
        #print(record.id)
        yield record


def tree_specify(tree, labelfmt='{gene}_{sp}', dots=True,
                 ensembl_version=ENSEMBL_VERSION, inplace=True):
    if not inplace:
        tree = tree.copy()
    seen_labels = {}
    for leaf in tree.iter_leaves():
        sp = ultimate_seq2sp(leaf.name, ensembl_version)
        if dots:
            sp = sp.replace(' ', '.')
        leaf.name = labelfmt.format(gene=leaf.name, sp=sp)
        leaf.add_feature('S', sp)

        if leaf.name in seen_labels:
            leaf.name += '.%d' % seen_labels[leaf.name]
            seen_labels[leaf.name] += 1
        else:
            seen_labels[leaf.name] = 1
    if not inplace:
        return tree


def specify(inputfile, outfile, labelfmt='{gene}_{sp}', format="fasta", dots=True,
            ensembl_version=ENSEMBL_VERSION):

    if format == 'nwk':
        import ete3
        tree = ete3.Tree(inputfile, format=1)
        tree_specify(tree, labelfmt, dots, ensembl_version)
        tree.write(outfile=outfile, format=1)
    else:
        SeqIO.write(iter_specify(inputfile, labelfmt, format, dots, ensembl_version),
                    outfile, format)
        #SeqIO.write(iter_specify(SeqIO.parse(inputfile, format),
        #            labelfmt,  dots, ensembl_version)

#def tree_specify(inputfile, outfile, labelfmt=

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('inputfile')
    parser.add_argument('outfile') 
    parser.add_argument('-l', '--labelfmt', default="{gene}_{sp}",
                        help='[%(default)s]')
    parser.add_argument('-f', '--format', default="fasta",
                        help='[%(default)s]')
    parser.add_argument('-e', '--ensembl-version', type=int,
                        default=ENSEMBL_VERSION, help='[%(default)s]')

    args = parser.parse_args()
    specify(**vars(args))
