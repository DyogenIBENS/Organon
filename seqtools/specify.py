#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
"""Append the species name to a gene name in a sequence alignment.
Can convert ensembl gene names to species, as well as assembly names
(ex: loxAfr3)."""

from sys import stdout
import os.path as op
import re
import argparse
from Bio import SeqIO
from genomicustools.identify import ultimate_seq2sp, SP2GENEID


ENSEMBL_VERSION = 87

KEYS = ('gene', 'sp', 'shortsp')  # Sorted by dependency order
DEFAULT_IN_FMT = r'(?P<gene>.*)'  # Use the named group feature.
# Here we capture the whole string as the gene name.
DEFAULT_OUT_FMT = TREEBEST_OUT_FMT = '{gene}_{sp}'
ALE_OUT_FMT = '{sp}_{gene}'

special_short2long = {
        'Scere': 'Saccharomyces cerevisiae',
        'Celeg': 'Caenorhabditis elegans',
        'Dmela': 'Drosophila melanogaster',
        'Csav' : 'Ciona savignyi',
        'Hsa'  : 'Homo sapiens'
        }
special_long2short = {
       'Saccharomyces cerevisiae': 'Scere',
       'Caenorhabditis elegans'  : 'Celeg',
       'Drosophila melanogaster' : 'Dmela',
       'Ciona savignyi'          : 'Csav',
       'Homo sapiens'            : 'Hsa'
        }

def identify_sp(infos, ensembl_version=ENSEMBL_VERSION):
    try:
        return ultimate_seq2sp(infos['gene'], ensembl_version)
    except KeyError:
        return special_short2long.get(infos['shortsp'],
                    convert_gene2sp('ENS'+infos['shortsp'].upper()+'G', ensembl_version))

def identify_shortsp(infos, ensembl_version):
    sp = infos['sp']
    return special_long2short.get(sp,
                                  SP2GENEID[ensembl_version][sp][3:6].capitalize())

identifiers = dict(sp=identify_sp, shortsp=identify_shortsp)


def parse_label(label, input_fmt=DEFAULT_IN_FMT, ensembl_version=ENSEMBL_VERSION):
    # NOTE: this is impossible to ask for an output gene label from a species label.
    infos = re.compile(input_fmt).match(label).groupdict()
    if 'sp' not in infos and 'sp_under' in infos:
        infos['sp'] = infos.pop('sp_under').replace('_', ' ')
    elif 'sp' not in infos and 'sp_dot' in infos:
        infos['sp'] = infos.pop('sp_dot').replace('.', ' ')
    for key in KEYS[1:]:
        if key not in infos:
            infos[key] = identifiers[key](infos, ensembl_version)
    return infos


REGEX_INT = re.compile(r'\d+')


def parse_simple_sed_subst(command: str):
    """Make a substitution function from a sed-like expression.
    The pattern can be any valid python re regex."""
    if command.startswith('s'):
        splitter = command[1]  # Usually '/'
    _, pattern, repl, modifiers = command.split(splitter)
    flags = False
    if 'Ii' in modifiers:
        flags |= re.I
    if 'mM' in modifiers:
        flags | re.M
    regex = re.compile(pattern, flags)
    match_counts = REGEX_INT.search(modifiers)
    count = 1
    if 'g' in modifiers:
        count = 0
    elif match_counts:
        count = int(match_counts.group(0))

    def substituter(string): return regex.sub(repl, string, count)
    return substituter

DEFAULT_TRANSFORMS = dict(sp=lambda string: string.replace(' ', '.'))
# NOTE that CAFE v5 (Hahn et al) does not accept dots in species names...


# ~~> genomicustools.identify
def makelabel(elements, get_label, input_fmt=DEFAULT_IN_FMT, label_fmt=DEFAULT_OUT_FMT,
            transforms=DEFAULT_TRANSFORMS,
            ensembl_version=ENSEMBL_VERSION):
    seen_labels = {}
    for elem in elements:
        infos = parse_label(get_label(elem), input_fmt, ensembl_version)
        for key, transform in transforms.items():
            infos[key] = transform(infos[key])
        label = label_fmt.format(**infos)

        if label in seen_labels:
            seen_labels[label] += 1
            label += '.%d' % seen_labels[label] - 1
        else:
            seen_labels[label] = 1
        yield elem, label, infos


def seq_specify(inputfile, file_fmt="fasta",
                input_fmt=DEFAULT_IN_FMT, label_fmt=DEFAULT_OUT_FMT,
                transforms=DEFAULT_TRANSFORMS,
                ensembl_version=ENSEMBL_VERSION):

    sequences = SeqIO.parse(inputfile, file_fmt)

    for record, label, _ in makelabel(sequences, lambda record: record.id,
                                    input_fmt, label_fmt, transforms,
                                    ensembl_version):
        record.id = label
        record.description = ''
        #print(record.id)
        yield record


def tree_specify(tree, input_fmt=DEFAULT_IN_FMT, label_fmt=DEFAULT_OUT_FMT,
                 transforms=DEFAULT_TRANSFORMS,
                 ensembl_version=ENSEMBL_VERSION, inplace=True):
    if not inplace:
        tree = tree.copy()
    for leaf, label, infos in makelabel(tree.iter_leaves(), lambda leaf: leaf.name,
            input_fmt, label_fmt, transforms, ensembl_version):
        leaf.name = label
        leaf.add_feature('S', infos['sp'])

    if not inplace:
        return tree


EXT2FMT = {'.fa': 'fasta',
           '.phy': 'phylip-relaxed',
           '.nw': 'nwk',
           '.nwk': 'nwk',
           '.newick': 'nwk',
           '.tree': 'nwk'}


def specify(inputfile, outfile, input_fmt=DEFAULT_IN_FMT,
            label_fmt=DEFAULT_OUT_FMT, transform=None, file_fmt=None,
            ensembl_version=ENSEMBL_VERSION):

    if file_fmt is None:
        ext = op.splitext(inputfile)[1]
        file_fmt = EXT2FMT[ext]
    transforms = DEFAULT_TRANSFORMS
    if transform is not None:
        for arg in transform:
            key, val = arg.split('=', 1)
            transforms[key] = parse_simple_sed_subst(val)

    if file_fmt == 'nwk':
        import ete3
        tree = ete3.Tree(inputfile, format=1)
        tree_specify(tree, input_fmt, label_fmt, transforms, ensembl_version)
        if outfile is stdout:
            print(tree.write(format=1))
        else:
            tree.write(outfile=outfile, format=1, format_root_node=True)
    else:
        SeqIO.write(seq_specify(inputfile, file_fmt, input_fmt, label_fmt, transforms,
                                ensembl_version),
                    outfile, file_fmt)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('inputfile')
    parser.add_argument('outfile', nargs='?', default=stdout)
    parser.add_argument('-i', '--input-fmt', default=DEFAULT_IN_FMT,
                        help='[%(default)s]')
    parser.add_argument('-l', '--label-fmt', default=DEFAULT_OUT_FMT,
                        help='[%(default)s]')
    parser.add_argument('-t', '--transform', action='append',
                        help='[%(default)s]')
    parser.add_argument('-f', '--file-fmt',
                        help='[automatic detection of nwk/fasta/phylip]')
    parser.add_argument('-e', '--ensembl-version', type=int,
                        default=ENSEMBL_VERSION, help='[%(default)s]')

    args = parser.parse_args()
    specify(**vars(args))
