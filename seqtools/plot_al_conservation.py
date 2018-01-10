#!/usr/bin/env python3
# -*- coding: utf-8 -*-


from sys import stderr, stdin

import os.path

from itertools import combinations_with_replacement

import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from Bio import AlignIO, Alphabet
import argparse


ext2fmt = {'.fa':    'fasta',
           '.fasta': 'fasta',
           '.mfa':   'fasta',
           '.phy':   'phylip-relaxed'}

codons = [''.join(codon) for codon in combinations_with_replacement('ACGT', 3)]
codons.append('---')

CodonAlphabet = Alphabet.Alphabet()
Alphabet.letters = codons
Alphabet.size = 1


def filename2format(filename):
    _, ext = os.path.splitext(filename)
    return ext2fmt[ext]


def al2list(align_col):
    """Convert an alignment object to a list of strings"""
    return [str(record.seq) for record in align_col]


def proportions(values):
    """ Compute proportion of each possible value.
    return dictionary with proportion of each value.
    values: list of elements."""
    value_set = set(values)
    L = len(values)
    return {val: float(values.count()) / L}


def np_proportions(values):
    # Convert values to integers
    vunique, vcounts = np.unique(values, return_counts=True)
    return vunique, vcounts.astype(float) / len(values)


def entropy(values, na=None):
    value_prop = proportions(values)
    return -sum(p*np.log2(p) for val, p in value_prop.items() if val != na)


def np_entropy_subfunc(value_unique, value_prop, na=None):
    if na is None:
        na_pos = np.ones(len(value_unique)).astype(bool)
    else:
        na_pos = value_unique == na

    return - (value_prop * np.log2(value_prop))[~na_pos].sum()


def np_entropy(values, na=None):
    value_unique, value_prop = np_proportions(values)
    return np_entropy_subfunc(value_unique, value_prop, na)


def al_stats(align, nucl=False):
    """Compute gap proportion and entropy value for each position of the alignment"""
    al_len = align.get_alignment_length()
    if al_len % 3:
        print("Not a codon alignment!", file=stderr)

    gap = '-'     if nucl else '---'
    npos = al_len if nucl else al_len // 3
    step = 1      if nucl else 3

    gap_prop = np.zeros(npos)
    al_entropy = np.ones(npos) * np.NaN
    is_gap = np.zeros((len(align), npos))

    for i in range(npos):
        values = al2list(align[:, (step*i):(step*(i+1))])
        is_gap[:,i] = np.array(values) == gap
        value_unique, value_prop = np_proportions(values)
        gap_prop[i]   = value_prop[np.argmax(value_unique == gap)] or 0
        al_entropy[i] = np_entropy_subfunc(value_unique, value_prop, gap)

    print(is_gap.shape)
    return gap_prop, al_entropy, is_gap


def plot_al_stats(gap_prop, al_entropy, is_gap, outfile=None):
    """"""
    if outfile is None:
        plt.switch_backend('Qt4Agg')

    fig, axes = plt.subplots(4, sharex=True)

    x = np.arange(len(gap_prop))
    axes[0].bar(x, gap_prop)
    axes[1].bar(x, al_entropy)
    axes[2].bar(x, (1-gap_prop) * (1 - al_entropy))
    axes[3].imshow(is_gap)
    axes[3].autoscale(False)

    axes[0].set_ylabel("Proportion of gaps (G)")
    axes[1].set_ylabel("Entropy (H)")
    axes[2].set_ylabel("Score : (1 - G)*(1 - H)")
    axes[3].set_ylabel("Alignment")
    axes[-1].set_xlabel("Residue position")
    
    if outfile is None:
        plt.show()
    else:
        fig.savefig(outfile)


def main(infile, outfile=None, format=None, nucl=False):
    align = AlignIO.read(infile, format=(format or filename2format(infile.name)))
    gap_prop, al_entropy, is_gap = al_stats(align, nucl)
    plot_al_stats(gap_prop, al_entropy, is_gap, outfile)



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('infile', nargs='?', default=stdin,
                        type=argparse.FileType('r'))
    parser.add_argument('-o', '--outfile')
    parser.add_argument('-f', '--format', help='Force format usage.' \
                        ' Can be any format accepted by Bio.alignIO')
    parser.add_argument('-n', '--nucl', action='store_true',
                        help='Process per nucleotide position (instead of codon)')
    
    
    args = parser.parse_args()
    main(**vars(args))
