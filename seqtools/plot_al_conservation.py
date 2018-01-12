#!/usr/bin/env python3
# -*- coding: utf-8 -*-


from sys import stderr, stdin

import os.path

from itertools import product

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

CODONS = [''.join(codon) for codon in product(*['ACGT']*3)]
NACODON = '---'
CODON2INT = {codon:i for i,codon in enumerate([NACODON] + CODONS)}
NUCL2INT  = {'-': 0, 'A': 1, 'C': 2, 'G': 3, 'T': 4}

# Does not work...
# Reading the alignment and issuing align[0][0] still yields a single nucleotide.
codonalphabet = Alphabet.Alphabet()
codonalphabet.letters = [NACODON] + CODONS
codonalphabet.size = 3


def al2array(align, nucl=False):
    """"""
    nseq = len(align)
    al_len = align.get_alignment_length()
    npos = al_len if nucl else al_len // 3
    step = 1      if nucl else 3
    #alarray = np.array([['']*npos]*len(align))
    
    return np.array([[seq[(j*step):((j+1)*step)] for j in range(npos)]
                                                 for seq in al2list(align)])


def al2int(align, nucl=False):
    """Converts a matrix of categorical values to integers"""
    alarray = al2array(align, nucl)
    converter_dict = NUCL2INT if nucl else CODON2INT
    return category2int(alarray, converter_dict)


def category2int(array, converter_dict):
    """Convert an array of categorical values using a dictionary"""
    return np.vectorize(converter_dict.__getitem__)(array)

def freq_matrix(vint, minlength=65):
    """Convert matrix of integers to frequencies (per column)"""
    counts = np.stack([np.bincount(vint[:,i], minlength=minlength) \
                        for i in range(vint.shape[1])], axis=-1)
    return counts.astype(float) / np.alen(vint)

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
    if al_len % 3 and not nucl:
        print("Not a codon alignment!", file=stderr)

    gap = '-'     if nucl else '---'
    npos = al_len if nucl else al_len // 3
    step = 1      if nucl else 3

    alint = al2int(align, nucl)

    #gap_prop = np.zeros(npos)
    #al_entropy = np.ones(npos) * np.NaN

    #for i in range(npos):
    #    values = al2list(align[:, (step*i):(step*(i+1))])
    #    value_unique, value_prop = np_proportions(values)
    #    #gap_prop[i]   = value_prop[np.argmax(value_unique == gap)] or 0
    #    al_entropy[i] = np_entropy_subfunc(value_unique, value_prop, gap)

    #gap_prop = (alint == 0).sum(axis=0) / alint.shape[0]
    al_freqs = freq_matrix(alint, minlength=(5 if nucl else 65))
    gap_prop = al_freqs[0,:]
    # Convert zeros to ones so that x * log(x) = 0
    al_freqs = al_freqs[1:,:]
    al_freqs[al_freqs == 0] = 1
    al_entropy = - (al_freqs * np.log2(al_freqs)).sum(axis=0)

    print(alint.shape)
    return gap_prop, al_entropy, alint


def plot_al_stats(gap_prop, al_entropy, alint, seqlabels=None, outfile=None):
    """"""
    if outfile is None:
        #try:
        #    plt.switch_backend('Qt5Agg')
        #except ImportError:
            plt.switch_backend('TkAgg')

    #try:
    #    alcmap = plt.get_cmap('tab20', alint.max() - 1)
    #except ValueError:

    nvalues = alint.max()
    print(nvalues)
    alcmap = plt.get_cmap('Dark2', nvalues)

    masked_al = np.ma.array(alint, mask=(alint==0))

    fig, axes = plt.subplots(4, sharex=True)

    x = np.arange(len(gap_prop))
    axes[0].step(x, gap_prop)
    axes[1].bar(x, al_entropy, width=1)
    axes[2].bar(x, (1-gap_prop) * (1 - al_entropy), width=1)
    
    #axes[3].imshow(is_gap, cmap='binary_r', aspect='auto') #, interpolation='gaussian')
    axes[3].imshow(masked_al, cmap=alcmap, aspect='auto') #, interpolation='gaussian')

    axes[0].set_ylabel("Proportion of gaps (G)")
    axes[1].set_ylabel("Entropy (H)")
    axes[2].set_ylabel("Score : (1 - G)*(1 - H)")
    axes[3].set_ylabel("Alignment")
    if seqlabels is not None:
        axes[3].set_yticks(np.arange(alint.shape[0]))
        axes[3].set_yticklabels(seqlabels, fontsize='xx-small')
    axes[-1].set_xlabel("Residue position")
    
    if outfile is None:
        plt.show()
    else:
        fig.savefig(outfile)


def main(infile, outfile=None, format=None, nucl=False):
    align = AlignIO.read(infile, format=(format or filename2format(infile.name)))
    seqlabels = [record.name for record in align]
    gap_prop, al_entropy, alint = al_stats(align, nucl)
    plot_al_stats(gap_prop, al_entropy, alint, seqlabels, outfile)



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
