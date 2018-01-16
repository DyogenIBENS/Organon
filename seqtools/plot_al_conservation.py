#!/usr/bin/env python3
# -*- coding: utf-8 -*-


from sys import stderr, stdin

import os.path

from itertools import product, combinations

import numpy as np
import matplotlib as mpl
mpl.use('Agg', warn=False)
import matplotlib.pyplot as plt
from Bio import AlignIO, Align, Alphabet
import argparse


ext2fmt = {'.fa':    'fasta',
           '.fasta': 'fasta',
           '.mfa':   'fasta',
           '.phy':   'phylip-relaxed'}

CODONS = [''.join(codon) for codon in product(*['ACGT']*3)]
NACODON = '---'
CODON2INT = {codon:i for i,codon in enumerate([NACODON] + CODONS)}
NUCL2INT  = {'-': 0, 'A': 1, 'C': 2, 'G': 3, 'T': 4}
STOPS = ['TAA', 'TAG', 'TGA']

# Does not work...
# Reading the alignment and issuing align[0][0] still yields a single nucleotide.
codonalphabet = Alphabet.Alphabet()
codonalphabet.letters = [NACODON] + CODONS
codonalphabet.size = 3

def make_unif_codon_dist():
    all_values = [NACODON] + CODONS
    uniform_codon_dist = 1 - np.diagflat(np.ones(len(all_values)))
    for i,j in combinations(range(len(all_values)), 2):
        # sum of nucleotide differences
        cod0, cod1 = all_values[i], all_values[j]
        pair_dist = (np.array(list(cod0)) != np.array(list(cod1))).sum()
        uniform_codon_dist[i, j] = uniform_codon_dist[j, i] = pair_dist
    
    for stop in STOPS:
        uniform_codon_dist[:, all_values.index(stop)] = np.NaN
        uniform_codon_dist[all_values.index(stop), :] = np.NaN

    return uniform_codon_dist

UNIF_CODON_DIST = make_unif_codon_dist()

#def read_codon_matrix(filename):
#    # ~/install/jprime-extra/arvecodon.txt


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


def count_matrix(vint, minlength=65):
    """column-wise matrix of counts of integers"""
    assert vint.dtype == int
    return np.stack([np.bincount(vint[:,i], minlength=minlength) \
                        for i in range(vint.shape[1])], axis=-1)

def freq_matrix(vint, minlength=65):
    """Convert matrix of integers to frequencies (per column)"""
    return count_matrix(vint, minlength).astype(float) / np.alen(vint)

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

# Not used in `main` anymore
def np_entropy_subfunc(value_unique, value_prop, na=None):
    if na is None:
        na_pos = np.full(len(value_unique), True)
    else:
        na_pos = value_unique == na

    return - (value_prop * np.log2(value_prop))[~na_pos].sum()

# Not used in `main`
def np_entropy(values, na=None):
    value_unique, value_prop = np_proportions(values)
    return np_entropy_subfunc(value_unique, value_prop, na)


def pairs_score(vint, indexpairs, dist_mat=UNIF_CODON_DIST):
    nrows, ncols = vint.shape
    pairs = [np.stack((vint[i,:], vint[j,:])) for i, j in indexpairs]
    pair_scores = np.array([[dist_mat[tuple(p[:,i])] for i in range(ncols)] for p in pairs])
    print(pair_scores)
    return pair_scores.sum(axis=0)

def sp_score(vint, dist_mat=UNIF_CODON_DIST):
    """Compute the sum of all pairs score, column-wise"""
    return pairs_scores(vint, combinations(range(nrows), 2), dist_mat)

def split_score(vint, split, dist_mat=UNIF_CODON_DIST):
    assert len(split) == 2
    return pairs_score(vint, product(*split), dist_mat)

def cov(X,Y, axis=0):
    mx = X.mean(axis=axis)
    my = Y.mean(axis=axis)
    return ((X - mx) * (Y - my)).mean(axis=axis)


def pearson_coeff(X, Y, axis=0):
    """Column-wise pearson coefficient"""
    coeff = cov(X, Y, axis)
    coeff /= (X.std(axis=axis) * Y.std(axis=axis))
    return coeff


def comp_parts(alint, compare_parts=None):
    if compare_parts is None:
        return None

    parts = [[int(i) for i in part.rstrip(')').lstrip('(').split(',')]
             for part in compare_parts.split(';')]

    assert len(parts) == 2

    alparts = [np.stack([alint[i,:] for i in part]) for part in parts]

    freq_mat1, freq_mat2 = [freq_matrix(alpart) for alpart in alparts]
    count_mat1, count_mat2 = [count_matrix(alpart) for alpart in alparts]

    manh_dist = np.abs(freq_mat1 - freq_mat2).sum(axis=0)
    pearson_c = pearson_coeff(count_mat1, count_mat2, axis=0)
    split_sc = split_score(alint, parts)
    print(split_sc.shape)

    return np.stack((manh_dist, split_sc,))


def al_stats(align, nucl=False):
    """Compute gap proportion and entropy value for each position of the alignment"""
    al_len = align.get_alignment_length()
    if al_len % 3 and not nucl:
        print("Not a codon alignment!", file=stderr)

    gap = '-'     if nucl else '---'
    npos = al_len if nucl else al_len // 3
    step = 1      if nucl else 3

    alint = al2int(align, nucl)
    alint = alint[:, alint.sum(axis=0) > 0]

    #gap_prop = np.zeros(npos)
    #al_entropy = np.full(npos, np.NaN)

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
    # scipy.stats.entropy
    al_freqs[al_freqs == 0] = 1
    al_entropy = - (al_freqs * np.log2(al_freqs)).sum(axis=0)

    print(alint.shape)
    return gap_prop, al_entropy, alint


def plot_al_stats(gap_prop, al_entropy, alint, dist_array=None, seqlabels=None,
                  outfile=None):
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

    nplots = 4 if dist_array is None else 5
    fig, axes = plt.subplots(nplots, sharex=True, figsize=(15,10))

    x = np.arange(len(gap_prop))
    axes[0].step(x, gap_prop, where='post')
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

    if dist_array is not None:
        ax5 = axes[4].twinx()
        l4 = axes[4].step(x, dist_array[0,:], color='#1f77b4', where='post',
                     alpha=0.7, label='manhattan dist')
        l5 = ax5.step(x, dist_array[1,:], color='#ff7f0e', where='post',
                     alpha=0.7, label='pearson corr')
        axes[4].legend(fontsize='xx-small')
        ax5.legend(fontsize='xx-small')
    
    if outfile is None:
        plt.show()
    else:
        fig.savefig(outfile)


def main(infile, outfile=None, format=None, nucl=False, records=None,
         slice=None, compare_parts=None):
    align = AlignIO.read(infile, format=(format or filename2format(infile.name)))
    if records:
        records = [int(r) for r in ','.split(records)]
        align = Align.MultipleSeqAlignment([align[r] for r in records])
    if slice:
        slstart, slend = [int(pos) for pos in slice.split(':')]
        align = align[:,sltart:slend]

    seqlabels = [record.name for record in align]
    gap_prop, al_entropy, alint = al_stats(align, nucl)
    dist_array = comp_parts(alint, compare_parts)
    #print(dist_array)
    plot_al_stats(gap_prop, al_entropy, alint, dist_array, seqlabels, outfile)



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('infile', nargs='?', default=stdin,
                        type=argparse.FileType('r'))
    parser.add_argument('-o', '--outfile')
    parser.add_argument('-f', '--format', help='Force format usage.' \
                        ' Can be any format accepted by Bio.alignIO')
    parser.add_argument('-n', '--nucl', action='store_true',
                        help='Process per nucleotide position (instead of codon)')
    parser.add_argument('-r', '--records',
                        help='select records (coma-sep list of 0-based integers)')
    parser.add_argument('-s', '--slice',
                        help='select positions (start:end). 0-based, end excluded')
    parser.add_argument('-c', '--compare-parts',
                        help='Plot the per-column correlation between two ' \
                             'groups of data, e.g "(0,1);(2,3,4)"')
    
    
    args = parser.parse_args()
    main(**vars(args))
