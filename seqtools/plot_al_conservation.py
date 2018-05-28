#!/usr/bin/env python3
# -*- coding: utf-8 -*-


"""Display a multiple sequence alignment (nucleotide or codon), and the 
conservation scores along it."""


from sys import stderr, stdin
#from time import clock # to display execution time

import os.path
import re

from itertools import product, combinations

import numpy as np
import matplotlib as mpl
mpl.use('Agg', warn=False)
import matplotlib.pyplot as plt
from Bio import AlignIO, Align, Alphabet
try:
    import argparse_custom as argparse
except ImportError:
    import argparse

from plottools import stackedbar

#try:
#    mpl.style.use('softer')
#except OSError:
# From Seaborn
mpl.rcParams['axes.prop_cycle'] = mpl.cycler('color',
                ['4C72B0', '55A868', 'C44E52', '8172B2', 'CCB974', '64B5CD'])
mpl.rcParams['grid.alpha'] = 0.5
mpl.rcParams['grid.linestyle'] = '--'
#    pass

mpl.rcParams['axes.grid'] = True
mpl.rcParams['axes.grid.axis'] = 'x'

# Change all black to dark grey
grey10 = '#1a1a1a'
for param, paramval in mpl.rcParamsDefault.items():
    if paramval == 'k':
        mpl.rcParams[param] = grey10


ext2fmt = {'.fa':    'fasta',
           '.fasta': 'fasta',
           '.mfa':   'fasta',
           '.phy':   'phylip-relaxed'}


CODONS = [''.join(codon) for codon in product(*['ACGT']*3)]
NACODON = '---'
#NCODONS = 
CODON2INT = {codon:i for i,codon in enumerate([NACODON] + CODONS)}
NUCL2INT  = {'-': 0, 'A': 1, 'C': 2, 'G': 3, 'T': 4}
STOPS = ['TAA', 'TAG', 'TGA']

# Does not work...
# Reading the alignment and issuing align[0][0] still yields a single nucleotide.
codonalphabet = Alphabet.Alphabet()
codonalphabet.letters = [NACODON] + CODONS
codonalphabet.size = 3

def make_unif_codon_dist():
    """Return a uniform codon distance matrix:
    each different nucleotide counts as one difference.
    Stop codons have a NaN distance.
    Invalid codons (containing N nucleotides) as well.
    """
    all_values = [NACODON] + CODONS
    size = len(all_values) + 1  # +1 because of invalid codons e.g 'NNN'
    uniform_codon_dist = 1 - np.diagflat(np.ones(size))
    for i,j in combinations(range(size - 1), 2):
        # sum of nucleotide differences
        cod0, cod1 = all_values[i], all_values[j]
        pair_dist = (np.array(list(cod0)) != np.array(list(cod1))).sum()
        uniform_codon_dist[i, j] = uniform_codon_dist[j, i] = pair_dist
    
    for stop in STOPS:
        stop_indices = all_values.index(stop)
        uniform_codon_dist[:, stop_indices] = np.NaN
        uniform_codon_dist[stop_indices, :] = np.NaN

    # Any extra invalid codon gets a NaN distance (e.g codons containing N's).
    uniform_codon_dist[:, size-1] = np.NaN
    uniform_codon_dist[size-1, :] = np.NaN

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


def al2int(align, nucl=False, allow_N=False):
    """Converts a matrix of categorical values to integers"""
    alarray = al2array(align, nucl)
    converter_dict = NUCL2INT if nucl else CODON2INT
    return category2int(alarray, converter_dict, allow_N)

def category2int(array, converter_dict, allow_N=False):
    """Convert an array of categorical values using a dictionary"""
    if not allow_N:
        ufunc = converter_dict.__getitem__
    else:
        #ufunc = lambda residue: converter_dict.get(residue, np.NaN)
        Ncodon_regex = re.compile('[NATCG]+$')
        Ncodon_int = max(converter_dict.values()) + 1
        def ufunc(residue):
            try:
                return converter_dict[residue]
            except KeyError:
                if Ncodon_regex.match(residue):
                    return Ncodon_int
                else:
                    raise

    return np.vectorize(ufunc)(array)


def count_matrix(vint, minlength=66):
    """column-wise matrix of counts of integers"""
    assert vint.dtype == int
    return np.stack([np.bincount(vint[:,i], minlength=minlength) \
                        for i in range(vint.shape[1])], axis=-1)

def freq_matrix(vint, minlength=66):
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


def freqs2entropy(freqs):
    """Compute the entropy of columns from the frequency matrix.

    Return a vector of entropy of the same length as the number of columns."""

    freqs = freqs.copy()
    freqs[freqs == 0] = 1
    return - (freqs * np.log2(freqs)).sum(axis=0)


def pairs_score(vint, indexpairs, dist_mat=UNIF_CODON_DIST):
    nrows, ncols = vint.shape
    pairs = [np.stack((vint[i,:], vint[j,:])) for i, j in indexpairs]
    pair_scores = np.array([[dist_mat[tuple(p[:,i])] for i in range(ncols)] for p in pairs])
    #print(pair_scores)
    return pair_scores.sum(axis=0)

def sp_score(vint, dist_mat=UNIF_CODON_DIST):
    """Compute the sum of all pairs score, column-wise"""
    return pairs_score(vint, combinations(range(vint.shape[0]), 2), dist_mat)

def part_sp_score(vint, split, dist_mat=UNIF_CODON_DIST):
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
    split_sc = part_sp_score(alint, parts)
    print(split_sc.shape)

    return np.stack((manh_dist, split_sc,))


def al_stats(align, nucl=False, allow_N=False):
    """Compute gap proportion and entropy value for each position of the alignment"""
    al_len = align.get_alignment_length()
    if al_len % 3 and not nucl:
        print("Not a codon alignment!", file=stderr)

    gap = '-'     if nucl else '---'
    npos = al_len if nucl else al_len // 3
    step = 1      if nucl else 3

    alint = al2int(align, nucl, allow_N)
    alint = alint[:, alint.sum(axis=0) > 0]

    #gap_prop = np.zeros(npos)
    #al_entropy = np.full(npos, np.NaN)

    #for i in range(npos):
    #    values = al2list(align[:, (step*i):(step*(i+1))])
    #    value_unique, value_prop = np_proportions(values)
    #    #gap_prop[i]   = value_prop[np.argmax(value_unique == gap)] or 0
    #    al_entropy[i] = np_entropy_subfunc(value_unique, value_prop, gap)

    #gap_prop = (alint == 0).sum(axis=0) / alint.shape[0]
    al_freqs = freq_matrix(alint, minlength=(6 if nucl else 66))
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
    #nvalues = 
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
        axes[4].step(x, dist_array.T, where='post', alpha=0.7)
        axes[4].legend(('pearson_corr', 'manhattan dist'), fontsize='xx-small')
    
    if outfile is None:
        plt.show()
    else:
        fig.savefig(outfile)


PLOTS = ['gap_prop', 'entropy', 'gap_entropy', 'al', 'sp_score']
COMP_PLOTS = ['manhattan', 'pearson', 'split_pairs_score']


class AlignPlotter(object):

    #colorcycle = plt.rcParams['axes.prop_cycle']

    plot_properties = {'al': {'title': 'Global scoring (all sequences)',
                              'ylabel': 'Alignment'},
                       'gap': {'ylabel': 'Proportion of gaps\n(G)'},
                               #'ylim': (-0.05,1.05)},
                       'entropy': {'ylabel': 'Entropy\n(H)'},
                       'gap_entropy': {'ylabel': 'Gap-entropy score:\n(1-H)*(1-G)'},
                       'sp': {'ylabel': 'SP score\n(Sum-of-pair differences)'},
                       'manh': {'title': 'Difference scoring between parts',
                                'ylabel': 'manhattan distance'},
                       'pearson': {'ylabel': "Pearson's correlation coefficient"},
                       'part_sp': {'ylabel': 'sum of pair differences'}}

    default_step = [('step', {'where': 'mid', 'alpha': 0.65})]
    default_bar = [('bar', {'width': 1})]

    plot_funcs = {'al':          [('imshow', {'aspect': 'auto'})],
                  #'al':          [('pcolormesh',   {'edgecolors': 'None'}),
                  #                ('invert_yaxis', {})],
                  #'gap':         default_step, #default_bar, "stackplot"
                  'gap':         [(stackedbar, {'width': 1, 'edgecolor': 'none'})],
                  'entropy':     default_bar,
                  'gap_entropy': default_bar,
                  'sp':          default_bar,
                  'manh':        default_step,
                  'pearson':     default_step,
                  'part_sp':     default_step}


    def __init__(self, alint, seqlabels=None, valid_range=(1,64)):
        """Initialize the Plotter instance from a matrix of integers.
        """
        self.alint = alint
        self.x = np.arange(self.alint.shape[1]) # X values for plotting
        
        self.malint = np.ma.array(alint, mask=(alint==0))
        self.seqlabels = seqlabels
        self.plotlist = ['al']
        self.plotdata = {'al': (self.malint,)}
        #print(np.unique(self.malint))
        #print(np.unique(self.malint - self.malint.max()))
        self.valid_range = valid_range

        #cmap_size = valid_range[1] - valid_range[0] # If len(valid_range) > 1
        # Dark2, tab20b, tab20, Set1
        cmap = plt.get_cmap('tab20b', valid_range[1] - valid_range[0] + 1) #, cmap_size)
        cmap.set_over((0.2, 0.2, 0.2))
        #cmap.set_under('k')
        # a norm is needed to set bounds!!!
        
        # Directly use integer values to get colors in cmap: Same speed but less clear.
        #self.plotdata = {'al': (self.malint - valid_range[0],)} # With NoNorm
        #norm = mpl.colors.NoNorm(0, valid_range[1] - valid_range[0])

        bounds = np.arange(valid_range[0]-0.5, valid_range[1] + 0.5)
        norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
        self.plot_funcs['al'][0][1].update(cmap=cmap, norm=norm)
        self.plot_properties['al']['yticks'] = np.arange(alint.shape[0])
        #self.plot_properties['al']['xlim'] = (0, alint.shape[1]) #pcolormesh
        self.plot_properties['al']['xlim'] = (-0.5, alint.shape[1]-0.5) #imshow

        if seqlabels is not None:
            self.plot_funcs['al'].append(('set_yticklabels',
                                          {'labels': seqlabels,
                                           'fontsize': 'xx-small'}))
        
    @classmethod
    def fromfile(cls, infile, format=None, nucl=False, allow_N=False,
                 slice=None, records=None, recordsfile=None, treefile=None):
        align = AlignIO.read(infile,
                             format=(format or
                                     filename2format(getattr(infile, 'name',
                                                             infile)))
                            )

        al_len = align.get_alignment_length()
        if al_len % 3 and not nucl:
            print("Not a codon alignment!", file=stderr)

        if records or recordsfile or treefile:
            if records:
                records = [int(r) for r in records.split(',')]
            else:
                recnames = [rec.name for rec in align]
                if recordsfile:
                    with open(recordsfile) as recf:
                        records = [recnames.index(line.rstrip()) for line in recf]
                elif treefile:
                    #try: alternative with ete3
                    from Bio.Phylo import read as phyloread
                    tree = phyloread(treefile, 'newick')
                    records = [recnames.index(leaf.name) for leaf in tree.get_terminals()]

            align = Align.MultipleSeqAlignment([align[r] for r in records])


        
        if slice:
            slstart, slend = [int(pos) for pos in slice.split(':')]
            if not nucl:
                slstart *= 3
                slend   *= 3
            if not (0 <= slstart < slend <= al_len):
                raise IndexError("Slice indices (%s,%s) out of bounds." \
                                 % (slstart, slend))
            align = align[:,slstart:slend]
        else:
            slstart, slend = None, None

        seqlabels = [record.name for record in align]
        # Convert matrix of codon/nucleotide strings to matrix of integers.
        alint = al2int(align, nucl, allow_N)
        # remove columns being 100% gaps
        alint = alint[:, alint.sum(axis=0) > 0]

        # Init and setup the class instance
        valid_range = (1, 4) if nucl else (1, 64)
        instance = cls(alint, seqlabels, valid_range)
        instance.fromfile_params = {'infile': infile, 'format': format,
                                    'slice': (slstart, slend),
                                    'records': records}
        instance.align = align
        instance.al_len = al_len
        instance.nucl = nucl
        instance.slstart = slstart
        return instance


    def measure(self):
        """Compute column-wise global measures: residue frequencies, entropy, 
        gap proportion, entropy-gap score, SP-score"""

        max_valid = self.valid_range[1]
        al_freqs = freq_matrix(self.alint, minlength=(max_valid+2))
        gap_prop = al_freqs[0,:]
        invalid_prop = al_freqs[(max_valid+1):, :].sum(axis=0)
        if invalid_prop.any():
            #print(al_freqs[(max_valid+1):, :])
            self.gap_prop = np.stack((gap_prop,
                                      invalid_prop), axis=0)
                                      #gap_prop + invalid_prop), axis=1)
            self.plot_funcs['gap'].append(('legend',
                                           {'labels': ('Gaps only',
                                                       #'Gaps + 
                                                       'Invalid values (N)')}))
        else:
            # change from 1D to 2D. Here, same as array.reshape((array.shape[0], 1))
            self.gap_prop = gap_prop[np.newaxis]

        # Convert zeros to ones so that x * log(x) = 0
        self.freqs = al_freqs[1:,:]
        # scipy.stats.entropy
        self.entropy = freqs2entropy(self.freqs)
        #self.gap_entropy = (1 - self.entropy) * (1 - self.gap_prop[:, 1])
        self.gap_entropy = (1 - self.entropy) * (1 - self.gap_prop.sum(axis=0))
        self.sp_score = sp_score(self.alint)
        
        self.plotlist.extend(('gap', 'entropy', 'gap_entropy', 'sp'))
        self.plotdata.update(gap=(self.x, self.gap_prop,),
                             entropy=(self.x, self.entropy,),
                             gap_entropy=(self.x, self.gap_entropy,),
                             sp=(self.x, self.sp_score,))


    def annotate_parts(self, ax):
        """Draw arrows with annotation to show 'compare_parts' on the
        alignment plot."""
        x0 = ax.get_xlim()[1]
        xtext = x0 * 1.05
        for part_number, part in enumerate(self.parts, start=1):
            ytext = np.mean(part)
            ax.annotate(str(part_number), xy=(x0, 0), xytext=(xtext, ytext))
            for row in part:
                ax.annotate("", xy=(x0, row), xytext=(xtext, ytext),
                            arrowprops={'arrowstyle': '->'})



    def comp_parts(self, compare_parts):
        """Compute column-wise metrics of difference between two parts (groups
        of records) of the dataset:
            - manhattan distance between frequencies
            - pearson coefficient between counts
            - sum of pair differences between any possible pair of distinct parts
        """
        if isinstance(compare_parts, str):
            self.parts = []
            for part in compare_parts.split(';'):
                part = part.rstrip(')').lstrip('(')
                parti = []
                for subpart in part.split(','):
                    try:
                        parti.append(int(subpart))
                    except ValueError:
                        parti.extend(range(*(int(i) for i in subpart.split(':'))))
                self.parts.append(parti)
        else:
            self.parts = compare_parts

        assert len(self.parts) == 2

        self.part_al = [np.stack([self.alint[i,:] for i in part]) \
                        for part in self.parts]

        self.part_freqs = fmat1, fmat2  = [freq_matrix(alpart) for alpart in self.part_al]
        self.part_counts = cmat1, cmat2 = [count_matrix(alpart) for alpart in self.part_al]

        self.part_manh_dist = np.abs(fmat1 - fmat2).sum(axis=0)
        self.part_pearson_c = pearson_coeff(cmat1, cmat2, axis=0)
        self.part_sp_score = part_sp_score(self.alint, self.parts)

        self.plotlist.extend(('manh', 'pearson', 'part_sp'))
        self.plotdata.update(manh=(self.x, self.part_manh_dist,),
                             pearson=(self.x, self.part_pearson_c,),
                             part_sp=(self.x, self.part_sp_score,))

        self.plot_funcs['al'].append((self.annotate_parts, {}))

    def makefig(self, figwidth=16, plotlist=None):

        plotlist = plotlist if plotlist is not None else self.plotlist
        nplots = len(plotlist)
        fig, axes = plt.subplots(nplots, sharex=True,
                                 figsize=(figwidth, 4*nplots), squeeze=False)
        for plot, ax in zip(plotlist, axes[:,0]):
            plotfuncname, plot_kwargs = self.plot_funcs[plot][0]
            try:
                plotfunc = getattr(ax, plotfuncname)
            except TypeError:
                plotfunc = plotfuncname
                # ax argument required:
                plot_kwargs.update(ax=ax)

            plotdata = self.plotdata[plot]

            plotfunc(*plotdata, **plot_kwargs)

            for extrafuncname, extra_kwargs in self.plot_funcs[plot][1:]:
                if isinstance(extrafuncname, str):
                    extrafunc = getattr(ax, extrafuncname)
                else:
                    extrafunc = lambda **kw: extrafuncname(ax, **kw)
                extrafunc(**extra_kwargs)

            plot_prop = self.plot_properties[plot]
            ax.set(**plot_prop)

        if self.slstart:
            # Update the displayed xticklabels if start position > 0
            slstart = self.slstart // 3 if not self.nucl else self.slstart
                
            ax.set_xticklabels([(slstart + int(xt)) for xt in \
                                ax.get_xticks()])

        self.fig, self.axes = fig, axes


    def display(self, outfile=None):
        if outfile is None:
            plt.show(block=True)
            #print(clock())
        else:
            self.fig.savefig(outfile, bbox_inches='tight')



def main_old(infile, outfile=None, format=None, nucl=False, allow_N=False,
             records=None, slice=None, compare_parts=None):
    align = AlignIO.read(infile, format=(format or filename2format(infile.name)))
    if records:
        records = [int(r) for r in ','.split(records)]
        align = Align.MultipleSeqAlignment([align[r] for r in records])
    if slice:
        slstart, slend = [int(pos) for pos in slice.split(':')]
        align = align[:,slstart:slend]

    seqlabels = [record.name for record in align]
    gap_prop, al_entropy, alint = al_stats(align, nucl, allow_N)
    dist_array = comp_parts(alint, compare_parts)
    plot_al_stats(gap_prop, al_entropy, alint, dist_array, seqlabels, outfile)


def main(infile, outfile=None, format=None, nucl=False, allow_N=False,
         records=None, recordsfile=None, treefile=None, slice=None, compare_parts=None,
         compare_only=False, figwidth=16, plotlist=None):
    
    if not outfile:
        plt.switch_backend('TkAgg')

    align_plot = AlignPlotter.fromfile(infile, format, nucl, allow_N, slice,
                                       records, recordsfile, treefile)
    if not compare_only:
        align_plot.measure()
    if compare_parts:
        align_plot.comp_parts(compare_parts)
    if plotlist is not None:
        plotlist = plotlist.split(',')
    align_plot.makefig(figwidth, plotlist)
    align_plot.display(outfile)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('infile', nargs='?', default=stdin,
                        type=argparse.FileType('r'))
    parser.add_argument('-o', '--outfile')
    parser.add_argument('-f', '--format', help='Force format usage.' \
                        ' Can be any format accepted by Bio.alignIO')
    parser.add_argument('-n', '--nucl', action='store_true',
                        help='Process per nucleotide position (instead of codon)')
    parser.add_argument('-N', '--allow-N', action='store_true',
                        help='Tolerate "N" nucleotides as NA value.')
    rec_g = parser.add_mutually_exclusive_group()
    rec_g.add_argument('-r', '--records',
                        help='select records (coma-sep list of 0-based integers)')
    rec_g.add_argument('-R', '--recordsfile',
                        help='select records from a file (one record name per line)')
    rec_g.add_argument('-t', '--treefile',
                        help='select records from a newick tree (reorder the sequences)')
    parser.add_argument('-s', '--slice',
                        help='select positions (start:end). '\
                             '0-based, end excluded, in number of '\
                             'codons/nucleotides.')
    parser.add_argument('-c', '--compare-parts',
                        help='Plot the per-column correlation between two ' \
                             'groups of data, e.g "(0,1);(2,3,4)"')
    parser.add_argument('-C', '--compare-only', action='store_true',
                        help='Do not display global scores')
    parser.add_argument('-w', '--figwidth', default=16, type=float,
                        help='Figure width (inches) [%(default)s]')
    parser.add_argument('-p', '--plotlist', default='al,gap',
                        help='comma-sep list of plots. Valid values are: '\
                             'al,gap,entropy,gap_entropy,sp,manh,pearson,'
                             'part_sp. [%(default)s]')
    
    args = parser.parse_args()

    main(**vars(args))
