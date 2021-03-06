#!/usr/bin/env python3
# -*- coding: utf-8 -*-


"""Display a multiple sequence alignment (nucleotide or codon), and the 
conservation scores along it."""


from sys import stdin
#from time import clock # to display execution time

import os.path
import re

from itertools import product, combinations

import numpy as np
import matplotlib as mpl
try:
    mpl.use('Agg', warn=False)
except TypeError:
    # keyword warn disappeared in recent matplotlib
    mpl.use('Agg')
import matplotlib.pyplot as plt
from Bio import AlignIO, Align, Alphabet
try:
    import argparse_custom as argparse
except ImportError:
    import argparse

from functools import reduce

from seqtools.IUPAC import gaps, unknown, nucleotides, ambiguous
from datasci.graphs import stackedbar, plottree
from dendro.bates import rev_dfw_descendants

import logging
logger = logging.getLogger(__name__)

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


# ~~> arrayal
STOPS = ['TAA', 'TAG', 'TGA']
CODONS = [''.join(codon) for codon in product(*[nucleotides]*3)
          if ''.join(codon) not in STOPS]
NACODON = '---'
#NCODONS = 
CODON2INT = {codon:i for i,codon in enumerate([NACODON] + CODONS + STOPS)}
#NUCL2INT  = {'-': 0, 'A': 1, 'C': 2, 'G': 3, 'T': 4}
NUCL2INT = {symbol: i for i, symbol in enumerate(gaps[0] + nucleotides)}

# Does not work...
# Reading the alignment and issuing align[0][0] still yields a single nucleotide.
codonalphabet = Alphabet.Alphabet()
codonalphabet.letters = [NACODON] + CODONS + STOPS
codonalphabet.size = 3


# ~~> arrayal
def make_unif_codon_dist():
    """Return a uniform codon distance matrix:
    each different nucleotide counts as one difference.
    Stop codons have a NaN distance.
    Invalid codons (containing N nucleotides) as well.
    """
    all_values = [NACODON] + CODONS
    size = len(all_values) + len(STOPS) + 1  # +1 because of invalid codons e.g 'NNN'
    uniform_codon_dist = 1 - np.diagflat(np.ones(size))
    for i,j in combinations(range(size - len(STOPS) - 1), 2):
        # sum of nucleotide differences
        cod0, cod1 = all_values[i], all_values[j]
        pair_dist = (np.array(list(cod0)) != np.array(list(cod1))).sum()
        uniform_codon_dist[i, j] = uniform_codon_dist[j, i] = pair_dist
    
    # Stops and invalid codons (e.g codons containing N's) get a NaN distance.
    uniform_codon_dist[:, -4:] = np.NaN
    uniform_codon_dist[-4:, :] = np.NaN

    return uniform_codon_dist

UNIF_CODON_DIST = make_unif_codon_dist()

#def read_codon_matrix(filename):
#    # ~/install/jprime-extra/arvecodon.txt


# ~~> arrayal
def al2list(align_col):
    """Convert an alignment object to a list of strings"""
    return [str(record.seq) for record in align_col]

def iter_strseq(align):
    """Iterate over an alignment object: yield sequences as strings"""
    for record in align:
        yield str(record.seq)


# ~~> arrayal
def al2array(align, nucl=False):
    """"""
    nseq = len(align)
    al_len = align.get_alignment_length()
    npos = al_len if nucl else al_len // 3
    step = 1      if nucl else 3
    #alarray = np.array([['']*npos]*len(align))
    
    return np.array([[seq[(j*step):((j+1)*step)] for j in range(npos)]
                                                 for seq in iter_strseq(align)])


# ~~> arrayal
def al2int(align, nucl=False, allow_N=False):
    """Converts a matrix of categorical values to integers"""
    alarray = al2array(align, nucl)
    converter_dict = NUCL2INT if nucl else CODON2INT
    return category2int(alarray, converter_dict, allow_N)


# ~~> arrayal
def category2int(array, converter_dict, allow_N=False):
    """Convert an array of categorical values using a dictionary"""
    if not allow_N:
        ufunc = converter_dict.__getitem__
    else:
        #ufunc = lambda residue: converter_dict.get(residue, np.NaN)
        Ncodon_regex = re.compile(r'[' + unknown + nucleotides + ''.join(ambiguous) + ']+$')
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


# ~~> arrayal
def count_matrix(vint, minlength=66):
    """column-wise matrix of counts of integers"""
    #assert np.issubdtype(vint.dtype, int  # 
    assert vint.dtype in (int, np.int64, np.int32, np.int16)
    return np.stack([np.bincount(vint[:,i], minlength=minlength) \
                        for i in range(vint.shape[-1])], axis=-1)

#count_matrix = np.vectorize(lambda array_1D: np.bincount(array_1D, minlength=minlength)
# Nope because this processes columns

# ~~> arrayal
def freq_matrix(vint, minlength=66):
    """Convert matrix of integers to frequencies (per column)"""
    return count_matrix(vint, minlength).astype(float) / np.alen(vint)

# ~~> arrayal
def presence_matrix(vint, minlength=66):
    """Convert a sequence of integers into a boolean matrix of presence of each value.
    Zero encodes gaps. >minlength encodes ambiguous/invalid residues."""

    assert len(vint.shape) == 2 and vint.shape[0] == 1, \
            "Can't input more than one sequence " + str(vint.shape)
    count_mat = count_matrix(vint, minlength)
    # Convert gaps (integer zero) and 'N' (integer minlength-1) to unknown: True everywhere.
    presence = count_mat[1:(minlength-1), :].astype(bool)
    # Broadcast row of gaps on all rows of `presence`.
    presence |= count_mat[0, :].astype(bool)
    # Broadcast row of 'N' on all rows of `presence`.
    presence |= count_mat[(minlength-1), :].astype(bool)

    return presence

    ### WTF
def presence_matrice(vint, minlength=66):
    assert len(vint.shape) == 2 and vint.shape[0] == 1, \
            "Can't input more than one sequence " + str(vint.shape)
    count_mat = count_matrix(vint, minlength)
    presence = count_mat.astype(bool)
    # Convert gaps (integer zero) and 'N' (integer minlength-1) to unknown: True everywhere.
    presence[0, :] = True
    return presence

def filename2format(filename):
    _, ext = os.path.splitext(filename)
    return ext2fmt[ext]


def proportions(values):
    """DEPRECATED.
    Compute proportion of each possible value.
    return dictionary with proportion of each value.
    values: list of elements."""
    value_set = set(values)
    L = len(values)
    return {val: float(values.count()) / L}

def entropy(values, na=None):
    """DEPRECATED"""
    value_prop = proportions(values)
    return -sum(p*np.log2(p) for val, p in value_prop.items() if val != na)


# ~~> arrayal/alvector
def np_proportions(values):
    # Convert values to integers
    vunique, vcounts = np.unique(values, return_counts=True)
    return vunique, vcounts.astype(float) / len(values)


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

    Return a vector of entropy of the same length as the number of columns.
    """
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

# ~~> datasci.stats?
def cov(X,Y, axis=0):
    mx = X.mean(axis=axis)
    my = Y.mean(axis=axis)
    return ((X - mx) * (Y - my)).mean(axis=axis)


# ~~> datasci.stats?
def pearson_coeff(X, Y, axis=0):
    """Column-wise pearson coefficient"""
    coeff = cov(X, Y, axis)
    coeff /= (X.std(axis=axis) * Y.std(axis=axis))
    return coeff


def comp_parts(alint, compare_parts=None):
    """DEPRECATED"""
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


# ~~> arrayal/evol
class Parsimony(object):
    #def parsimony_score(alint, tree, seqlabels, minlength=66, get_children=None,
    #                parts=None):
    """
    Computes the column-wise parsimony score based on the provided tree.
    
    parts: [tuples of indices] e.g. [(0,1,2)] represents leaves 0,1,2 as a clade.
    The outgroup must NOT be added into parts.
    """
    def __init__(self, alint, tree, seqlabels, minlength=66, get_children=None,
                 parts=None):

        self.alint = alint
        self.tree = tree
        self.seqlabels = seqlabels
        assert alint.shape[0] == len(seqlabels)
        self.minlength = minlength
        # Bio.Phylo by default.
        self.get_children = (lambda tree, node: node.clades) if get_children is None else get_children
        try:
            self.root = tree.clade
        except AttributeError:
            self.root = tree
        self.iter_tree = list(rev_dfw_descendants(tree, self.get_children, include_leaves=True,
                                        queue=[self.root]))

        # Holds the currently seen nodes, and their parsimony score and sequence.
        self.process_sequences = {}
        self.anc_states = {}  # These could be stored as the same object.

        # column-wise parsimony score
        self.score = np.zeros(alint.shape[1], dtype=int)
        
        self.parts = [set(p) for p in parts] if parts is not None else []
        self.node_to_part = {leaf_n: i for i,p in enumerate(self.parts)
                        for leaf_n in p}
        self.part_scores = [self.score.copy() for p in self.parts]
        self.part_branch_nbs = [0] * len(self.parts)

    def rootwards(self):
        # Index of encountered leaves/sequences
        leaf_nb = self.alint.shape[0] - 1

        process_sequences = self.process_sequences
        parts = self.parts
        merged_parts = [dict() for _ in parts]  # {node: set(seen_part_leaves)}
        # When merged_parts[i][node] == parts[i]: we found the MRCA of part.
        branch_nb = 0
        for parent, children in self.iter_tree:
            logger.debug('* %r -> %s', parent, children)
            if len(children) > 2:
                logger.warning("%d > 2 children at %r", len(children), parent)
            if not children:
                # Parent is a leaf. Obtain the sequence.
                assert parent.name == self.seqlabels[leaf_nb], \
                    "The alignment is not ordered as the tree. Seq %d: %s != leaf %s" \
                        % (leaf_nb, self.seqlabels[leaf_nb], parent.name)
                process_sequences[parent] = presence_matrix(self.alint[np.newaxis,leaf_nb,:], self.minlength)
                try:
                    p = self.node_to_part[parent] = self.node_to_part[leaf_nb]
                    merged_parts[p][parent] = set((leaf_nb,))
                    logger.info('Assigned leaf %d %r to part %d', leaf_nb, parent, p)
                except KeyError:  # When there is no part for this leaf.
                    pass

                leaf_nb -= 1
            else:
                # .pop(ch) ?
                try:
                    children_seqs = [process_sequences[ch] for ch in children]
                except KeyError as err:
                    #logger.debug('Processed sequences:\n%s',
                    #             '\n'.join('%r: %s' % (node, pseq.astype(np.int32))
                    #                       for node, pseq in process_sequences.items()))
                    logger.error('parent = %r; leaf_nb = %s', parent, leaf_nb)
                    raise
                children_inter = reduce(np.logical_and, children_seqs)
                children_union = reduce(np.logical_or, children_seqs)
                #print(children_inter, children_union, sep='\n')
                # Add one to each column where a substitution is needed
                empty_inter = ~children_inter.any(axis=0)

                # The new nucleotide set is the intersection if it's not empty,
                # otherwise the union
                process_sequences[parent] = children_inter
                process_sequences[parent][:, empty_inter] = children_union[:, empty_inter]
                #process_sequences[parent] = np.where(empty_inter,
                #                                     children_union,
                #                                     children_inter)

                if parts:
                    children_parts = list(set((self.node_to_part.get(ch) for ch in children)))
                    if len(children_parts) == 1:
                        # Still inside the same clade. This must be exclusive (monophyletic).
                        #assert all descendant leaves are in the same part.
                        p = children_parts[0]
                        try:
                            merged_parts[p][parent] = set.union(*(merged_parts[p][ch] for ch in children))
                            self.node_to_part[parent] = p
                            logger.info('Assigned node %r to part %d', parent, p)
                            # Just update the intra-clade score.
                            self.part_scores[p] += empty_inter
                            self.part_branch_nbs[p] += len(children)
                            # Still inside one clade, so skip updating the *global* score.
                            continue
                        except TypeError:
                            # p is None because no child was found in node_to_part
                            pass  # So we update the outgroup score...
                    else:
                        # We leave one or more clades.
                        # 1. assert that at least one is monophyletic.
                        part_oldest_node = []  # Monophyly check variable
                        # 2. TODO: Map substitutions to the descendent clades
                        for ch in children:
                            try:
                                p = self.node_to_part[ch]
                            except KeyError:
                                # this is the outgroup, skip. The global score should be updated though
                                continue
                            if merged_parts[p][ch] == parts[p]:
                                # The child is the exclusive MRCA (It contains exactly the part species).
                                # We update the score of the clade:
                                # NOTE: this can't be done on the first postorder traversal...
                                # Because we need to know the parent state to orient the substitution!

                                # Uncomment this to add the stem substitutions:
                                #part_scores[p] += empty_inter
                                #part_branch_nbs[p] += 1  # Add the stem branch
                                part_oldest_node.append(ch)

                        if not part_oldest_node:
                            raise ValueError("The given parts do not make monophyletic clades."
                                             + ';'.join('(%s)' % (
                                                     ','.join(self.seqlabels[i] for i in p))
                                                     for p in parts))
                        # Here, the current score will be updated (see just below)

                branch_nb += len(children)  # Except if children & merged_parts
                self.score += empty_inter
            #print(process_sequences[parent])

        # Number of branches of the **unrooted** tree:
        branch_nb -= 1

        if parts:
            return (self.score.astype(float), *self.part_scores), (branch_nb, *self.part_branch_nbs)

        self.branch_nb = branch_nb
        return self.score.astype(float) / branch_nb
    
    def leafwards(self):
        # Retraverse leafwards to get the ancestral states. Must have traversed rootwards first.
        self.score_leafward = np.zeros(self.alint.shape[1], dtype=int)
        anc_states = self.anc_states
        anc_states[self.root] = self.process_sequences.pop(self.root)
        for parent, children in reversed(self.iter_tree):
            parent_state = anc_states[parent]
            for ch in children:
                child_possibles = self.process_sequences.pop(ch)
                branch_inter = np.logical_and(parent_state, child_possibles)  # No substitution needed.
                branch_union = np.logical_or(parent_state, child_possibles)  # choice needed
                non_empty_inter = branch_inter.any(axis=0)
                self.score_leafward += (~non_empty_inter).sum()  # How many sites with empty intersection (i.e. change)

                # If non empty intersection, the child nucleotides are intersected with those of the parent
                # otherwise, unchanged.
                anc_states[ch] = child_possibles
                anc_states[ch][:, non_empty_inter] = branch_inter[:, non_empty_inter]
        #TODO: change self.part_scores
        return self.score_leafward / self.branch_nb

    def __call__(self):
        self.rootwards()
        return self.leafwards()

# For backward compatibility:
def parsimony_score(*args, **kwargs):
    return Parsimony(*args, **kwargs).leafwards()

# ~~> arrayal
def get_position_stats(align, nucl=False, allow_N=False):
    """Compute gap proportion and entropy value for each position of the alignment"""
    al_len = align.get_alignment_length()
    if al_len % 3 and not nucl:
        logger.error("Not a codon alignment!")

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

    logger.info(alint.shape)
    return gap_prop, al_entropy, alint


#from dendro.any import biophylo
def get_items_biophylo(tree, nodedist):
    return [(child, child.branch_length) for child in nodedist[0].clades]

def get_label_biophylo(tree, node):
    return node.name


# ~~> arrayal
def reorder_al(align, records=None, record_order=None):
    """
    param: `align`: Bio.Align object
    param: `record_order`: sequence of integers
    param: `records`: sequence of record names
    """
    # Mutually exclusive args
    if not (record_order is None) ^ (records is None):
        raise ValueError('`record_order` and `records` are mutually exclusive. Provide exactly one.')

    if records is not None:
        recnames = [rec.name for rec in align]
        record_order = [recnames.index(recordname) for recordname in records]

    return Align.MultipleSeqAlignment([align[r] for r in record_order])


def plot_al_stats(gap_prop, al_entropy, alint, dist_array=None, seqlabels=None,
                  outfile=None):
    """DEPRECATED."""
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
    logger.info(nvalues)
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
        axes[3].set_yticklabels(seqlabels, fontsize='x-small')
    axes[-1].set_xlabel("Residue position")

    if dist_array is not None:
        axes[4].step(x, dist_array.T, where='post', alpha=0.7)
        axes[4].legend(('pearson_corr', 'manhattan dist'), fontsize='x-small')
    
    if outfile is None:
        plt.show()
    else:
        fig.savefig(outfile)


def annotate_summary(ax, values):
    fmt = '%g' if len(np.asarray(values).shape) <= 1 else '%s'
    with np.printoptions(precision=4):
        summary = 'Mean: {0}\nMedian: {0}\nStd-dev: {0}'.format(fmt) % (
                    np.nanmean(values, axis=-1),  # axis=-1
                    np.nanmedian(values, axis=-1),
                    np.nanstd(values, axis=-1))
        #FIXME: no effect of precision when shape of arrays have dimension 1.
    ax.text(0.995, 0.99, summary, transform=ax.transAxes,
            verticalalignment='top', horizontalalignment='right')


PLOTS = ['gap_prop', 'entropy', 'gap_entropy', 'al', 'sp_score']
COMP_PLOTS = ['manhattan', 'pearson', 'split_pairs_score']

def al_colorbar(ax, cax=None, orientation='vertical', **kwargs):
    kwargs = {'shrink':0.05, **kwargs}
    im = ax.images[-1]
    cbar = plt.colorbar(im, ax=(ax if cax is None else None),
                        cax=cax,
                        orientation=orientation, extend='both', **kwargs)
    if cax is not None and orientation == 'horizontal':
        cax_pos = cax.get_position()
        cax.set_position(cax_pos.shrunk_to_aspect(1/20))
    cbar.set_label('residues')
    #ticks = cbar.get_ticks()
    #cbar.set_ticks()
    #nticks = len(ticks)
    logger.debug('cbar lim= %s; ylim= %s; xlim= %s', cbar.get_clim(),
                 cbar.ax.get_ylim(), cbar.ax.get_xlim())
    logger.debug('norm.vmax = %g; norm.N = %g; cmap.N = %g',
                 im.norm.vmax, im.norm.N, im.cmap.N)
    if im.norm.vmax > 5:
        sep = '\n' if orientation=='horizontal' else ' '
        cmax = im.cmap.N
        ticklabels = []
        prev_codon_col = im.cmap(im.norm(0))
        ticks = []
        logger.debug('codon_col 0: %s', prev_codon_col)
        for i, codon in enumerate(CODONS+STOPS, start=1):  #CODON2INT.items()
            codon_col = im.cmap(im.norm(i))
            #logger.debug('codon %d %r col %s; != %s', i, codon, codon_col,
            #             codon_col != prev_codon_col)
            if codon_col == prev_codon_col:
                ticklabels[-1].append(codon)
            else:
                ticks.append(i-1)
                try:
                    ticklabels[-1] = sep.join(ticklabels[-1])
                except IndexError:
                    pass
                ticklabels.append([])
                prev_codon_col = codon_col
        ticks.append(i-1)
        if not ticklabels[-1]:
            ticklabels[-1].append(codon)
        ticklabels[-1] = sep.join(ticklabels[-1])
        logger.debug('Setting %d ticks %s, %d labels %s', len(ticks), ticks,
                     len(ticklabels), ticklabels)
        cbar.set_ticks(ticks)
        cbar.set_ticklabels(ticklabels)  # [NACODON] + CODONS
        logger.debug('cbar lim= %s; ylim= %s; xlim= %s', cbar.get_clim(),
                     cbar.ax.get_ylim(), cbar.ax.get_xlim())
        #TODO
        cbar.ax.text(1, -0.05, gaps[0]*3, fontsize='x-small', clip_on=False) # under
        #cbar.ax.annotate(gaps[0]*3, (1, -0.05), xytext=(2, 0),
        #                 textcoords='offset points',
        #                 fontsize='x-small', annotation_clip=False) # under
    else:
        cbar.set_ticks(np.arange(len(nucleotides))+0.5)
        cbar.set_ticklabels(nucleotides)
    cbar.ax.tick_params(labelsize='x-small')
    return cbar


class AlignPlotter(object):

    #colorcycle = plt.rcParams['axes.prop_cycle']

    # FIXME: these should be instance attributes. Otherwise, when modified in 
    # an instance, they are changed for each newly created instance...
    @classmethod
    def set_plots(cls):
        plot_properties = {'al': {'title': 'Global scoring (all sequences)'},
                                  #'ylabel': 'Alignment'},
                           'gap': {'ylabel': 'Proportion of gaps\n(G)'},
                                   #'ylim': (-0.05,1.05)},
                           'pars': {'ylabel': 'Parsimony score\n(min number of changes per branch)'},
                           'entropy': {'ylabel': 'Entropy\n(H)'},
                           'gap_entropy': {'ylabel': 'Gap-entropy score:\n(1-H)*(1-G)'},
                           'sp': {'ylabel': 'SP score\n(Sum-of-pair differences)'},
                           'manh': {'title': 'Difference scoring between parts',
                                    'ylabel': 'manhattan distance'},
                           'pearson': {'ylabel': "Pearson's correlation coefficient"},
                           'part_sp': {'ylabel': 'sum of pair differences'},
                           'part_pars': {'ylabel': 'parsimony score within parts'},
                           'tree': {}}

        default_step = ('step', {'where': 'mid', 'alpha': 0.65})
        default_bar = ('bar', {'width': 1})
        default_stacked = (stackedbar, {'width': 1, 'edgecolor': 'none'})

        #from collections import namedtuple
        #funcArgs = namedtuple('funcArgs', 'func args')

        # { plot: [list of tuples(function, kwargs)] }
        plot_funcs = {'al':          [('imshow', {'aspect': 'auto'})], #(al_colorbar, {})],
                      #'al':          [('pcolormesh',   {'edgecolors': 'None'}),
                      #                ('invert_yaxis', {})],
                      #'gap':         default_step, #default_bar, "stackplot"
                      'gap':         [default_stacked],
                      'pars':        [default_bar, (annotate_summary,)],
                      'entropy':     [default_bar, (annotate_summary,)],
                      'gap_entropy': [default_bar, (annotate_summary,)],
                      'sp':          [default_bar],
                      'manh':        [default_step],
                      'pearson':     [default_step],
                      'part_sp':     [default_step],
                      'part_pars':   [default_stacked,
                                      (annotate_summary,),
                                      ('legend', {'bbox_to_anchor': (0,0,0.5,1)})],
                      'tree':        [(plottree, {'get_items': get_items_biophylo,
                                                  'get_label': get_label_biophylo,
                                                  'label_params': {'fontsize': 'x-small'}}),
                                      ('tick_params', {'labelright': False}),
                                      ('grid', {'b': False})]}
        return plot_properties, plot_funcs


    def __init__(self, alint, seqlabels=None, valid_range=(1,64), tree=None,
                 topology_only=False):
        """Initialize the Plotter instance from a matrix of integers.
        
        - alint: alignment provided as a numpy 2D array of integers.
        - tree: tree in Bio.Phylo format.
        """
        self.plot_properties, self.plot_funcs = self.set_plots()
        self.alint = alint
        self.x = np.arange(self.alint.shape[1]) # X values for plotting
        logger.debug('Alignment shape: %s', alint.shape)
        
        self.malint = np.ma.array(alint, mask=(alint==0))
        self.seqlabels = seqlabels
        self.plotlist = ['al']
        self.plotdata = {'al': (self.malint,)}
        #print(np.unique(self.malint))
        #print(np.unique(self.malint - self.malint.max()))
        self.valid_range = valid_range
        self.tree = tree
        if tree is not None:
            self.plotlist.insert(0, 'tree')
            self.plotdata['tree'] = (tree,)
            self.plot_funcs['tree'][0][1].update(topology_only=topology_only)

        #cmap_size = valid_range[1] - valid_range[0] # If len(valid_range) > 1
        # Dark2, tab20b, tab20, Set1
        #FIXME: valid_range for codons should exclude the stops: (1,61)
        cmap = plt.get_cmap('tab20b', valid_range[1] - valid_range[0] + 1) #, cmap_size)
        cmap.set_over((0.2, 0.2, 0.2))
        cmap.set_under('w')
        # a norm is needed to set bounds!!!
        
        # Directly use integer values to get colors in cmap: Same speed but less clear.
        #self.plotdata = {'al': (self.malint - valid_range[0],)} # With NoNorm
        #norm = mpl.colors.NoNorm(0, valid_range[1] - valid_range[0])

        bounds = np.arange(valid_range[0]-0.5, valid_range[1] + 0.6)
        norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
        # Maybe a ListedColormap is better.
        self.plot_funcs['al'][0][1].update(cmap=cmap, norm=norm)
        self.plot_properties['al']['yticks'] = np.arange(alint.shape[0])
        #TODO: Add the *last* xtick showing alignment length.
        #self.plot_properties['al']['xlim'] = (0, alint.shape[1]) #pcolormesh
        self.plot_properties['al']['xlim'] = (-0.5, alint.shape[1]-0.5) #imshow
        #TODO: Add plt.colorbar for residue colors.

        if seqlabels is not None:
            self.plot_funcs['al'].append(('set_yticklabels',
                                          {'labels': seqlabels,
                                           'fontsize': 'x-small'}))

    @classmethod
    def fromfile(cls, infile, format=None, nucl=False, allow_N=False,
                 ungap=True, slice=None, records=None, recordsfile=None,
                 treefile=None, topology_only=False):
        align = AlignIO.read(infile,
                             format=(format or
                                     filename2format(getattr(infile, 'name',
                                                             infile)))
                            )

        al_len = align.get_alignment_length()
        if al_len % 3 and not nucl:
            logger.error("Not a codon alignment!")

        tree = None

        if records or recordsfile or treefile:
            record_order = None
            if records:
                record_order = [int(r) for r in records.split(',')]
                records = None
            elif recordsfile:
                with open(recordsfile) as recf:
                    records = [line.rstrip() for line in recf]
            elif treefile:
                #try:
                #    import ete3
                #    tree = ete3.Tree(treefile, format=1)
                #    get_leaves = tree.get_leaves
                #except ImportError:
                from Bio.Phylo import read as phyloread
                tree = phyloread(treefile, 'newick')
                get_leaves = tree.get_terminals
                records = [leaf.name for leaf in get_leaves()]

            align = reorder_al(align, records, record_order)

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

        # ungap: remove columns being 100% gaps
        if ungap:
            alint = alint[:, alint.sum(axis=0) > 0]

        # Init and setup the class instance
        valid_range = (1, 4) if nucl else (1, 64)
        instance = cls(alint, seqlabels, valid_range, tree, topology_only)
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
        self.plot_funcs['entropy'][1] = self.plot_funcs['entropy'][1] \
                                     + ({'values': self.entropy},) 
        self.plot_funcs['gap_entropy'][1] = self.plot_funcs['gap_entropy'][1] \
                                     + ({'values': self.gap_entropy},) 

        if self.tree is not None:
            if hasattr(self.tree, 'children'):
                get_children = lambda tree, node: node.children
            else:
                get_children = lambda tree, node: node.clades
            self.parsimony_score = parsimony_score(self.alint, self.tree,
                                                   seqlabels=self.seqlabels,
                                                   minlength=(max_valid+2),
                                                   get_children=get_children)
            self.plotlist.insert(self.plotlist.index('al')+1, 'pars')
            self.plotdata['pars'] = (self.x, self.parsimony_score)
            # annotate_summary:
            self.plot_funcs['pars'][1] += ({'values': self.parsimony_score},)


    def annotate_parts(self, ax):
        """Draw arrows with annotation to show 'compare_parts' on the
        alignment plot."""
        #x0, x1 = ax.get_xlim()
        #xtext = x1 + (x1-x0) * 1.05
        x1, xtext = 1, 1.02
        hlines = []
        for part_number, part in enumerate(self.parts, start=1):
            ytext = np.mean(part) / self.alint.shape[0]
            ax.annotate(str(part_number), xy=(x1, ytext), xytext=(xtext, ytext),
                        xycoords='axes fraction', alpha=0.8)
            for row in part:
                ax.annotate("", xy=(x1, (row+0.5)/self.alint.shape[0]),
                            xytext=(xtext, ytext),
                            xycoords='axes fraction',
                            arrowprops={'arrowstyle': '->', 'alpha':0.8},
                            alpha=0.8)
            # if part is contiguous:
            row0, row_end = min(part), max(part)+1
            if set(part) == set(range(row0, row_end)):
                if row0 not in hlines:
                    ax.axhline(row0-0.5, linestyle='--', color='#191919', alpha=0.5)
                    hlines.append(row0)
                if row_end not in hlines:
                    ax.axhline(row_end-0.5, linestyle='--', color='#191919', alpha=0.5)
                    hlines.append(row_end)


    def comp_parts(self, compare_parts):
        """Compute column-wise metrics of difference between two parts (groups
        of records) of the dataset:
            - manhattan distance between frequencies
            - pearson coefficient between counts
            - sum of pair differences between any possible pair of distinct parts
        """
        ### TODO: add a parsimony comparison: 1 if null intersection.
        if isinstance(compare_parts, str):
            compare_parts = [part.rstrip(')').lstrip('(').split(',')
                             for part in compare_parts.split(';')]
        # Convert to the proper type
        self.parts = []
        for part in compare_parts:
            parti = []
            for subpart in part:
                try:
                    parti.append(int(subpart))
                except ValueError:
                    try:
                        parti.extend(range(*(int(i) for i in subpart.split(':'))))
                    except ValueError:
                        label_reg = re.compile(subpart)
                        parti.extend((i for i,lab in enumerate(self.seqlabels)
                                      if label_reg.match(lab)))
            self.parts.append(parti)

        # The complementary sequences.
        out_part = list(set(range(self.alint.shape[0]))
                        .difference(set().union(*(self.parts))))
        all_parts = self.parts + [out_part]

        #assert len(self.parts) == 2

        self.part_al = [self.alint[part,:] for part in all_parts[:2]]

        self.part_freqs = fmat1, fmat2  = [freq_matrix(alpart) for alpart in self.part_al]
        self.part_counts = cmat1, cmat2 = [count_matrix(alpart) for alpart in self.part_al]

        self.part_manh_dist = np.abs(fmat1 - fmat2).sum(axis=0)
        self.part_pearson_c = pearson_coeff(cmat1, cmat2, axis=0)
        self.part_sp_score = part_sp_score(self.alint, all_parts[:2])

        self.plotlist.extend(('manh', 'pearson', 'part_sp'))
        self.plotdata.update(manh=(self.x, self.part_manh_dist,),
                             pearson=(self.x, self.part_pearson_c,),
                             part_sp=(self.x, self.part_sp_score,))

        self.plot_funcs['al'].append((self.annotate_parts, {}))

        if self.tree is not None: # and 'part_pars' in self.plotlist:
            if hasattr(self.tree, 'children'):
                get_children = lambda tree, node: node.children
            else:
                get_children = lambda tree, node: node.clades
            self.part_pars_scores, part_branch_nbs = parsimony_score(self.alint, self.tree,
                                                seqlabels=self.seqlabels,
                                                minlength=(self.valid_range[1]+2),
                                                get_children=get_children,
                                                parts=self.parts)
            logger.debug('Finished partial parsimony scores')
            self.plotlist.append('part_pars')
            #self.plotdata['part_pars'] = tuple((arr for p_sc in self.part_pars_scores
            #                                    for arr in (self.x, p_sc)))
            self.plotdata['part_pars'] = (self.x, self.part_pars_scores)

            #self.plot_funcs['part_pars'][0][1].update(alpha=0.6)
            self.plot_funcs['part_pars'][2][1].update(labels=['outgroup']+
                                                      list(range(1, 1+len(self.parts))))
            self.plot_funcs['part_pars'][1] += ({'values': self.part_pars_scores},)


    def makefig(self, plotlist=None, figwidth=16, plotheight=4):

        plotlist = plotlist if plotlist is not None else self.plotlist
        nplots = len(plotlist)
        fig = plt.figure(figsize=(figwidth, plotheight*nplots))

        logger.debug('Nplots: %s, %s', nplots, plotlist)
        #fig, axes = plt.subplots(nplots, sharex=True,
        #                         figsize=(figwidth, 4*nplots), squeeze=False)
        
        if 'tree' in plotlist:
            rows = nplots - 1
            cols = 5  # tree plot spans 1 col, alignment spans 4.
            colspan = cols - 1
            gridspec_kw = {'width_ratios': [1, colspan]}
            # And hide alignment labels
            #self.plot_funcs['al'].append(('tick_params', {'label_left'}))
        else:
            rows = nplots
            cols = 1
            colspan = 1
            gridspec_kw = {}
        logger.debug("subplots: %d, %d" % (rows, cols))
        pos = (0,0)
        #ax = fig.add_subplot(rows, cols, pos, gridspec_kw=gridspec_kw)
        ax = plt.subplot2grid((rows,cols), pos, colspan=1, fig=fig)
        axes = [ax]

        for plot in plotlist:
            logger.debug('Plot=%r; pos=%s', plot, pos)
            logger.debug('len(axes): %d', len(axes))

            plotfuncname, plot_kwargs = self.plot_funcs[plot][0]
            try:
                plotfunc = getattr(ax, plotfuncname)
            except TypeError:
                plotfunc = plotfuncname
                # ax argument required:
                plot_kwargs.update(ax=ax)

            plotdata = self.plotdata[plot]

            plotfunc(*plotdata, **plot_kwargs)

            for i, (extrafunc, extra_kwargs) in enumerate(self.plot_funcs[plot][1:],
                                                          start=1):
                try:
                    if isinstance(extrafunc, str):
                        extrafunc = getattr(ax, extrafunc)
                    else:
                        extra_kwargs.update(ax=ax)
                        #extrafunc = lambda **kw: extrafuncname(ax, **kw)
                    extrafunc(**extra_kwargs)
                except BaseException as err:
                    err.args += ('extra func %d %r' % (i, extrafunc),)
                    raise

            plot_prop = self.plot_properties[plot]
            ax.set(**plot_prop)

            # Next ax:
            if len(axes) < nplots:
                if plot == 'tree':
                    #pos += 1
                    #ax = fig.add_subplot(rows, cols, pos, sharey=ax)
                    pos = (pos[0], pos[1]+1)
                    ax = plt.subplot2grid((rows,cols), pos, colspan=colspan, fig=fig, sharey=ax, autoscalex_on=False)
                else:
                    #pos += 2 if 'tree' in plotlist else 1
                    #ax = fig.add_subplot(rows, cols, pos, sharex=ax)
                    if plot == 'al' and rows>1 and cols>1:
                        # Plot the colorbar below the tree
                        cax = plt.subplot2grid((rows,cols),(pos[0]+1, 0), fig=fig,
                                               aspect=20)
                        al_colorbar(ax, cax)#, orientation='horizontal')
                    pos = (pos[0]+1, pos[1])
                    logger.debug('rows,cols: (%d,%d); pos: %s; colspan: %d',
                                 rows, cols, pos, colspan)
                    ax = plt.subplot2grid((rows,cols), pos, colspan=colspan, fig=fig, sharex=ax, autoscalex_on=False)
                axes.append(ax)

        # On the last plot:
        if self.slstart:
            # Update the displayed xticklabels if start position > 0
            slstart = self.slstart // 3 if not self.nucl else self.slstart
                
            ax.set_xticklabels([(slstart + int(xt)) for xt in \
                                ax.get_xticks()])

        fig.tight_layout()
        self.fig, self.axes = fig, axes


    def display(self, outfile=None, **kwargs):
        logger.debug('Displaying')
        if outfile is None:
            plt.show(block=True)
            #print(clock())
        else:
            self.fig.savefig(outfile, bbox_inches='tight', **kwargs)



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
    gap_prop, al_entropy, alint = get_position_stats(align, nucl, allow_N)
    dist_array = comp_parts(alint, compare_parts)
    plot_al_stats(gap_prop, al_entropy, alint, dist_array, seqlabels, outfile)


def plot_al_conservation(infile, format=None, nucl=False, allow_N=False,
         ungap=True, records=None, recordsfile=None, treefile=None, slice=None,
         topology_only=False, compare_parts=None,
         compare_only=False, plotlist=None, figwidth=16, plotheight=4):

    align_plot = AlignPlotter.fromfile(infile, format, nucl, allow_N, ungap,
                                       slice, records, recordsfile, treefile,
                                       topology_only)
    if logger.getEffectiveLevel() <= logging.DEBUG:
        del align_plot.plot_funcs['tree'][1]  # tickparams(labelright=False)
    if not compare_only:
        align_plot.measure()
    if compare_parts:
        align_plot.comp_parts(compare_parts)
    if plotlist is not None:
        plotlist = plotlist.split(',')
        if set(plotlist).intersection(('pars', 'tree', 'part_pars')):
            assert treefile is not None, \
                    "A treefile must be given to compute parsimony score or plot a tree."
            # plotlist.remobe('pars') ## 'tree'
            if 'tree' not in plotlist and 'al' in plotlist:
                plotlist.insert(plotlist.index('al'), 'tree')
    align_plot.makefig(plotlist, figwidth, plotheight)
    return align_plot


def main():
    logging.basicConfig(format="%(levelname)s:%(funcName)s:%(message)s")
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
    parser.add_argument('-g', '--no-ungap', dest='ungap', action='store_false',
                        help='Do not remove columns of gaps.')
    rec_g = parser.add_mutually_exclusive_group()
    rec_g.add_argument('-r', '--records',
                        help='select records (coma-sep list of 0-based integers)')
    rec_g.add_argument('-R', '--recordsfile',
                        help='select records from a file (one record name per line)')
    rec_g.add_argument('-t', '--treefile',
                        help='select records from a newick tree (reorder the sequences)')
    parser.add_argument('-T', '--topology-only', action='store_true',
                        help='Plot unit branch lengths instead of the real ones.')
    parser.add_argument('-s', '--slice',
                        help='select positions (start:end). '\
                             '0-based, end excluded, in number of '\
                             'codons/nucleotides.')
    parser.add_argument('-c', '--compare-parts',
                        help='Plot the per-column correlation between two ' \
                             'groups of data, e.g "(0,1);(2,3,4)". The last '\
                             'part is calculated automatically. slices and '\
                             'regex (matching sequence names) allowed')
    parser.add_argument('-C', '--compare-only', action='store_true',
                        help='Do not display global scores')
    parser.add_argument('-W', '--figwidth', default=16, type=float,
                        help='Figure width (inches) [%(default)s]')
    parser.add_argument('-H', '--plotheight', default=1.8, type=float,
                        help='Individual plot height (inches) [%(default)s]')
    parser.add_argument('-p', '--plotlist', default='al,gap,pars',
                        help='comma-sep list of plots. Valid values are: '\
                             'al,gap,pars,entropy,gap_entropy,sp,manh,pearson'
                             'part_sp,part_pars. [%(default)s]')
    parser.add_argument('-v', '--verbose', action='count', default=0)
    
    args = parser.parse_args()

    if args.verbose > 1:
        logger.setLevel(logging.DEBUG)
    elif args.verbose > 0:
        logger.setLevel(logging.INFO)
    delattr(args, 'verbose')

    outfile = args.outfile
    delattr(args, 'outfile')
    if not outfile:
        plt.switch_backend('TkAgg')

    plot_al_conservation(**vars(args)).display(outfile)


if __name__ == '__main__':
    main()
