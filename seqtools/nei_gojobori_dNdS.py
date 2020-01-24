#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Compute the dN/dS ratio between 2 sequences using the Nei & Gojobori 1986 method:
unweighted pathway analysis.
"""


from __future__ import print_function, division


from sys import stderr
import argparse as ap
from itertools import permutations
from math import log, sqrt, exp, nan
import random

# Needed for calling YN00 (PAML)
import os.path as op
import tempfile
from shutil import rmtree
import subprocess as sp
import re

import logging
logger = logging.getLogger(__name__)


stop = '*'

# Might have a difference in my genetic code with PAML's genetic codes.
genetic_code = {
    'TTT': 'Phe', 'TTC': 'Phe', 'TTA': 'Leu', 'TTG': 'Leu',
    'TCT': 'Ser', 'TCC': 'Ser', 'TCA': 'Ser', 'TCG': 'Ser',
    'TAT': 'Tyr', 'TAC': 'Tyr', 'TAA': stop,  'TAG': stop,
    'TGT': 'Cys', 'TGC': 'Cys', 'TGA': stop,  'TGG': 'Trp',
    'CTT': 'Leu', 'CTC': 'Leu', 'CTA': 'Leu', 'CTG': 'Leu',
    'CCT': 'Pro', 'CCC': 'Pro', 'CCA': 'Pro', 'CCG': 'Pro',
    'CAT': 'His', 'CAC': 'His', 'CAA': 'Gln', 'CAG': 'Gln',
    'CGT': 'Arg', 'CGC': 'Arg', 'CGA': 'Arg', 'CGG': 'Arg',
    'ATT': 'Ile', 'ATC': 'Ile', 'ATA': 'Ile', 'ATG': 'Met',
    'ACT': 'Thr', 'ACC': 'Thr', 'ACA': 'Thr', 'ACG': 'Thr',
    'AAT': 'Asn', 'AAC': 'Asn', 'AAA': 'Lys', 'AAG': 'Lys',
    'AGT': 'Ser', 'AGC': 'Ser', 'AGA': 'Arg', 'AGG': 'Arg',
    'GTT': 'Val', 'GTC': 'Val', 'GTA': 'Val', 'GTG': 'Val',
    'GCT': 'Ala', 'GCC': 'Ala', 'GCA': 'Ala', 'GCG': 'Ala',
    'GAT': 'Asp', 'GAC': 'Asp', 'GAA': 'Glu', 'GAG': 'Glu',
    'GGT': 'Gly', 'GGC': 'Gly', 'GGA': 'Gly', 'GGG': 'Gly'}



def read_seq(alfile):
    """Return a list of (label, sequence)."""
    sequences = []
    current_seq = None
    with open(alfile) as f:
        line0 = next(f)
        if line0.startswith('>'):
            seq_label = line[1:].strip()
            current_seq = ''
            for line in f:
                if line.startswith('>'):
                    if current_seq is not None:
                        sequences.append((seq_label, current_seq))
                    seq_label = line[1:].strip()
                    current_seq = ''
                else:
                    try:
                        current_seq += line.rstrip().replace(' ', '')
                    except TypeError as err:
                        err.args = ((err.args[0]
                                     + '. Received sequence before label: wrong file format? (need Phylip/Fasta)',)
                                    + err.args[1:])
                        raise
            sequences.append((seq_label, current_seq))
        else:
            try:
                ns, ls = [int(x) for x in line0.split()]  # try Phylip format (sequential)
            except ValueError:
                err.args = ((err.args[0] + '. Unrecognized format (need Phylip/Fasta)',)
                            + err.args[1:])
                raise
            for line in f:
                sequences.append(line.strip().split(maxsplit=1))
                assert len(sequences[-1][1]) == ls, "Wrong sequence length"
            assert len(sequences) == ns, "Wrong number of sequences"
    return sequences


def frac_synonymous(codon, i):
    aa = genetic_code[codon]
    tot_subst = 0
    tot_syn = 0
    for other_nucleotide in set('ATCG').difference((codon[i],)):
        other_codon = codon[:i] + other_nucleotide + codon[(i+1):]
        if genetic_code[other_codon] is not stop:
            tot_subst += 1
            if genetic_code[other_codon] == aa:
                tot_syn += 1
    return tot_syn/tot_subst


def nei_gojobori(seq1, seq2):
    """Return sequence-wise raw values for computing the dN/dS: L, S1, S2, Nobs, Sobs
    
    L: length of both sequences (nucleotides)
    S1: number of synonymous sites from 1st sequence to 2d sequence (computed by NG86)
    S2: number of synonymous sites from 2d sequence to 1st (computed by NG86)
    Nobs: number of observed non-synonymous differences
    Sobs: number of observed synonymous differences.

    Unweighted pathway method (Nei & Gojobori, 1986).
    """
    assert len(seq1) == len(seq2)
    L = len(seq1)
    assert L % 3 == 0, "Number of nucleotides is not a multiple of 3."

    S1, S2 = 0, 0  # Number of sequence sites potentially synonymous.
    Sobs = 0  # Number of observed synonymous differences.
    Nobs = 0
    for k in range(0, L, 3):
        #print('%d: ' % k, file=stderr)
        codon1 = seq1[k:(k+3)]
        codon2 = seq2[k:(k+3)]
        if codon1 == '---' or codon2 == '---':
            L -= 1
            continue

        # n: number of potential nonsynonymous sites per codon (fraction of
        #    possible 1 nucleotide substitution being non synonymous)
        # s: number of potential synonymous sites per codon = 3 - n
        S1 += sum(frac_synonymous(codon1, i) for i in range(3))
        S2 += sum(frac_synonymous(codon2, i) for i in range(3))
        #print('Frac synonymous 1:', *(frac_synonymous(codon1, i) for i in range(3)), file=stderr)
        #print('Frac synonymous 2:', *(frac_synonymous(codon2, i) for i in range(3)), file=stderr)
        
        # Observed differences
        changed = [i for i, (nucl1, nucl2) in enumerate(zip(codon1, codon2)) if nucl1 != nucl2]
        ndiff = len(changed)
        if ndiff == 1:
            sobs = (genetic_code[codon1] == genetic_code[codon2])
        elif ndiff > 1:
            # Enumerate possible pathways of single-nucleotide changes:
            # For 2 differences, there are 2 pathways, depending on which
            # nucleotide changes first;
            # for 3 differences, there are 6 pathways.
            paths_sobs = []
            paths_nobs = []
            for order in permutations(changed):
                path_sobs = 0
                prev_step_codon = list(codon1)
                for i in order:
                    step_codon = list(prev_step_codon)  # copy the list.
                    step_codon[i] = codon2[i]
                    if genetic_code[''.join(step_codon)] is stop:
                        break
                    path_sobs += (genetic_code[''.join(prev_step_codon)]
                                  == genetic_code[''.join(step_codon)])
                    prev_step_codon = step_codon
                else:
                    paths_sobs.append(path_sobs)
            sobs = sum(paths_sobs) / len(paths_sobs)
        if ndiff:
            Sobs += sobs
            Nobs += ndiff - sobs

    return L, S1, S2, Nobs, Sobs
    #mean_S = (S1+S2)/2
    #return Nobs / (L - mean_S), Sobs / mean_S


def pNpS(L, S1, S2, Nobs, Sobs):
    mean_S = (S1+S2)/2
    return Nobs / (L - mean_S), Sobs / mean_S

def dNdS(L, S1, S2, Nobs, Sobs):
    mean_S = (S1+S2)/2
    pN, pS = Nobs / (L - mean_S), Sobs / mean_S
    return jukes_cantor(pN), jukes_cantor(pS)


def jukes_cantor(p):
    """Given p observed differences, return the evolutionary distance."""
    try:
        return -3/4 * log(1 - 4/3 * p)
    except ValueError as err:
        if "math domain error" in err.args[0]:
            err.args += ("p=%g" % p,)
            print(('WARNING: p=%g >= 3/4 in Jukes-Cantor correction '
                   '(infinite distance)') % p,
                    file=stderr)
            #return nan
        else:
            raise


def inv_jukes_cantor(d):
    return (1 - exp(- d*4/3)) * 3/4


def variance_of_divergence():
    pass


def resample(seq):
    """sample with replacement"""
    return [random.choice(seq) for position in seq]


# No, I won't be using Numpy
def variance(x):
    n = len(x)
    m = sum(x)/n
    return sum((xi - m)**2 for xi in x)/(n-1)*n


def standard_deviation(x):
    return sqrt(variance(v))


def bootstrap_dNdS(seq1, seq2):
    L1, L2 = len(seq1), len(seq2)
    seq1_codons = [seq1[i:(i+3)] for i in range(0, L1, 3)]
    seq2_codons = [seq2[i:(i+3)] for i in range(0, L2, 3)]
    
    all_dS = []
    all_dN = []
    for k in range(1000):
        newseq1 = ''.join(resample(seq1_codons))
        newseq2 = ''.join(resample(seq2_codons))
        dN, dS = dNdS(*nei_gojobori(newseq1, newseq2))
        if dN is not None and dS is not None:
            all_dN.append(dN)
            all_dS.append(dS)
    
    return variance(all_dN), variance(all_dS)


def detail_nei_gojobori(seq1, seq2):
    """Print dNdS computations with intermediate values"""
    L, S1, S2, Nobs, Sobs = nei_gojobori(seq1, seq2)
    print('L=%d S1=%g S2=%g Nobs=%d Sobs=%d' % (L,S1,S2,Nobs,Sobs))

    mean_S = (S1+S2)/2
    pN, pS = Nobs / (L - mean_S), Sobs / mean_S
    print('pN=%g, pS=%g' % (pN, pS))

    dN, dS = jukes_cantor(pN), jukes_cantor(pS)
    try:
        ratio = '%g' % (dN/dS)
    except ZeroDivisionError:
        ratio = nan
    print('dN=%g, dS=%g;  dN/dS = %s' %(dN, dS, ratio))

    #var_dN, var_dS = bootstrap_dNdS(seq1, seq2)

    #Z = (dN-dS) / sqrt(var_dN + var_dS)
    #print('Standard score Z(dN - dS) = %g' % Z)


def write_phy(filehandle, *seqs):
    filehandle.write('\t%d %d\n' % (len(seqs), len(seqs[0].replace('\n', ''))))
    for i, seq in enumerate(seqs, start=1):
        filehandle.write('seq%d  %s\n' % (i, seq.replace('\n', '')))


#@contextlib.contextmanager
def set_tmp_files_yn00(seq1, seq2):
    # Create the seq files
    #with tempfile.TemporaryDirectory(prefix='yn00_') as tmploc:
    tmploc = tempfile.mkdtemp(prefix='yn00_')
    with tempfile.NamedTemporaryFile('w', suffix='.phy', dir=tmploc,
                                     delete=False) as tmpalfile:
        tmpbasename,_ = op.splitext(op.basename(tmpalfile.name))
        write_phy(tmpalfile, seq1, seq2)

    tmproot = op.join(tmploc, tmpbasename)
    with open(tmproot + '.ctl', 'w') as tmpctlfile:
        tmpctlfile.write("""
seqfile = {tmpal}  * sequence data file name
outfile = {tmproot}.out  * main result file
verbose = 0     * 1: detailed output (list sequences), 0: concise output
icode = 0       * 0:universal code; 1:mammalian mt; 2-10:see below
weighting = 0   * weighting pathways between codons (0/1)?
commonf3x4 = 0  * use one set of codon freqs for all pairs (0/1)?
""".format(tmpal=tmpalfile.name, tmproot=tmproot))
    return tmploc, tmproot+'.ctl', tmproot+'.out'


def nei_gojobori_yn00(seq1=None, seq2=None, ctl=None, out=None, cleanup=True):
    noseq = seq1 is None or seq2 is None
    nofile = ctl is None or out is None
    logger.debug('noseq=%s nofile=%s cleanup=%s', noseq, nofile, cleanup)
    assert (not noseq and nofile) or \
            (noseq is None and not nofile), \
            "(seq1, seq2) and (ctl, out) are mutually exclusive"
    if nofile:
        tmploc, ctl, out = set_tmp_files_yn00(seq1, seq2)

    dNdS, dN, dS = call_nei_gojobori_yn00(ctl, out)

    if nofile and cleanup:
        logger.debug('Removing tmp dir %r' % tmploc)
        rmtree(tmploc)
    return dNdS, dN, dS


def call_nei_gojobori_yn00(ctl, out):
    """Use the yn00 program from PAML (Z. Yang)"""
    sp.check_call(['yn00', ctl])

    reg_nsls = re.compile(r'ns\s*=\s*(\d+)\tls\s*=(\d+)')
    reg_ynmethod = re.compile(re.escape("(B) Yang & Nielsen (2000) method"))
    reg_ngmethod = re.compile(re.escape("Nei & Gojobori 1986. dN/dS (dN, dS)"))
    reg_result = re.compile(r'\w+\s+(-?\d+\.\d+)\s*\((-?\d+\.\d+)\s+(-?\d+\.\d+)\)')
    with open(out) as outfile:
        for line in outfile:
            if reg_ynmethod.match(line):
                break
            if reg_nsls.match(line):
                ns, ls = reg_nsls.match(line).groups()
            if reg_ngmethod.match(line):
                for i in range(5):
                    nextline = next(outfile)
                    if nextline.startswith('seq1'):
                        seq1_line = nextline.rstrip()
                    elif nextline.startswith('seq2'):
                        seq2_line = nextline.rstrip()
                        break
    dNdS, dN, dS = [float(x) for x in reg_result.match(seq2_line).groups()]
    return dNdS, dN, dS


def main():
    parser = ap.ArgumentParser(description=__doc__)
    parser.add_argument('alfile',
                        help='Alignment file containing 2 coding sequences of nucleotides (fasta/phylip).')
    args = parser.parse_args()
    (_, seq1), (_, seq2) = read_seq(args.alfile)
    detail_nei_gojobori(seq1, seq2)


if __name__ == '__main__':
    logging.basicConfig()
    main()
