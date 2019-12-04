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
from math import log, sqrt
import random


stop = '*'

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



def read_fasta(alfile):
    """Return a list of (label, sequence)."""
    sequences = []
    current_seq = None
    with open(alfile) as f:
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
                                 + '. Received sequence before label: wrong file format? (need Fasta)',)
                                + err.args[1:])
                    raise
    sequences.append((seq_label, current_seq))
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
    """Return dN, dS by the unweighted pathway method (Nei & Gojobori, 1986)."""
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


def dNdS(L, S1, S2, Nobs, Sobs):
    mean_S = (S1+S2)/2
    pN, pS = Nobs / (L - mean_S), Sobs / mean_S
    return jukes_cantor(pN), jukes_cantor(pS)


def jukes_cantor(p):
    """Given p observed differences, return the evolutionary distance."""
    return -3/4 * log(1 - 4/3 * p)


def variance_of_divergence():
    pass


def resample(seq):
    """sample with replacement"""
    return [random.choice(seq) for _ in seq]


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
        all_dN.append(dN)
        all_dS.append(dS)
    
    return variance(all_dN), variance(all_dS)


def main():
    parser = ap.ArgumentParser(description=__doc__)
    parser.add_argument('alfile',
                        help='Alignment file containing 2 coding sequences of nucleotides (fasta).')
    args = parser.parse_args()
    (_, seq1), (_, seq2) = read_fasta(args.alfile)
    L, S1, S2, Nobs, Sobs = nei_gojobori(seq1, seq2)

    mean_S = (S1+S2)/2
    pN, pS = Nobs / (L - mean_S), Sobs / mean_S
    print('pN=%g, pS=%g' % (pN, pS))

    dN, dS = jukes_cantor(pN), jukes_cantor(pS)
    try:
        ratio = '%g' % (dN/dS)
    except ZeroDivisionError:
        ratio = None
    print('dN/dS = %s  (dN=%g, dS=%g)' %(ratio, dN, dS))

    var_dN, var_dS = bootstrap_dNdS(seq1, seq2)

    Z = (dN-dS) / sqrt(var_dN + var_dS)
    print('Standard score Z(dN - dS) = %g' % Z)


def test():
    assert round(frac_synonymous('TTT', 2), 2) == 0.33
    assert round(frac_synonymous('TTT', 1), 2) == 0
    assert round(frac_synonymous('TTT', 0), 2) == 0
    assert round(frac_synonymous('TGT', 2), 2) == 0.5
    assert round(frac_synonymous('TGT', 1), 2) == 0
    assert round(frac_synonymous('TGT', 0), 2) == 0
    assert nei_gojobori('TTT', 'GTA') == (3, 1/3, 1, 1.5, 0.5), nei_gojobori('TTT', 'GTA')
    assert nei_gojobori('TTG', 'AGA') == (3, 2/3, 1/2+1/3, 9/4, 3/4), nei_gojobori('TTG', 'AGA')

    # Ensembl 93 ENSG00000186260 ENSPANG00000020098
    seq1 = ("ATGGATCACACAGGGGCGATAGACACCGAGGATGAAGTGGGACCTTTAGCCCATCTTGCTCCAAGTCCTCAGAGTGAAGCTGTGGCTCATGAATTCCAGGAACTCTCCTTGCAGTCCAGTCAAAACTTACCCCCTCTGAACGAAAGGAAAAATGTGCTCCAGCTGAGGCTGCAACAAAGGAGGACGAGAGAACAACTAGTGGACCAGGGCATCATGCCACCTTTGAAGAGCCCAGCGGCATTCCATGAACAGATAAAAAGCTTGGAACGAGCCAGAACTGAAAACTTTTTGAAACACAAGATTCGGAGTCGACCAGATCGTTCTGAACTTGTCAGGATGCACATTTTAGAAGAAACATTTGCAGAGCCATCCCTGCAGGCTACTCAGATGAAGTTGAAAAGAGCTCGACTAGCAGATGATCTGAATGAAAAGATTGCTCAAAGACCTGGTCCTATGGAGCTGGTAGAGAAAAACATCCTTCCTGTGGACTCCAGTGTTAAAGAAGCAATTATAGGCGTTGGGAAGGAGGACTATCCCCACACTCAGGGCGATTTCTCATTTGATGAAGACAGCAGTGACGCTTTGTCTCCGGACCAGCCTGCGAGTCAGGAGTCACAGGGGTCAGCCGCGTCCCCAAGTGAGCCAAAAGTTAGTGAATCGCCATCTCCTGTGACTACAAACACTCCAGCGCAGTTTGCTTCAGTGTCCCCAACAGTTCCTGAATTCTTGAAAACTCCTCCAACTGCAGATCAGCCTCCCCCACGGCCTGCAGCTCCTGTCCTCCCCACAAACACTGTGTCCTCAGCAAAGCCTGGCCCAGCACTGGTGAAGCAAAGCCATCCCAAGAATCCAAATGACAAACACCGTAGCAAAAAGTGCAAAGATCCCAAACCACGGGTAAAGAAGTTAAAGTACCACCAATACATTCCACCAGATCAGAAGGGTGAGAAGAATGAGCCGCAGATGGACTCTAACTACGCCCGCCTGCTCCAGCAGCAGCA"
            "GCTGTTCCTGCAACTGCAGATCCTGAGTCAGCAGAAGCAGCACTACAACTACCAGACCATCCTGCCTGCACCATTCAAGCCACTCAATGACAAAAATAGTAACAGTGGGAATTCAGCTTTGAACAATGCCACACCTAACACACCAAGACAGAATACATCTACTCCTGTGAGAAAGCCAGGACCTCTGCCTTCTAGCCTGGATGACTTAAAGGTATCAGAACTGAAGACAGAACTGAAGTTAAGGGGTCTGCCAGTGTCAGGCACCAAACCGGACCTCATTGAGCGCCTAAAACCCTACCAGGAAGTGAACAGCAGCGGCCTTGCTGCTGGGGGCATCGTGGCAGTGTCATCATCAGCCATTGTCACCAGTAACCCAGAAGTCACTGTGGCCTTGCCGGTTACAACACTACACAACACTGTGACTAGCTCAGTCTCTACTCTCAAGGCAGAATTGCCACCTACAGGAACCAGCAACGCAACCCGTGTGGAAAATGTTCATTCCCCTCTGCCCATTTCACCATCTCCCTCCGAACAGTCCAGTCTCAGTACTGATGACACAAACATGGCAGACACTTTCACCGAGATTATGACCATGATGTCGCCTTCACAGTTCTTGAGTTCATCTCCTTTGAGAATGACAAATAATGAAGACAGTCTGAGTCCCACCAGCAGCACTCTGTCAAACCTGGAACTGGATGCAGCCGAAAAGGATCGCAAGCTTCAGGAGAAAGAGAAGCAAATCGAAGAGCTGAAGAGGAAACTGGAACAAGAGCAGAAGCTCGTGGAAGTGCTGAAAATGCAACTTGAGGTTGAAAAACGAGGGCAGCAGCAGCGGCCCCTGGAAGCCCAGCCCAGTGCCCCAGGTCATTCTGTC---AAGTCAGATCAGAAGCACGGCAGCCTTGGCTCCTCCATCAAAGATGAGGCCTCACTCCCTGACTGCTCCAGCTCCAGGCAGCCCATCCCAGTAGCCAGCCACGCTGTAGGCCAGCCCGTCTCTA"
            "CAGGTGGCCAGACCCTTGTTGCCAAAAAGGCTGTAGTTATCAAGCAAGAGGTCCCTGTGGGCCAGGCAGAGCAGCAGAGTGTCGTCTCGCAGTTTTATGTGAGTTCCCAGGGACAGCCACCGCCTGCTGTTGTTGCTCAGCCCCAGGCTTTACTGACCACGCAGACTGCTCAGCTGCTGCTCCCAGTGTCCATCCAGGGCTCGAGTGTCACCTCAGTGCAACTCCCTGTAGGCAGCCTCAAACTCCAGACTTCACCACAAGCAGGAATGCAGACTCAGCCTCAGATAGCAACTGCTGCACAAATACCAACTGCTGCCTTGGCCTCAGGCTTGGCCCCAACTGTACCTCAGACACAAGACACGTTCCCGCAGCATGTGCTCAGTCAGCCTCAACAAGTCAGAAAGGTTTTCACAAACTCAGCATCATCAAATACAGTTCTTCCATATCAGAGACATCCTGCCCCAGCTGTCCAGCAGCCCTTTATCAATAAGGCCTCCAACAGTGTTCTTCAATCCAGAAATGCTCCGCTTCCATCCCTGCAAAATGGACCTAACACACCCAACAAGCCTAGTTCACCCCCGCCACCCCAGCAATTTGTCGTCCAGCACTCTCTATTTGGGAGTCCAGTCGCCAAGACAAAAGATCCCCCCCGCTATGAGGAGGCCATCAAGCAGACACGCAGCACACAGGCCCCTCTGCCAGAGATTTCCAACGCTCACAGTCAGCAGATGGATGACCTCTTTGATATCCTCATTAAGAGTGGAGAGATCTCCCTCCCCATAAAAGAAGAACCTTCTCCTATTTCCAAAATGAGACCAGTGACAGCCAGCATCACCACAATGCCAGTGAATACAGTGGTGTCCCGGCCACCACCCCAAGTCCAAATGGCACCACCTGTATCTTTAGAACCTATGGGCAGTTTATCTGCCAGCTTAGAGAACCAACTAGAAGCTTTCTTGGATGGAACTTTACCCTCAGCCAATGAAATTCCTCCACTACAA"
            "AGCAGCAGTGAAGACAGAGAGCCCTTCTCTCTGATCGAGGACCTCCAGAATGATCTGCTGAGTCACTCAGGTATGCTGGACCATTCACACTCACCCATGGAGACTTCCGAGACCCAGTTTGCTGCAGGTACTCCCTGTCTGTCTCTCGACCTGTCAGACTCAAACTTGGACAACATGGAGTGGTTGGACATTACCATGCCCAACTCCTCTTCAGGACTCACTCCTCTCAGCACCACCGCGCCGAGCATGTTCTCTGCTGACTTTCTAGACCCACAGGACCTACCGCTGCCATGGGAC")

    seq2 = ("ATGGATCACACAGGGGCGATAGACACCGAGGATGAAGTGGGACCTTTAGCCCATCTTGCTCCGAGTCCTCAGAGTGAAGCTGTGGCTCATGAATTCCAGGAACTCTCCTTGCAGTCCAGTCAGAACTTACCCCCTCTGAACGAAAGGAAAAATGTGCTCCAGCTGAGGCTGCAACAAAGGAGGACGAGAGAACAACTAGTGGACCAGGGCATCATGCCACCTTTGAAGAGCCCAGCGGCATTCCATGAACAGATAAAAAGCTTGGAACGAGCCAGAACCGAAAACTTTTTGAAACACAAGATTCGGAGTCGACCAGATCGTTCTGAACTTGTCAGGATGCACATTTTAGAAGAAACATTTGCAGAGCCATCCCTGCAGGCTACTCAGATGAAGTTGAAAAGAGCTCGACTAGCAGATGATCTGAATGAAAAGATTGCTCAAAGACCTGGCCCTATGGAGCTGGTAGAGAAAAACATCCTTCCTGTGGACTCCAGTGTTAAAGAAGCAATTATAGGCGTTGGGAAGGAGGACTATCCCCACACTCAGGGCGATTTCTCATTTGATGAAGACAGCAGTGACGCTTTGTCTCCGGACCAGCCTGCGAGTCAGGAGTCACAGGGGTCAGCCGCGTCCCCAAGTGAGCCAAAAGTTAGTGAATCGCCATCTCCTGTGACTACAAACACTCCAGCCCAGTTTGCTTCAGTGTCCCCAACAGTTCCTGAATTCTTGAAAACTCCTCCAACTGCAGATCAGCCTCCCCCTCGGCCTGCAGCTCCTGTCCTCCCCACAAACACTGTGTCCTCAGCAAAGTCTGGCCCAGCGCTGGTGAAGCAAAGCCATCCCAAGAATCCAAATGACAAACACCGTAGCAAAAAGTGCAAAGATCCCAAACCACGGGTAAAGAAGTTAAAATACCACCAATACATTCCACCAGATCAGAAGGGTGAGAAGAACGAGCCGCAGATGGACTCCAACTACGCCCGCCTGCTCCAGCAGCAGCA"
            "GCTGTTCCTACAGCTGCAGATCCTGAGTCAGCAGAAGCAGCACTACAACTACCAGACCATCCTGCCTGCACCATTCAAGCCACTCAATGACAAAACTAGTAACAGTGGGAATTCAGCTTTGAACAATACCACACCTAACACACCAAGACAGAATACATCTGCTCCTGTGAGAAAGCCAGGACCTCTGCCTTCTAGCCTGGATGACTTAAAGGTGTCAGAACTGAAGACAGAACTGAAGTTAAGGGGTCTGCCAGTGTCAGGCACCAAACCAGACCTCATTGAGCGCCTGAAACCCTACCAGGAAGTGAACAGCAGCGGCCTTGCTGCTGGGGGCATCGTGGCAGTGTCATCGTCAGCCATTGTCACCAGTAACCCAGAAGTCACTGTGGCCTTGCCGGTTACAACACTACACAACACTGTGACTAGCCCAGTCTCTACTCTCAAGGCAGAATTGCCATCTACAGGAACCAGCAACGCAGCCCGTGTGGAAAATGTTCATTCCCCTCTGCCCATTTCACCATCTCCCTCTGAACAGTCCAGTCTCAGTACCGATGACACAAATATGGCAGACACTTTCACCGAGATTATGACCATGATGTCTCCTTCACAGTTCTTGAGTTCATCTCCTTTGAGAATGACAAATAATGAAGACAGTCTGAGTCCTACCAGCAGCACTCTGTCGAACCTGGAACTGGATGCAGCCGAGAAGGATCGCAAGCTTCAGGAGAAAGAGAAGCAAATCGAAGAACTGAAGAGGAAACTGGAACAAGAGCAGAAGCTCGTGGAAGTGCTGAAAATGCAACTTGAGGTTGAAAAACGAGGGCAGCAACAGCGGCCCCTGGAACCCCAGCCCAGTGCCCCAAGTCATTCTGTCAACAAGTCAGATCAGAAGCACAGCAGCCTTGGCTCCTCCATCAAAGACGAGGCCTCACTACCCGACTGCTCCAGCTCCAGGCAGCCCATCCCAGTAGCCAGCCACACTGTAGGCCAGCCTGTCTCTA"
            "CAGTTGGCCAGACCCTTGTTGCCAAAAAGGCTGTAGTTATCAAGCAAGAGGTCCCTGTGGCCCAGGCAGAGCAGCAGAGTGTCGTCTCGCAGTTTTATGTGAGTTCCCAGGGACAGCCACCGCCTGCTGTTGTTGCTCAGCCCCAGGCTTTACTGACCACGCAGACTGCTCAGCTCCTGCTCCCAGTGTCCATCCAGGGCTCGAGTGTCACCTCAGTGCAACTCCCTGTAGGCAGCCTCAAACTCCAGACTTCACCACAAGCAGGAATGCAGACTCAGCCTCAGATAGCAACTGCTGCACAAATACCAACTGCTGCCTTGGCCTCAGGCTTGGCCCCAGCTGTACCTCAGACACAAGACACGTTCCCACAGCATGTGCTCAGTCAGCCTCAACAAGTCAGAAAGGTTTTCACAAACTCAGCA---TCAAATACAGTTCTTCCATATCAGAGACATCCTGCTCCAGCTGTCCAGCAGCCCTTTATCAATAAGGCCTCCAACAGCGTTCTTCAATCCAGAAATGCTCCGCTTCCATCCCTGCAAAATGGACCTAACACACCCAACAAGCCTAGTTCACCCCCGCCACCCCAGCAATTTGTCGTCCAGCACTCTCTATTTGGGAGCCCAGTCGCCAAGACAAAAGATCCCCCCCGCTATGAGGAGGCCATCAAGCAAACACGCAGCACACAGGCCCCCCTGCCAGAGATTTCCAACGCACACAGTCAGCAGATGGATGACCTCTTTGATATACTCATTAAGAGTGGAGAGATCTCCCTCCCCATAAAAGAAGAACCTTCTCCTATTTCCAAAATGAGACCAGTGACAGCCAGCATCACCACAATGCCAGTGAATACAGTGGTGTCCCGGCCACCACCCCAAGTCCAAATGGCACCACCTGTATCTTTAGAACCTATGGGCAGTTTATCTGCCAGCTTAGAGAACCAACTAGAAGCTTTCTTGGATGGAACTTTACCCTCAGCCAATGAAATTGCTCCATTACAA"
            "AGCAGCAGTGAAGACAGAGAGCCCTTCTCTCTGATCGAGGACCTCCAGAATGACCTGCTGAGTCACTCAGGTATGCTGGACCATTCACACTCACCCATGGAGACTTCCGAGACCCAGTTTGCTGCAGGTACTCCCTGTCTGTCTCTCGACCTGTCAGACTCAAACTTGGATAACATGGAGTGGTTGGACATTACCATGCCCAACTCCTCTTCAGGACTCACTCCTCTCAGCACCACCACCCCAAGCATGTTCTCCGCCGACTTTCTAGACCCACAGGACCTACCACTACCATGGGAC")

    result = nei_gojobori(seq1, seq2)
    L, S1, S2, Nobs, Sobs = result
    mean_S = (S1+S2)/2
    pN, pS = Nobs / (L - mean_S), Sobs / mean_S
    dN, dS = dNdS(*result)
    # Values given by 'codeml'
    assert round(dN, 4) == 0.0065 and round(dS, 4) == 0.0607, "pN=%g, pS=%g, dN=%g, dS=%g" % (pN, pS, dN, dS)


if __name__ == '__main__':
    main()
