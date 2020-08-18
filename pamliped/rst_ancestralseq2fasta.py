#!/usr/bin/env python3
# -*- coding: utf-8 -*-


"""Convert the ancestral sequences output of `codeml` (PAML) into a fasta file."""

import re


def get_seq(lines):
    line = next(lines)
    labtree_re = re.compile(r"tree with node labels for Rod Page's TreeView$")
    while not labtree_re.match(line):
        line = next(lines)

    line = next(lines)

    labtree = line

    line = next(lines)

    anc_node_re = re.compile(r'Nodes (\d+) to (\d+) are ancestral$')
    while not anc_node_re.match(line):
        line = next(lines)

    anc_nodes = [int(x) for x in anc_node_re.match(line).groups()]

    line = next(lines)

    anc_seqs = []
    codon_prob_re = re.compile('([ACGT]{3})\(([01]\.\d+)\)')
    start_anc_re = re.compile(r'Prob distribution at node (\d+), by site')

    for node in range(anc_nodes[0], anc_nodes[1]):
        while not start_anc_re.match(line):
            line = next(lines)

        node_num = int(start_anc_re.match(line).group(1))

        line = next(lines)
        
        while not re.match('^\s+site\s+Freq\s+Data$', line):
            line = next(lines)

        line = next(lines)
        line = next(lines)

        anc_seq = ''
        while re.match(r'\s+\d+\s+\d+', line):
            codons_str = line.split(':')[1].split()
            codons = []
            for cod_str in codons_str:
                cod, prob = codon_prob_re.match(cod_str).groups()
                codons.append((cod, float(prob)))

            anc_seq += max(codons, key=lambda x: x[1])[0]

        anc_seqs.append((node_num, anc_seq))

    return anc_seqs




