#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Module centralizing multiple phylogenetic tree parsers"""

from sys import stdin
import os.path as op


def read_multinewick(lines, stripchars='\r\n'):
    newick = ''
    for line in lines:
        parts = line.strip(stripchars).split(';')
        newick += parts.pop(0)
        while parts:
            yield newick + ';'
            newick = parts.pop(0)
    if newick.strip(stripchars):
        yield newick + ';'


def iter_from_phyltree(treefile, *args, **kwargs):
    from LibsDyogen import myPhylTree
    yield myPhylTree.PhylogeneticTree(treefile, *args, **kwargs)

def iter_from_prottree(treefile, *args, **kwargs):
    from LibsDyogen import myProteinTree
    for tree in myProteinTree.loadTree(treefile, *args, **kwargs):
        yield tree

def iter_from_ete3(treefile, *args, **kwargs):
    import ete3
    if treefile is not stdin and not op.exists(treefile):
        yield ete3.Tree(treefile, *args, **kwargs)
        raise StopIteration

    f = treefile if treefile is stdin else open(treefile)
    try:
        for newick in read_multinewick(f):
            try:
                yield ete3.Tree(newick, *args, **kwargs)
            except ete3.parser.newick.NewickError as err:
                err.args = (err.args[0] + 'ERROR with treefile %r ...' % treefile[:50],)
                raise
    finally:
        if f is not stdin:
            f.close()

def iter_from_biophylo(treefile, *args, **kwargs):
    from Bio import Phylo
    for tree in Phylo.parse(treefile, *args, **kwargs):
        yield tree

def iter_from_skbio(treefile, *args, **kwargs):
    import skbio
    raise NotImplementedError()


parserchoice = {'phyltree': iter_from_phyltree,
                'PhylTree': iter_from_phyltree,
                'Phyltree': iter_from_phyltree,
                'prottree': iter_from_prottree,
                'ProtTree': iter_from_prottree,
                'Prottree': iter_from_prottree,
                'ete3':     iter_from_ete3,
                'Ete3':     iter_from_ete3,
                'skbio':    iter_from_skbio,
                'biophylo': iter_from_biophylo}
