#!/usr/bin/env python3
# -*- coding: utf-8 -*-


"""Phylogenetic trees can be represented as string in various formats:
    - Newick (many variants) -> Ete3 is a good tool for this.
    - Nexus
    - PhyloXML
      
Also the embedding of comments may vary.
- NHX in newick trees
- another format in nexus trees such as produced by BEAST.

The Bio.Phylo.convert doesn't properly convert the comments from BEAST Nexus files.
"""

import re
from collections import OrderedDict
from Bio import Phylo


def beast_comment_parser(text):
    assert not text or text.startswith('[&') and text.endswith(']')
    variables = OrderedDict()
    #reg_interval = re.compile('^{-?\d+')
    if text is None:
        return variables

    text = re.sub(r'^\[&|\]$', '', text)
    while text:
        try:
            var, text = text.split('=', maxsplit=1)
        except ValueError as err:
            err.args += ('from "%s"' % text,)
            raise

        if text.startswith('{'):
            end = text.find('}')
            value = [float(x) for x in text[1:end].split(',')]
            text = text[(end+1):].lstrip(',')
        else:
            try:
                value, text = text.split(',', 1)
            except ValueError:
                value = text
                text = ''
            try:
                value = float(value)
            except ValueError:
                pass
        variables[var] = value

    return variables


def NHX_comment_formatter(variables: dict):
    if not variables:
        return ''
    text = '&&NHX:'  # Bio.Phylo does not want '[' or ']' at the start/end.
    for var, value in variables.items():
        text += var + '='
        if isinstance(value, (list, tuple)):
            text += '{' + '..'.join(str(x) for x in value) + '}'
        else:
            text += str(value)
        text += ':'
    
    return text.rstrip(':')


def nexus2nhx(infile, outfile):
    """Designed to work with Beast/treeannotator output"""
    with open(outfile, 'w') as out:
        for tree in Phylo.parse(infile, 'nexus'):
            for node in tree.get_terminals() + tree.get_nonterminals():
                node.comment = NHX_comment_formatter(
                                    beast_comment_parser(node.comment))
            Phylo.write(tree, out, 'newick')

## TODO: tests:
# 1. Check if the resulting newick can be parsed by Ete3/Biopython/Newick-utils.

