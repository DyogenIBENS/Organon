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

from Bio.Nexus import Nexus
from Bio.Phylo import NewickIO
from dendro.bates import dfw_pairs_generalized
from dendro.any import BioNexus as bionexus_methods
from dendro.converters import BioNexusTrees_to_BioPhylo



def parse_newick2nexus(handle):
    nx = Nexus.Nexus()
    for tree in Phylo.parse(handle, 'newick'):
        nx.trees.append(Nexus.Tree(tree.format('newick')))  # Sorry, I'm relying on serialization into a newick string.
    return nx

def index_nexus_treenodes(nx):
    # See:
    #   nxtree.all_ids()
    #   nxtree.get_taxa()
    #   nxtree.id
    translate = {}
    for ntree in nx.trees:
        for node_id, node in ntree.chain.items():
            nodename = node.data.taxon
            if node_id in translate:
                raise ValueError('Duplicate node id found: %s %s' % (node_id, nodename))
            translate[node_id] = nodename
    nx.translate = translate  # Overwrite current translate

# Rewrite of Bio.Phylo.NexusIO:
# Structure of a Nexus BayesTraits file
NEX_TEMPLATE = """\
#NEXUS
Begin Trees;
 translate
  %(translate)s;
  %(trees)s
End;
"""
# I don't use this block
NEX_TAXA_BLOCK = """Begin Taxa;
 Dimensions NTax=%(count)d;
 TaxLabels %(labels)s;
End;
"""

# 'index' starts from 1; 'tree' is the Newick tree string
TREE_TEMPLATE = "Tree tree%(index)d=%(tree)s"

def write_nexus_trees_to_bayestraits(nx, handle, **kwargs):
    """Modified from Bio.Phylo.NexusIO.write():

    add a translate block converting leaf names to integers.
    """
    trees = BioNexusTrees_to_BioPhylo(nx.trees)
    writer = NewickIO.Writer(trees)
    nexus_trees = [
        TREE_TEMPLATE % {"index": idx + 1, "tree": nwk}
        for idx, nwk in enumerate(
            writer.to_strings(plain=False, plain_newick=True, **kwargs)
        )
    ]
    translate = ["%d %s" % id_name for id_name in nx.translate.items()]

    # Unused in my output format (BayesTraits) + why aren't they unique?
    tax_labels = [taxon for nt in nx.trees for taxon in nt.get_taxa()]
    #tax_labels = [str(x.name) for x in chain(*(t.get_terminals() for t in trees))]

    text = NEX_TEMPLATE % {
        "count": len(tax_labels),
        "labels": " ".join(tax_labels),
        "trees": "\n".join(nexus_trees),
        "translate": ",\n ".join(translate)
    }
    handle.write(text)
    return len(nexus_trees)


def nhx2bayestraits(infile, outfile, **nwk_kwargs):
    """
    Phylo.convert(infile, 'newick', outfile, 'nexus') does not reindex labels,
    as required by some programs (BayesTraits).
    """
    nx = parse_newick2nexus(infile)
    index_nexus_treenodes(nx)
    with open(outfile, 'w') as out:
        write_nexus_trees_to_bayestraits(nx, out, **nwk_kwargs)


# Checks:
# 1. load the resulting nexus
# nx2 = Nexus.Nexus(nexfile)
# 2. plot the converted tree:
# plottree(nx2.trees[0], bionexus_methods.get_items, bionexus_methods.get_label, nx2.trees[0].node(nx2.trees[0].root))
