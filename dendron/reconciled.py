#!/usr/bin/env python3
# -*- coding: utf-8 -*-


"""Tools to handle "reconciled" phylogenetic gene trees

Reconciliation: mapping of a gene tree on a species tree.

Uses `ete3.Tree` objects.

This tools recognize different ways of annotating the species at each node:
    - TreeBest: the `S` tag in NHX comments.
    - custom (Ensembl gene trees): the start of the node name, e.g
                                   HomoPanENSGT0039000...
"""


from genomicustools.identify import ultimate_seq2sp
import logging
#logging.basicConfig(format='%(levelname)s:%(module)s l.%(lineno)d:%(funcName)s:%(message)s')
logger = logging.getLogger(__name__)
#ch = logging.StreamHandler()
#ch.setFormatter(logging.Formatter("%(levelname)s:%(message)s"))
#logger.addHandler(ch)


# For types
#import re, ete3


ENSEMBL_VERSION = 85


#def get_taxon(node: ete3.TreeNode, ancgene2sp : re._pattern_type,
#              ensembl_version : int = ENSEMBL_VERSION):
def get_taxon(node, ancgene2sp, ensembl_version=ENSEMBL_VERSION):
    """from a gene name in my newick gene trees, find the taxon:
        either:
            - node is a leaf (e.g ENSMUSG00...)
            - node is internal (e.g Mus.musculusENSGT...)"""
    if node.is_leaf():
        taxon = ultimate_seq2sp(node.name, ensembl_version)
    else:
        try:
            taxon = ancgene2sp.match(node.name).group(1).replace('.', ' ')
        except AttributeError:
            raise ValueError("Can not match species name in %r" % node.name)
    return taxon


def get_taxon_treebest(node, *args):
    """Get the taxon of a gene node, using a treebest output tree.
    
    *args are not used, they are here for compatibility with `get_taxon`"""
    try:
        return node.S.replace('.', ' ')
    except AttributeError as err:
        err.args += (node.name, node)
        raise


def infer_gene_event(node, taxon, children_taxa):
    """Use taxon information to tell whether a gene tree node (Ete3 format) is:
    - a leaf,
    - a speciation,
    - or a duplication.
    
    Use TreeBest annotation if present.
    
    param: `children_taxa` must be a set (because the number of *uniq* elements
           is used).
    
    This is the original version from `genetree_drawer.py` (2018/11/19).
    """
    #if node.is_leaf():
    #import ete3
    #assert isinstance(node, ete3.TreeNode)
    #assert isinstance(taxon, str)
    #assert isinstance(children_taxa, set)

    if not children_taxa:
        return 'leaf'

    children_names = [ch.name for ch in node.children]
    treebest_isdup = getattr(node, 'D', None)

    if (treebest_isdup == 'Y') or taxon in children_taxa:
       #(len(children_taxa) == 1 and taxon in children_taxa): #\
                    #and len(node.children) > 1 \
       #or \
       #(len(node.children) == 1 and taxon in children_taxa):
        # Should rather check if taxon in children_taxa.
        if treebest_isdup == 'N':
            logger.warning("The node %r -> %s "
                  "is marked as a *speciation* but looks like a duplication"
                  " after taxon information.",
                  node.name, children_names)
        elif treebest_isdup == 'Y' and taxon not in children_taxa:
            logger.warning("The node %r -> %s "
                    "is marked as a *duplication* but is ambiguous: "
                    "intermediate speciation nodes are missing.",
                  node.name, children_names)

        if len(children_taxa) > 1:  # or treebest_isdup == 'Y'
            logger.warning("The node %r -> %s "
                  "is a duplication + a speciation. Not "
                  "truly reconciled tree.",
                  node.name, children_names)
        elif len(node.children) == 1:
            logger.warning("The node %r -> %s "
                  " is a duplication/speciation with one "
                  "descendant.",
                  node.name, children_names)
        return 'dup'
    else:
        ### FIXME
        if len(children_taxa) == 1:
            msg = "The node %r -> %s"
            if len(node.children) > 1:
                msg += " is a duplication + a speciation. Not " \
                       "truly reconciled tree."
                #return 'dup'
            #else:
            #    msg += " is a duplication/speciation with one " \
            #           "descendant."
                #if not node.is_root(): event = 'dup'
                # TODO: check this in the duplication block above.
                logger.warning(msg, node.name, children_names)

        return 'spe'

# TODO:
# infer_gene_event_treebest
# infer_gene_event_taxa
# infer_gene_events that warns if mismatches between the 2.

def get_children_ete3(node, *args):
    return node.children

def get_children_dyogenprot(node, tree):
    """return children from a node in a LibsDyogen.myProteinTree format."""
    return [ch for ch,_ in tree.data.get(node, [])]

def get_node_name_ete3(node, *args):
    return node.name

def get_node_name_dyogenprot(node, tree):
    nodeinfo = tree.info[node]
    return nodeinfo.get('gene_name',
                        (nodeinfo['taxon_name'] +
                         nodeinfo.get('family_name', '').split('/')[0]))

def infer_gene_event_taxa(node, taxon, children_taxa,
                          get_children=get_children_ete3,
                          get_name=get_node_name_ete3, *args):

    """Use taxon information to tell whether a gene tree node is:
    - a leaf,
    - a speciation,
    - or a duplication.
    
    param: `children_taxa` must be a set (because the number of *uniq* elements
           is used).
    param: *args: extra arguments to be passed to `get_children`/`get_node_name`.
    """
    ### "ambiguous" or "dupspe" value should be returned in case of doubt.

    if not children_taxa:
        return 'leaf'

    children = get_children(node, *args)
    nodename = get_name(node, *args)
    
    if taxon in children_taxa:
       #(len(children_taxa) == 1 and taxon in children_taxa): #\
                    #and len(children) > 1 \
       #or \
       #(len(children) == 1 and taxon in children_taxa):
        # Should rather check if taxon in children_taxa.
        if len(children_taxa) > 1:
            logger.warning("The node %r -> %s "
                  "is a duplication + a speciation. Not "
                  "truly reconciled tree." %
                  (nodename, [get_name(ch, *args) for ch in children]))
        elif len(children) == 1:
            logger.warning("The node %r -> %s "
                  " is a duplication/speciation with one "
                  "descendant." %
                  (nodename, [get_name(ch, *args) for ch in children]))
        return 'dup'
    else:
        ### FIXME
        if len(children_taxa) == 1:
            msg = "The node %r -> %s"
            if len(children) > 1:
                msg += " is a duplication + a speciation. Not " \
                       "truly reconciled tree."
                #return 'dup'
            #else:
            #    msg += " is a duplication/speciation with one " \
            #           "descendant."
                #if not node.is_root(): event = 'dup'
                # TODO: check this in the duplication block above.
                logger.warning(msg, nodename, [get_name(ch, *args) for ch in children])

        return 'spe'
