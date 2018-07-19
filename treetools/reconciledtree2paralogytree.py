#!/usr/bin/env python3
# -*- coding: utf-8 -*-


""" Take a 'reconciled' gene tree and converts it to a 'paralogy tree':

- reconciled tree: gene tree with species mapped at nodes, i.e gene tree with
                   gene duplication information.

- paralogy tree: tree whose branch represent a **pair of paralogs**.
                 The aim is to make the display of duplication/deletion events easier.
"""


from sys import stdin
import re
import argparse
import ete3

from LibsDyogen import myPhylTree
from genetree_drawer import get_taxon, get_taxon_treebest, infer_gene_event

ENSEMBL_VERSION = 85
PHYLTREEFILE = "/users/ldog/glouvel/ws_alouis/GENOMICUS_SVN/data{0}/" \
                   "PhylTree.Ensembl.{0}.conf"


def buildparalogies(genetree, get_taxon, ancgene2sp, ensembl_version=ENSEMBL_VERSION):
    
    # List of trees representing paralogies
    paralogies = []

    def extendparalogy(paralogs, taxon, paralogy_node=None):
        
        # set of paralogs for each child speciation
        paralogs_at_speciation = {}
        
        # Empty paralogs if genes reached a speciation node.
        # Otherwise, replace the node by the duplication descendants and start a new paralogy
        while paralogs:
            paralog = paralogs.pop()

            children_taxa = [get_taxon(ch, ancgene2sp, ensembl_version) \
                             for ch in paralog.children]
            event = infer_gene_event(paralog, taxon, set(children_taxa))

            if event == 'dup':
                new_paralogy = ete3.TreeNode(name='-'.join(ch.name for ch in paralog.children))
                new_paralogy.add_feature("taxon", taxon)
                if paralogy_node is None:
                    # Root paralogy
                    paralogies.append(new_paralogy)
                else:
                    # Duplicated paralogy (stem from a set of paralogs)
                    paralogy_node.add_child(new_paralogy)

                #assert len(paralog.children) > 1
                extendparalogy(set(paralog.children), taxon, new_paralogy)

                paralogs.update(paralog.children)

            else:
                for child_taxon, speciated_paralog in zip(children_taxa, paralog.children):
                    child_taxon_paralogs = paralogs_at_speciation.setdefault(taxon, set())
                    child_taxon_paralogs.update(speciated_paralog)

        assert taxon not in paralogs_at_speciation, "Intermediate speciation nodes are missing."

        for child_taxon, speciated_paralogs in paralogs_at_speciation.item():
            # Speciated paralogy
            if len(speciated_paralogs) > 1:
                new_paralogy = ete3.TreeNode(name='-'.join(ch.name for ch in speciated_paralogs))
                new_paralogy.add_feature("taxon", child_taxon)
                if paralogy_node is None:
                    paralogies.append(new_paralogy)
                else:
                    paralogy_node.add_child(new_paralogy)
            else:
                new_paralogy = None
            extendparalogy(speciated_paralogs, child_taxon, new_paralogy)

    
    taxon = get_taxon(genetree, ancgene2sp, ensembl_version)
    extendparalogy(set((genetree,)), taxon)
    return paralogies



def main(inputnwk, outputnwk, ensembl_version=ENSEMBL_VERSION,
         phyltreefile=PHYLTREEFILE, treebest=False):

    phyltree = myPhylTree.PhylogeneticTree(phyltreefile.format(ensembl_version))

    ancgene2sp = re.compile(r'(' + r'root|'
                    + r'|'.join(list(phyltree.listSpecies) +
                                sorted(phyltree.listAncestr,
                                       key=lambda a: len(a),
                                       reverse=True)).replace(' ','\.')
                    + r')(.*)$')

    genetree = ete3.Tree(inputnwk or stdin.read(), format=1)

    if treebest: get_taxon = get_taxon_treebest
    
    for paralogy in buildparalogies(genetree, get_taxon, ancgene2sp, ensembl_version):
        out = paralogy.write(format=1, format_root_node=True, features=['taxon'], outfile=outputnwk)
        if out:
            print(out)

    

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('inputnwk', nargs='?')
    parser.add_argument('outputnwk', nargs='?')
    parser.add_argument('-e', '--ensembl-version', type=int, default=ENSEMBL_VERSION)
    parser.add_argument('-p', '--phyltreefile', default=PHYLTREEFILE)
    parser.add_argument('-t', '--treebest', action='store_true')
    
    args = parser.parse_args()
    main(**vars(args))
