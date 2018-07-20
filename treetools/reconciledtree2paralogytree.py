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
PHYLTREEFILE = "/users/ldog/glouvel/GENOMICUS{0}/PhylTree.Ensembl.{0}.conf"


def buildparalogies(genetree, get_taxon, ancgene2sp, ensembl_version=ENSEMBL_VERSION):
    
    # List of trees representing paralogies
    paralogies = []
    callnb = 0

    # Paralog_packs
    def extendparalogy(paralogs, taxon, paralogy_node=None, indent=0):
        nonlocal callnb
        callnb += 1

        # set of paralogs for each child speciation
        paralogs_at_speciation = {}
        
        # Empty paralogs if genes reached a speciation node.
        # Otherwise, replace the node by the duplication descendants and start a new paralogy
        print('    '*indent + '## Call: ', callnb)
        print('    '*indent + taxon, 'paralogs:', [p.name for p in paralogs],
                'node:', paralogy_node is not None)
        #for paralogs in paralog_packs:
        while paralogs:
            paralog = paralogs.pop() # if isinstance(paralog, tuple)

            children_taxa = [get_taxon(ch, ancgene2sp, ensembl_version) \
                             for ch in paralog.children]
            event = infer_gene_event(paralog, taxon, set(children_taxa))

            if event == 'dup':
                print('    '*indent + '- Dup;', end=' ')
                new_paralogy = ete3.TreeNode(name='-'.join(ch.name for ch in paralog.children))
                new_paralogy.add_feature("S", taxon)
                new_paralogy.add_feature("D", "Y")
                paralogs.add(tuple(paralog.children))
                print('    '*indent + 'updated paralogs:', [p.name for p in paralogs])

                if paralogy_node is None:
                    # Root paralogy
                    paralogies.append(new_paralogy)
                    print('    '*indent + 'Create new paralogy', new_paralogy.name)
                else:
                    # Duplicated paralogy (stem from a set of paralogs)
                    paralogy_node.add_child(new_paralogy)
                    print('    '*indent + 'dup-extend paralogy from', paralogy_node.name)

                #assert len(paralog.children) > 1
                extendparalogy(set(paralog.children), taxon, new_paralogy, indent+1)
                print('    '*indent + '## Back to parent call after dup.')

            else:
                print('    '*indent + '- Spe/Leaf;', end=' ')
                print('    '*indent + paralog.name, '->', children_taxa, [chp.name for chp in paralog.children])

                for child_taxon, speciated_paralog in zip(children_taxa, paralog.children):
                    child_taxon_paralogs = paralogs_at_speciation.setdefault(child_taxon, set())
                    child_taxon_paralogs.add(speciated_paralog)
                    print('    '*indent + 'speciated_paralog:', child_taxon, speciated_paralog.name)

        assert taxon not in paralogs_at_speciation, \
            "Intermediate speciation nodes are missing at: %s -> %s" % \
            (taxon, tuple(paralogs_at_speciation.keys()))

        print('    '*indent + 'paralogs at speciation:\n' + '    '*indent + ((',\n'+'    '*indent).join('%r: %s' % (cht, [p.name for p in chp]) for cht,chp in paralogs_at_speciation.items()) if paralogs_at_speciation else '{}'))

        #for child_taxon, speciated_paralogs in paralogs_at_speciation.items():
        while paralogs_at_speciation:
            child_taxon, speciated_paralogs = paralogs_at_speciation.popitem()
            # Speciated paralogy
            if len(speciated_paralogs) > 1:
                new_paralogy = ete3.TreeNode(name='-'.join(ch.name for ch in speciated_paralogs))
                new_paralogy.add_feature("S", child_taxon)
                new_paralogy.add_feature("D", "N")
                if paralogy_node is None:
                    # Means that there already was a creation from dup. Pass.
                    #import ipdb; ipdb.set_trace(context=1)
                    paralogies.append(new_paralogy)
                    print('    '*indent + 'Create new paralogy', new_paralogy.name)
                else:
                #if paralogy_node is not None:
                    paralogy_node.add_child(new_paralogy)
                    print('    '*indent + 'speciation-extend paralogy from', paralogy_node.name)
            else:
                print('    '*indent + 'Continue with no paralogy with', [p.name for p in speciated_paralogs])
                new_paralogy = None
            extendparalogy(speciated_paralogs, child_taxon, new_paralogy, indent+1)
            print('    '*indent + '## Back to parent call after spe.')

    
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

    #get_taxon = get_taxon_treebest if treebest else get_taxon
    
    for paralogy in buildparalogies(genetree, get_taxon, ancgene2sp, ensembl_version):
        out = paralogy.write(format=1, format_root_node=True, features=['S', 'D'], outfile=outputnwk)
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
