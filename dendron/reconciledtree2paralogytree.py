#!/usr/bin/env python3
# -*- coding: utf-8 -*-


""" Take a 'reconciled' gene tree and converts it to a 'paralogy tree':

- reconciled tree: gene tree with species mapped at nodes, i.e gene tree with
                   gene duplication information.

- paralogy tree: tree whose branch represent a **pair of paralogs**.
                 The aim is to make the display of duplication/deletion events easier.
"""


from sys import stdin, stderr
import re
import argparse
import ete3

from LibsDyogen import myPhylTree
from genetree_drawer import get_taxon, get_taxon_treebest, infer_gene_event
#from dendron.reconciled

ENSEMBL_VERSION = 85
PHYLTREEFILE = "/users/ldog/glouvel/GENOMICUS{0}/PhylTree.Ensembl.{0}.conf"


def pdbug(indent, *args, prefix='', sep=' ', end='\n', file=stderr):
    """Debugging print"""
    message = sep.join(str(arg) for arg in args)
    message = '    '*indent + prefix + message.replace('\n', '\n' + '    '*indent + prefix)
    print(message, end=end, file=file)


def buildparalogies(genetree, get_taxon, ancgene2sp,
                    ensembl_version=ENSEMBL_VERSION,
                    include_singleton_branches=False):

    callnb = 0  #debug
    indent = -1  #debug

    # List of trees representing paralogies
    paralogies = []

    if include_singleton_branches:
        def make_singleton_node(nodename, taxon):
            snode = ete3.TreeNode(name=nodename)
            snode.add_feature("S", taxon)
            snode.add_feature("P", False)  # whether it is a paralogy branch?
            return snode
    else:
        def make_singleton_node(nodename, taxon):
            return None

    # Paralog_packs: List of the 2 packs of genes on each 2 sides of the dup.
    # Can be of length>2 if the duplication is such.
    def extendparalogy(paralog_packs: list, taxon: str,
                       parent_paralogy=None) -> None:
        """'Speciate' the paralogy, i.e divide the `paralog_packs` into their
        descendant paralog_packs in each child_taxon, and call itself again on
        each descendant paralogy"""
        nonlocal callnb, indent  #debug
        callnb += 1  #debug
        indent += 1  #debug

        pdbug(indent, taxon, 'call:', callnb, prefix='# ')
        pdbug(indent, 'paralog_packs:', paralog_packs)

        # Create a node object representing the paralogy, if any. Attach to parent if any.
        para_size = len(paralog_packs)
        if para_size > 1:
            para_name = ';'.join('-'.join(ch.name for ch in pack) for pack in paralog_packs)
            new_paralogy = ete3.TreeNode(name=para_name)
            new_paralogy.add_feature("S", taxon)
            new_paralogy.add_feature("P", True)
            #new_paralogy.add_feature("D", ("Y" if getattr(parent_paralogy, 'S', None)==taxon else "N")) #optional
            #parent_paralogy.add_feature("D", )
        else:
            assert len(paralog_packs[0]) == 1
            new_paralogy = make_singleton_node(list(paralog_packs[0])[0].name, taxon)
            pdbug(indent, "No paralogy")
        
        if new_paralogy is not None:
            if parent_paralogy is None:
                paralogies.append(new_paralogy)
                pdbug(indent, 'Create new paralogy', new_paralogy.name)
            else:
                parent_paralogy.add_child(new_paralogy)
                pdbug(indent, ("dup" if parent_paralogy.D == "Y" else "spe") + '-extend paralogy from', parent_paralogy.name)
        #pdbug(indent, new_paralogy.get_tree_root() if new_paralogy else None, prefix='> ')

        # The descendant paralog_pack in each species after the speciation:
        paralog_packs_after_speciation = {}  # {ch: [set()] * para_size}
        #speciated_taxa = set()
        
        # Empty paralogs if genes reached a speciation node.
        # Otherwise, replace the node by the duplication descendants and start a new paralogy
        seen_paralogs = set() #debug

        for pack_i, paralog_pack in enumerate(paralog_packs):
            has_sub_paralogies = len(paralog_pack) > 1
            #pdbug(indent, pack_i, paralog_pack, prefix='* ')

            while paralog_pack:
                paralog = paralog_pack.pop()  # if isinstance(paralog_pack, set) else paralog = paralog_pack

                children_taxa = [get_taxon(ch, ancgene2sp, ensembl_version)
                                 for ch in paralog.children]
                event = infer_gene_event(paralog, taxon, set(children_taxa))

                #pdbug(indent, repr(paralog), event+':', paralog.children, prefix=' ** ')

                if event == 'dup':
                    #if para_size > 1 / new_paralogy is not None
                    if has_sub_paralogies:
                        paralog_pack.update(paralog.children)

                    pdbug(indent, 'Dup!')
                    if new_paralogy is not None:
                        new_paralogy.add_feature('D', 'Y')
                    extendparalogy([set((ch,)) for ch in paralog.children],
                                   children_taxa[0], new_paralogy)
                    indent -= 1  #debug
                    #also extendparalogy of the current paralogy (new_paralogy)

                else:
                    for child_taxon, speciated_paralog in zip(children_taxa, paralog.children):
                        speciated_paralog_packs = paralog_packs_after_speciation.setdefault(child_taxon, [set() for p in range(para_size)])
                        #if child_taxon not in paralog_packs_after_speciation:
                        #    paralog_packs_after_speciation[child_taxon] = [set() for i * para_size
                        speciated_paralog_packs[pack_i].add(speciated_paralog)
                        #pdbug(indent, ' ** Add speciated paralog:', repr(speciated_paralog),
                        #      '->', pack_i, child_taxon)

                        assert speciated_paralog not in seen_paralogs #debug
                        seen_paralogs.add(speciated_paralog) #debug

        assert taxon not in paralog_packs_after_speciation, \
            "Intermediate speciation nodes are missing at: %s -> %s" % \
            (taxon, tuple(paralogs_after_speciation.keys()))

        #if para_size > 1 or no dup:
        #if there hasn't been a complete update of the paralog_pack
        # if there was a dup, this is redundant.
        for child_taxon, speciated_paralog_packs in paralog_packs_after_speciation.items():
            speciated_paralog_packs = [pack for pack in speciated_paralog_packs if pack]

            redundant_nodes = (len(speciated_paralog_packs) > 1 and set.intersection(*speciated_paralog_packs))
            assert not redundant_nodes, paralog_packs_after_speciation

            pdbug(indent, 'Spe!')
            if new_paralogy is not None:
                new_paralogy.add_feature('D', 'N')
            extendparalogy(speciated_paralog_packs, child_taxon, new_paralogy)
            indent -= 1  #debug

    taxon = get_taxon(genetree, ancgene2sp, ensembl_version)
    extendparalogy([set((genetree,))], taxon)
    return paralogies



def main(inputnwk, outputnwk, ensembl_version=ENSEMBL_VERSION,
         phyltreefile=PHYLTREEFILE, treebest=False,
         include_singleton_branches=False):

    phyltree = myPhylTree.PhylogeneticTree(phyltreefile.format(ensembl_version))

    ancgene2sp = re.compile(r'(' + r'root|'
                    + r'|'.join(list(phyltree.listSpecies) +
                                sorted(phyltree.listAncestr,
                                       key=lambda a: len(a),
                                       reverse=True)).replace(' ','\.')
                    + r')(.*)$')

    genetree = ete3.Tree(inputnwk or stdin.read(), format=1)

    #get_taxon = get_taxon_treebest if treebest else get_taxon

    for paralogy in buildparalogies(genetree, get_taxon, ancgene2sp,
                                    ensembl_version,
                                    include_singleton_branches):
        out = paralogy.write(format=1, format_root_node=True, features=['S', 'D', 'P'], outfile=outputnwk)
        if out:
            print(out)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('inputnwk', nargs='?')
    parser.add_argument('outputnwk', nargs='?')
    parser.add_argument('-e', '--ensembl-version', type=int, default=ENSEMBL_VERSION)
    parser.add_argument('-p', '--phyltreefile', default=PHYLTREEFILE)
    parser.add_argument('-t', '--treebest', action='store_true')
    parser.add_argument('-s', '--include-singleton-branches',
                        action='store_true',
                        help='Include single gene nodes, but with a special mark.')

    args = parser.parse_args()
    main(**vars(args))
