#!/usr/bin/env python3
# -*- coding: utf-8 -*-


""" Take a 'reconciled' gene tree and converts it to a 'paralogy tree':

- reconciled tree: gene tree with species mapped at nodes, i.e gene tree with
                   gene duplication information.

- paralogy tree: tree whose branch represent a **pair of paralogs**.
                 The aim is to make the display of duplication/deletion events easier.
"""


from sys import stdin, stderr, stdout
import re
import argparse
import ete3

from LibsDyogen import myPhylTree
from genetree_drawer import get_taxon, get_taxon_treebest, infer_gene_event

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

        # Create a node object representing the paralogy, if any. Attach to parent if any.
        para_size = len(paralog_packs)
        if para_size > 1:
            para_name = '|'.join('-'.join(ch.name for ch in pack) for pack in paralog_packs)
            current_paralogy = ete3.TreeNode(name=para_name)
            current_paralogy.add_feature("S", taxon)
            current_paralogy.add_feature("P", True)
            #current_paralogy.add_feature("D", ("Y" if getattr(parent_paralogy, 'S', None)==taxon else "N")) #optional
            #parent_paralogy.add_feature("D", )
        else:
            #assert len(paralog_packs[0]) == 1  # Not sure if it should happen
            current_paralogy = make_singleton_node(list(paralog_packs[0])[0].name, taxon)
            pdbug(indent, "No paralogy")
        
        if current_paralogy is not None:
            if parent_paralogy is None:
                current_paralogy.add_feature('D', 'Y')
                paralogies.append(current_paralogy)
                pdbug(indent, 'Create new paralogy', current_paralogy.name)
            else:
                current_paralogy.add_feature('sub', True)  # It is a "sub-paralogy"
                parent_paralogy.add_child(current_paralogy)
                pdbug(indent, ("dup" if parent_paralogy.D == "Y" else "spe") + '-extend paralogy from', parent_paralogy.name)
        #pdbug(indent, current_paralogy.get_tree_root() if current_paralogy else None, prefix='> ')
        pdbug(indent, 'paralog_packs:', paralog_packs)

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
                    #if para_size > 1 / current_paralogy is not None
                    pdbug(indent, 'Dup!')
                    if current_paralogy is not None:
                        current_paralogy.add_feature('D', 'Y')
                        #if current_paralogy.P:
                        #    # Only update the pack when we are in a paralogy.
                        #    # If not, the current branch can be dropped, because
                        #    # each paralog is already going to be checked in 
                        #    # subsequent calls.
                    paralog_pack.update(paralog.children)
                    #if callnb == 32: import ipdb; ipdb.set_trace(context=1)
                    if not has_sub_paralogies:
                        extendparalogy([set((ch,)) for ch in paralog.children],
                                       children_taxa[0], current_paralogy)
                        if not getattr(parent_paralogy, 'P', False):
                            # It miraculously worked. I don't know why.
                            # This is needed with option `include_singleton_branches=True`:
                            # it avoids drawing a "duplicate" branch for singleton genes left from a new paralogy node.
                            break

                        extendparalogy(paralog_packs, children_taxa[0], current_paralogy)
                        # Seems to work okay, but I suspect some paralogs are
                        # checked several times.
                        #break
                        has_sub_paralogies = True


                    indent -= 1  #debug
                    #also extendparalogy of the current paralogy (current_paralogy)
                    #if has_sub_paralogies:
                    #    paralog_pack.update(paralog.children)
                        #extendparalogy
                    #has_sub_paralogies = True


                else:
                    # What if it is a leaf? The paralog is just popped out.
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
            if current_paralogy is not None:
                current_paralogy.add_feature('D', 'N')
            extendparalogy(speciated_paralog_packs, child_taxon, current_paralogy)
            indent -= 1  #debug

    taxon = get_taxon(genetree, ancgene2sp, ensembl_version)
    extendparalogy([set((genetree,))], taxon)
    return paralogies



def main(inputnwk, outputnwk, ensembl_version=ENSEMBL_VERSION,
         phyltreefile=PHYLTREEFILE, treebest=False,
         include_singleton_branches=False, include_singleton_root=False):

    phyltree = myPhylTree.PhylogeneticTree(phyltreefile.format(ensembl_version))

    ancgene2sp = re.compile(r'(' + r'root|'
                    + r'|'.join(list(phyltree.listSpecies) +
                                sorted(phyltree.listAncestr,
                                       key=lambda a: len(a),
                                       reverse=True)).replace(' ','\.')
                    + r')(.*)$')

    genetree = ete3.Tree(inputnwk or stdin.read(), format=1)

    #get_taxon = get_taxon_treebest if treebest else get_taxon

    with (open(outputnwk, 'w') if outputnwk else stdout) as out:
        for paralogy in buildparalogies(genetree, get_taxon, ancgene2sp,
                                        ensembl_version,
                                        include_singleton_branches):
            outtext = paralogy.write(format=1, format_root_node=True,
                                     features=['S', 'D', 'P'])
            out.write(outtext + '\n')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('inputnwk', nargs='?', help='Input file or `stdin` if None')
    parser.add_argument('outputnwk', nargs='?', help='Output file or `stdout` if None')
    parser.add_argument('-e', '--ensembl-version', type=int, default=ENSEMBL_VERSION)
    parser.add_argument('-p', '--phyltreefile', default=PHYLTREEFILE)
    parser.add_argument('-t', '--treebest', action='store_true')
    parser.add_argument('-s', '--include-singleton-branches',
                        action='store_true',
                        help='Include single gene nodes, but with a special mark.')
    parser.add_argument('-r', '--include-singleton-root',
                        action='store_true',
                        help='Keep the original root of the gene tree, even ' \
                             'if it\'s not a paralogy. Discard other singletons')

    args = parser.parse_args()
    main(**vars(args))
