#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from sys import stderr
import os
import os.path as op
from glob import glob
import re
import argparse
import ete3


def iter_splitgenes_ancgenes(filename):
    """Example file: ~/ws2/DUPLI_data85/split_genes_info-Rodentia-ancgenes.tsv"""
    with open(filename) as inf:
        for line in inf:
            ancgene, desc_genes, desc_prot = line.rstrip().split('\t')
            yield ancgene, desc_genes.split()


def load_genetreelist(genetreelistfile):
    """Do not use"""
    gt_list = []
    with open(genetreelistfile) as gtlf:
        for line in gtlf:
            if not line.startswith('#'):
                assert line.startswith('ENSGT00')
                gt_list.append(line.rstrip())
    return gt_list


def rewrite_all_genetrees(gt_list, split_genes):
    """Do not use"""
    tot_SG = 0

    for gt in gt_list:
        gt_nwk = op.join(alignments_dir, gt, gt + 'Clean.nwk')
        gt_newnwk = op.join(alignments_dir, gt, gt + 'CleanNoSG.nwk')

        tree = ete3.Tree(gt_nwk, format=1)
        for node in tree.traverse('preorder'):
            if node.name in split_genes:
                assert not node.is_root()
                assert set(node.get_leaf_names()) == split_genes[node.name]
                parent = node.up
                node.detach()
                tot_SG += 1
                while not parent.children:
                    node = parent
                    parent = node.up
                    node.detach()

        tree.write(format=1, format_root_node=True, outfile=gt_newnwk,
                   features=["reinserted"])
    print('Gene splits: %d' % tot_SG, file=stderr)

ANCESTORLIST = ["Afrotheria",
                "Carnivora",
                "Cetartiodactyla",
                "Chiroptera",
                "Insectivora",
                "Lagomorpha",
                "Marsupialia",
                "Neognathae",
                "Rodentia",
                "Simiiformes",
                "Xenarthra"]

def load_cladeof(templatefile="~/ws2/lewitusfamilies/clade_%s.txt"):
    cladeof = {}
    for clade in ANCESTORLIST:
        with open(op.expanduser(templatefile % clade.lower())) as cladef:
            for line in cladef:
                cladeof[line.rstrip()] = clade
    return cladeof


def find_src_files(split_ancgene, cladeof, alignments_dir, src_subtreedir,
                   end='*.nwk', exclude_end='_codeml.nwk'):
    ancestor_i = split_ancgene.index('ENSGT00')
    try:
        suffix_i = split_ancgene.index('.', ancestor_i)
    except ValueError:
        suffix_i = len(split_ancgene)
    clade = cladeof[split_ancgene[:ancestor_i]]
    genetree = split_ancgene[ancestor_i:suffix_i]

    candidate_pattern = op.join(alignments_dir,
                                genetree,
                                src_subtreedir,
                                clade + genetree + end)
    candidate_files = glob(candidate_pattern)
    src_files = [f for f in candidate_files \
        if not f.endswith(exclude_end) and \
           split_ancgene[ancestor_i:].startswith(\
                    op.splitext(op.basename(f))[0].replace(clade, '').split('_')[0])]
    assert len(src_files) == 1, "At genesplit %s\n" \
                    "glob failed with pattern %r: found %d files from %d candidates" \
                                % (split_ancgene, candidate_pattern,
                                    len(src_files), len(candidate_files)) + str(candidate_files)
    src_file, = src_files
    return src_file

def max_delete(node):
    while len(node.up.children) == 1:
        node = node.up
    node.detach()

def edit_split(subtree, split_ancgene, edit_func=max_delete):
    """Delete all descendants from a gene split."""
    splitgene_nodes = subtree.search_nodes(name=split_ancgene)
    assert len(splitgene_nodes) <= 1
    if len(splitgene_nodes) == 0:
        raise LookupError('%r not found.' % split_ancgene)
    splitgene_node, = splitgene_nodes
    return edit_func(splitgene_node)


def keep_1_leaf(node, which='dist', delete=False):
    """Keep only one descendant leaf.
    
    param: `which`: - "topo": keep the closest, by number of nodes;
                    - "dist": keep the closest, by distance;
                    - "alpha" : keep the first, alphabetically"""
    if which == 'alpha':
        kept_leaf = sorted(node.get_leaf_names())[0]
    elif which == 'topo':
        kept_leaf, _ = node.get_closest_leaf(topology_only=True)
    elif which == 'dist':
        kept_leaf, _ = node.get_closest_leaf(topology_only=False)
    node.prune([kept_leaf], preserve_branch_length=True)
    if delete:
        node.delete(prevent_nondicotomic=False, preserve_branch_length=True)
    #return kept_leaf.name


def main(SGlistfile, alignments_dir='.', src_subtreedir='subtreesCleanO2', 
         out_subtreedir='subtreesCleanO2noSG', action='del',
         source_field="ancgene"):
    """First search the corresponding source file and edit the tree."""

    edit_functions = {'del': max_delete,
                      'keep1': keep_1_leaf}

    cladeof = load_cladeof()

    outputted = []
    count_genesplits = 0

    if source_field == "ancgene":
        #search after the ancgene name.
        def iter_src_files():
            for split_ancgene, split_descendants in iter_splitgenes_ancgenes(SGlistfile):
                src_file = find_src_files(split_ancgene, cladeof, alignments_dir, src_subtreedir)
                yield src_file
    elif source_field == "id":
        raise NotImplementedError

    for src_file in iter_src_files():
        out_file = src_file.replace(src_subtreedir, out_subtreedir)
        out_dir = op.dirname(out_file)
        
        if out_file in outputted:
            src_file = out_file
        
        subtree = ete3.Tree(src_file, format=1)
        try:
            kept = edit_split(subtree, split_ancgene, edit_func=edit_functions[action])
        except LookupError:
            if out_file in outputted:
                print('%r not found in already edited tree %s' \
                        % (split_ancgene, src_file), file=stderr)
                continue
            else:
                raise LookupError('%r not found in %s' % (split_ancgene, src_file))

        count_genesplits += 1
        
        if not op.exists(out_dir):
            os.mkdir(out_dir)
        subtree.write(format=1, format_root_node=True, outfile=out_file,
                      features=['reinserted'])
        outputted.append(out_file)

    print('Found %d genesplits.' % count_genesplits, file=stderr)
    print('\n'.join(outputted))



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('SGlistfile')
    #parser.add_argument('genetreelistfile')
    parser.add_argument('alignments_dir', nargs='?', default='.')
    parser.add_argument('--out_subtreedir', default='subtreesCleanO2noSG')
    parser.add_argument('-a', '--action', default='del', choices=['del', 'keep1'])

    args = parser.parse_args()
    main(**vars(args))

