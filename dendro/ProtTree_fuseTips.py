#!/usr/bin/env python3


"""Fuse the given groups of subspecies by keeping **only one** representative sequence.

Example: 'Heterocephalus glaber female' and 'Heterocephalus glaber male',
    or 'Cricetulus griseus CHOK1GS' and 'Cricetulus griseus Crigri'

Rules:
    - if there exist one gene sequence for each group member (directly steming from the previous node), select a default one.
    - otherwise keep the single sequence.

# TODO: Between males and females, if there exist *different* sequences, with one being on the heterogametic sex chromosome, it should be counted as a duplication.
"""


from sys import stdout, setrecursionlimit
import os.path as op
import argparse as ap
from LibsDyogen import myProteinTree
from dendro.bates import dfw_descendants_generalized, iter_distleaves, iter_leaves
from genomicustools.identify import convert_gene2species

import logging
logger = logging.getLogger(__name__)
logging.basicConfig()
logger.setLevel(logging.INFO)

setrecursionlimit(20000)

# Example species2sub:
species2sub = {'Heterocephalus glaber': ['Heterocephalus glaber female',
                                         'Heterocephalus glaber male'],
               'Cricetulus griseus': ['Cricetulus griseus Crigri',
                                      'Cricetulus griseus CHOK1GS']}

# Ordered by priority: the fusion will keep the first encountered sequence.
species2seq = {'Heterocephalus glaber': ['ENSHGLG00000',
                                         'ENSHGLG00100'],
               'Cricetulus griseus': ['ENSCGRG00000',
                                      'ENSCGRG00001']}

def get_data(tree, nodedist):
    return tree.data.get(nodedist[0], [])

def get_children(tree, node):
    return [c for c,_ in tree.data.get(node, [])]


def fuse_subspecies(forest, species2seq, delete_distant_orthologs=False):
    n_fused = 0
    n_2ch = 0
    n_1ch = 0
    n_single = 0  # Sequence from the given redundant set, without ortholog in the tree.
    n_separated = 0  # Do not share a MRCA in the given species with its apparent orthologs.

    for tree in forest:

        info = tree.info
        data = tree.data
        # Backup 'tree_name' in case the root is deleted
        #tree_name = tree.info[tree.root]['tree_name']
        
        kept_children = set()  # Check variable.
        removed_children = set()
        edited_parents = set()
        for (parent,dist), childrendists in dfw_descendants_generalized(tree, get_data,
                                                queue=[(tree.root, 0)]):
            if info[parent]['Duplication'] != 0:
                continue

            parent_taxon = info[parent]['taxon_name']
            if parent_taxon in species2seq:
                
                assert 'gene_name' not in info[parent]
                data.pop(parent)
                assert len(childrendists) <= len(species2seq[parent_taxon]), \
                        "Too many descendant sequences at node %d (%s)" \
                        % (parent, parent_taxon)
                if len(childrendists) > 1:
                    # Test each sequence start, and quit at first match.
                    for seqstart in species2seq[parent_taxon]:
                        for ch_i, (ch, chdist) in enumerate(childrendists):
                            gene_names, gene_dists = zip(*(
                                (info[tip]['gene_name'], gdist)
                                for tip, gdist in iter_distleaves(tree, get_data,
                                                                  root=ch)
                                          ))

                            if len(gene_names) > 1:
                                assert info[ch]['Duplication'] != 0
                            
                            if gene_names[0].startswith(seqstart):
                                # No inner speciation should be possible.
                                assert all(gn.startswith(seqstart) for gn in gene_names)
                                break
                        else:
                            raise ValueError("%s not matched by any of %s" % \
                                             (gene_names, species2seq[parent_taxon]))
                    
                    #except KeyError as err:
                    #    err.args += ("at %s" % [ch for ch,_ in childrendists],)
                    #    raise
                    
                    kept_ch, kept_dist = ch, chdist
                    removed_ch = [c for c,_ in (childrendists[:ch_i]
                                                + childrendists[(i+1):])]
                    for rc in removed_ch:
                        info.pop(rc)
                        rdat = data.pop(rc, None)
                        if rdat:
                            logger.warning('Removed child %d had descendants: %s',
                                           rc, rdat)
                    removed_children.update(removed_ch)
                    n_2ch += 1

                else:
                    kept_ch, kept_dist = childrendists[0]
                    n_1ch += 1

                if kept_dist > 0:
                    logger.warning("%d %r sequence distance > 0!",
                                   kept_ch, info[kept_ch].get('gene_name'))
                #if 850 in removed_children:
                #    import ipdb; ipdb.set_trace()

                # Replace with the correct child.
                info[parent].update(info.pop(kept_ch))  # So that 'tree_name' or 'Bootstrap' fields are conserved.
                
                try:
                    data[parent] = data.pop(kept_ch)
                except KeyError as err:
                    #if err.args[0] != kept_ch:
                    #    err.args += ('parent = %s ; kept_ch = %s' % (parent, kept_ch),)
                    #    raise
                    pass
                #info[parent].update(taxon_name=parent_taxon)
                
                kept_children.add(kept_ch)
                edited_parents.add(parent)
                n_fused += 1
            else:
                for i, (ch,_) in enumerate(childrendists):
                    for sp, spseqs in species2seq.items():
                        for j, seq in enumerate(spseqs):
                            if info[ch].get('gene_name', '').startswith(seq):
                                logger.warning("Unexpected parent: %d (%s) at leaf %d %s",
                                               parent, parent_taxon,
                                               ch, info[ch]['gene_name'])
                                break
                        else:
                            continue
                                               
                        # If match, Need to check that there is no sister
                        # sequence in the other child
                        sister_childrendists = childrendists[:i] + childrendists[(i+1):]
                        sister_seqs = spseqs[:j] + spseqs[(j+1):]
                        has_orthologs = 0
                        for sister_ch,_ in sister_childrendists:
                            sister_matched = [info[tip]['gene_name']
                                              for tip in iter_leaves(tree,
                                                    get_children, [sister_ch])
                                              for sseq in sister_seqs
                                              if info[tip]['gene_name'].startswith(sseq)]
                            if sister_matched:
                                has_orthologs += len(sister_matched)
                                logger.warning('Potential subspecies orthologs found: '
                                               '%s', sister_matched)
                        
                        # But exclude those if we find paralogs
                        paralogs = [info[tip]['gene_name'] for tip in
                                        iter_leaves(tree, get_children, [ch])
                                    if info[tip]['gene_name'].startswith(seq)]
                        
                        if has_orthologs and not paralogs:
                            n_separated += 1
                            # Comment: there was none of those in Ensembl 93.
                        else:
                            n_single += 1
        assert not removed_children.intersection(info)
        assert not removed_children.intersection(data)
        assert not kept_children.intersection(info)
        assert not kept_children.intersection(data)
        assert len(edited_parents) == len(edited_parents.intersection(info))
        
        yield tree

    logger.info("\n%d fusions (%d from >2 sequences, %d from 1 sequence).\n"
                "%d singles (only one of the two subspecies was found)\n"
                "%d separated (1 or more distant orthologs in the sister "
                "subspecies were found)",
                n_fused, n_2ch, n_1ch, n_single, n_separated)

###TODO: dynamic prog from leaves to root, storing the subspecies orthologs, until treated.
def rootward_fusesub(tree, currentnode):
    currentdata = tree.data.get(currentnode)
    nextsubseqs = {}
    if not currentdata:
        gene_name = tree.info['gene_name']
        for sp, spseqs in species2seq:
            for seqstart in spseqs:
                if gene_name.startswith(seqstart):
                    sp_subseqs = nextsubseqs.setdefault(sp, {})
                    sp_subseqs[seqstart] = currentnode
                    break
            else:
                continue
            break
        return nextsubseqs

    # Fuse if compatible pairs
    for ch, dist in currentdata:
        for sp, subseq2node in rootward_fusesub(tree, ch):
            nextsubseqs[sp].update()

    return currentsubseqs


def main(ensembl_version, forestfile, outfile=None, delete_distant_orthologs=False):
    out = stdout if outfile is None else open(outfile, 'w')
    for fused_tree in fuse_subspecies(
            myProteinTree.loadTree(op.expanduser(forestfile % ensembl_version)),
            species2seq,
            delete_distant_orthologs):
        fused_tree.printTree(out)
    if outfile is not None:
        out.close()


if __name__ == '__main__':
    parser = ap.ArgumentParser(description=__doc__)
    parser.add_argument('ensembl_version', type=int)
    parser.add_argument('forestfile', default='~/GENOMICUS%d/tree.1.ensembl.bz2',
                        nargs='?')
    parser.add_argument('-o', '--outfile')
    args = parser.parse_args()
    main(**vars(args))

