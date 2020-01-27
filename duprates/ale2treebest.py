#!/usr/bin/env python3
# -*- coding: utf-8 -*-


"""Read output file from ALE (.uml_rec/.ml_rec: samples of reconciled gene trees), and 
reformat them as for TreeBest.

A new NHX tag is needed: T= for transfer.
"""

import argparse as ap
import re
import ete3
import logging
logger = logging.getLogger(__name__)


def parse_ALEoutput(aleoutputfile):
    #TODO: archiparse
    outputs = {}
    with open(aleoutputfile) as aleoutput:
        header = next(aleoutput).rstrip()
        assert header.startswith('#') and 'ALE' in header

        line = next(aleoutput).rstrip()
        while not line:
            line = next(aleoutput).rstrip()

        assert line.startswith('S:\t'), line
        # Found the species tree
        outputs['Stree'] = line[3:]

        line = next(aleoutput).rstrip()
        while not line:
            line = next(aleoutput).rstrip()

        assert line.startswith('Input ale from:'), line
        outputs['ale input'] = line.split('\t', 1)

        line = next(aleoutput).rstrip()
        assert line.startswith('>logl:'), line
        outputs['logl'] = float(line.split(':')[1])

        line = next(aleoutput).rstrip()
        assert line == 'rate of\t Duplications\tTransfers\tLosses', repr(line)

        line = next(aleoutput).rstrip()
        assert line.startswith('ML \t'), repr(line)
        outputs['rate_Duplications'], \
        outputs['rate_Transfers'], \
        outputs['rate_Losses'] = [float(r) for r in line.split('\t')[1:]]

        line = next(aleoutput).rstrip()
        while not line:
            line = next(aleoutput).rstrip()
        try:
            Nrec = int(re.compile(r'(\d+) reconciled G-s:').search(line).group(1))
        except AttributeError as err:
            err.args += ('at line %r' % line,)
            raise
        line = next(aleoutput).rstrip()
        while not line:
            line = next(aleoutput).rstrip()

        outputs['reconciliations'] = [line]
        for i in range(Nrec-1):
            outputs['reconciliations'].append(next(aleoutput).rstrip())

        line = next(aleoutput).rstrip()
        assert line.startswith('# of\t Duplications\tTransfers\tLosses\tSpeciations'), repr(line)
        line = next(aleoutput).rstrip()
        assert line.startswith('Total \t'), repr(line)
        outputs['total_Duplications'], \
        outputs['total_Transfers'], \
        outputs['total_Losses'], \
        outputs['total_Speciations'] = [float(n) for n in line.split('\t')[1:]]

        line = next(aleoutput).rstrip()
        while not line:
            line = next(aleoutput).rstrip()

        table = {}
        assert line.startswith('# of\t Duplications\tTransfers\tLosses\t'), repr(line)
        # .uml_rec: '# of\t Duplications\tTransfers\tLosses\tOriginations\tcopies\tsingletons\textinction_prob'
        # .ml_rec: '# of\t Duplications\tTransfers\tLosses\tcopies'
        table['header'] = line.split('\t')[1:]
        for line in aleoutput:
            fields = line.rstrip().split('\t')
            table[fields[1]] = [float(x) for x in fields[2:]]
        outputs['branch_stats'] = table

        return outputs

#TODO: quickfix, move this out.
from genomicustools import identify
init_match = {'Scere': 'Saccharomyces cerevisiae', 
              'Celeg': 'Caenorhabditis elegans',
              'Dmela': 'Drosophila melanogaster',
              'Csav': 'Ciona savignyi',
              'Hsa': 'Homo sapiens'}
init_match.update({gid[3:6].capitalize(): sp
                   for sp,gid in identify.SP2GENEID[93].items()
                   if sp not in init_match.values()})

def ale_species_numbers_2_names(ale_Stree, labelled_tree, init_match=None):
    """Numbers correspond to time indices and speciations.
    Here we suppose exact correspondence."""
    label_as_clade = {}
    clade_to_label = {}
    #TODO: check no duplicated node names.
    for lnode in labelled_tree.traverse('postorder'):
        name = lnode.name.split('|')[0]
        if lnode.is_leaf():
            label_as_clade[name] = set((name,))
            clade_to_label[(name,)] = name
        else:
            clade = set.union(*(label_as_clade[ch.name.split('|')[0]] for ch in lnode.children))
            assert name not in label_as_clade, repr(name)
            label_as_clade[name] = clade
            clade_to_label[tuple(sorted(clade))] = name
    
    ale_clades = {}
    matched = {}
    for node in ale_Stree.traverse('postorder'):
        if node.is_leaf():
            matched_name = init_match.get(node.name, node.name)
            assert matched_name in label_as_clade, "%r %s" % (matched_name, label_as_clade)
            ale_clades[node.name] = set((matched_name,))
            matched[node.name] = matched_name
        else:
            try:
                ale_clade = set.union(*(ale_clades[ch.name] for ch in node.children))
            except KeyError as err:
                err.args += ('From %s -> %s' % (node.name, [ch.name for ch in node.children]),)
                raise
            matched[node.name] = clade_to_label[tuple(sorted(ale_clade))]
            ale_clades[node.name] = ale_clade
    return matched


REGEX_TRANSFER = re.compile(r'(D@|T@|Tb@|@)(\d+)\|([^.@|]+)')
# D@  : duplication
# T@  : transfer out
# Tb@ : transfer back
# @   : transfer in.
# .T@ : transfer and loss of the remained copy.

# . expects BRANCH
# D@ T@ .T@ expect TIME|BRANCH with the BRANCH *left*
# @ expects TIME|BRANCH with BRANCH *reached*.
# .T@ expect TIME|BRANCH@TIME|BRANCH ... Because the single descendent is a transfer receiver.

REGEX_LEAFINFO = re.compile(r'([^.@]+)([.@].+)?$')
REGEX_EVENT = re.compile(r'|'.join((r'\.\d+',
                                    r'\.[A-Za-z0-9-]+$',  # non-numeric species name only possible at the end (leaf) 
                                    r'\.?[TD]@\d+\|[^.@|]+',
                                    r'\.?Tb@\d+\|[^.@|]+',
                                    r'@\d+\|[^.@|]+')))
REGEX_LEAFNAME = re.compile(r'([^_]+)(?:_(\d+))?_(.+)')


def alelabel_to_events(label, current, number_to_names):
    #events = label.split('.')
    parent_S = (None if current is None else current.S)
    n_new = 0  # DEBUG
    while label:
        logger.debug("label='%s'", label)
        event_match = REGEX_EVENT.match(label)
        event = event_match.group()
        if event.startswith('.'):
            event = event[1:]
        label = label[event_match.end():]
    #for event in events:
        #if event == '':
        #    continue
        try:
            current = current.add_child()
        except AttributeError:
            assert current is None, "parent node should be ete3.TreeNode or None."
            current = ete3.TreeNode()
        n_new += 1
        m = REGEX_TRANSFER.match(event)
        if not m:
            # Speciation and loss. i.e. transmission to one single descendant species.
            current.add_features(D=0)  # Speciation (and loss)
            branch = event
        else:
            typ, time, branch = m.groups()
            if typ == 'D@':
                current.add_features(D=2)
            elif typ == 'T@':
                # Transfer out
                # should be represented as dupli
                current.add_features(D=11, T=1)
            elif typ == '@':
                # Transfer into that branch:
                # we will represent it simply as a speciation node.
                current.add_features(D=0, T=-1)
            elif typ == 'Tb@':
                # Is it a gene divergence in non represented lineages?
                # Could be represented as Duplication or Speciation, does not matter.
                current.add_features(D=12, T=0)  # Should S==branch?
        current.add_feature('S', number_to_names.get(branch, branch))
    logger.debug('From %r: created %d new nodes leading to species %s',
                 parent_S, n_new, branch)
    return current, n_new


def ale2treebest(tree, numbers_to_names=None):
    if numbers_to_names is None:
        numbers_to_names = {}
    label = tree.name
    treeb, n_new = alelabel_to_events(label, None, numbers_to_names)
    
    alenode2treebestnode = {tree: treeb}  # Map the ale nodes to the corresponding
                               # most recent node in treebest format.
    for parent in tree.traverse('levelorder'):
        for node in parent.children:
            if node.is_leaf():
                leafname, label = REGEX_LEAFINFO.match(node.name).groups()
                try:
                    leafname_split = leafname.split('_')
                    #sp, genename = leafname.rsplit('_', 1)
                    #if len(leafname_split) == 2:
                    #    sp, genename = leafname_split
                    #elif len(leafname_split) == 3:
                    #    sp, count, genename = leafname_split
                    sp, count, genename = REGEX_LEAFNAME.match(leafname).groups()
                except ValueError as err:
                    err.args += (leafname, node.name)
                    raise
                logger.debug('* Leaf %r -> (%r, %r)', node.name, leafname, label)
                if label is None:
                    label = ''  # We feed alelabel_to_events with the correct info.
                label += '.' + sp
            else:
                logger.debug('* Node %r', node.name)
                label = node.name
            # A new current node is branched on top of the parent.
            current, n_new = alelabel_to_events(label, alenode2treebestnode[parent],
                                         numbers_to_names)
            current.name = genename if node.is_leaf() else current.S  # Unnecessary but probably convenient
            #TODO: spread node distances according to the time indices from ALE.
            node_dist_unset = current
            for i in range(n_new):
                node_dist_unset.dist = node.dist/n_new
                node_dist_unset = node_dist_unset.up
            #try:
            #    assert sp == current.S, 'Leaf at %s / %s' %(sp, current.S)
            #except NameError:
            #    pass
            alenode2treebestnode[node] = current
    return treeb


def all_rec_to_treebest(ale_rec_file, treebest_outfile, phyltreefile=None):
    aleout = parse_ALEoutput(ale_rec_file)
    Stree = ete3.Tree(aleout['Stree'], format=1)
    if phyltreefile:
        phyltree = ete3.Tree(phyltreefile, format=1)
        numbers_to_names = ale_species_numbers_2_names(Stree, phyltree, init_match)
    else:
        numbers_to_names = None

    with open(treebest_outfile, 'w') as out:
        for tree in aleout['reconciliations']:
            treeb = ale2treebest(tree, numbers_to_names)
            out.write(treeb.write(format=1, format_root_node=True,
                                  features=['S', 'D', 'T']) + '\n')
            # Sorry, format_root_node=True required for genetree_drawer.


def main():
    parser = ap.ArgumentParser(description=__doc__)
    parser.add_argument('ale_rec_file')
    parser.add_argument('treebest_outfile')
    parser.add_argument('phyltreefile', nargs='?')
    args = parser.parse_args()
    all_rec_to_treebest(**vars(args))

    
if __name__ == '__main__':
    logging.basicConfig()
    main()
