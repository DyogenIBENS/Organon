#!/usr/bin/env python3
# -*- coding: utf-8 -*-


from sys import stderr, setrecursionlimit, exit
import os.path as op
import argparse as ap
import pandas as pd
import multiprocessing as mp
from multiprocessing_logging import install_mp_handler
from functools import partial
from copy import deepcopy
from collections import defaultdict

from LibsDyogen import myPhylTree, myProteinTree, myFile
from lateralite.find_candidates import tree_extract_genecounts
from genomicustools.identify import convert_gene2species

import logging
logger = logging.getLogger(__name__)

ncores = 3


def get_robusts_from_prottree(dataset, prottreefile, ancestors, phyltrees):

    ###WARNING: dubious gene trees could be marked as "robust"!!

    dataset_robusts = pd.DataFrame(0, index=ancestors,
                                   columns=('all_crown', 'implicit',
                                            'robusts', 'Ngenes=Nspecies'),
                                   dtype=int)
                                   #name='%s_%s' % (dataset, edition))
    prottrees = list(myProteinTree.loadTree(op.expanduser(prottreefile)))
    for ancestor in ancestors:
        logger.info('%s: %s', dataset, ancestor)
        anc_genecounts, anc_spgenes, anc_branches = \
                tree_extract_genecounts(
                        deepcopy(prottrees),
                        ancestor,
                        phyltrees[dataset],
                        onlybasal=True)

        dataset_robusts.loc[ancestor, 'all_crown'] = anc_genecounts.shape[0]
        dataset_robusts.loc[ancestor, 'implicit'] = \
                (anc_branches != ancestor).all(axis=1).sum()
        dataset_robusts.loc[ancestor, 'robusts'] = \
                (anc_genecounts == 1).all(axis=1).sum()
        dataset_robusts.loc[ancestor, 'Ngenes=Nspecies'] = \
                (anc_genecounts.sum(axis=1) == len(phyltrees[dataset].species[ancestor])).sum()
    
    return dataset, dataset_robusts


def ancgenes_extract_genecounts(ancgenes_file, ensembl_version,
                                spset=set(('Homo sapiens',))):
    anc_genecounts = []
    anc_spgenes = []
    ancgenes = []
    with myFile.openFile(ancgenes_file, 'r') as f:
        for line in f:
            ancgene, *genes = line.split()
            ancgenes.append(ancgene)

            gcounts = defaultdict(int)
            spgenes = defaultdict(list)
            for g in genes:
                species = convert_gene2species(g, ensembl_version)
                gcounts[species] += 1
                if species in spset:
                    spgenes[species].append(g)

            anc_genecounts.append(gcounts)
            anc_spgenes.append(spgenes)

    return ancgenes, anc_genecounts, anc_spgenes


def get_robusts_from_ancgenes(ancgenes_tmpl, ancestors, phyltree, ensembl_version):
    dataset_robusts = pd.DataFrame(0, index=ancestors,
                                   columns=('all_stem', 'robusts', 'Ngenes=Nspecies'),
                                   dtype=int)
    for ancestor in ancestors:
        logger.info('%s: %s', ensembl_version, ancestor)

        ancgenes, anc_genecounts, _ = ancgenes_extract_genecounts(
                                        ancgenes_tmpl % phyltree.fileName[ancestor],
                                        ensembl_version)
        anc_genecounts = pd.DataFrame(anc_genecounts, index=ancgenes)\
                         .rename_axis('%s_ancestral_gene' % ancestor)
        #anc_spgenes = pd.DataFrame(anc_spgenes, index=ancgenes)\
        #                 .rename_axis('%s_ancestral_gene' % ancestor)
        dataset_robusts.loc[ancestor, 'all_stem'] = anc_genecounts.shape[0]
        dataset_robusts.loc[ancestor, 'robusts'] = \
                (anc_genecounts == 1).all(axis=1).sum()
        dataset_robusts.loc[ancestor, 'Ngenes=Nspecies'] = \
                (anc_genecounts.sum(axis=1) == len(phyltree.species[ancestor])).sum()

    return dataset_robusts


# DEFAULT VALUES
outbasename = 'default'
phyltreefiles = {'75_ed': '~/GENOMICUS75/PhylTree.Ensembl.75.conf',
                 '85_noed': '~/ws2/DUPLI_data85/PhylTree.TimeTree2018.Ensembl-like.nwk',
                 '93_noed': '~/ws2/DUPLI_data93/PhylTree.TimeTree201901.Ensembl93-like.nwk'}

prottrees = {'75_ed': '~/GENOMICUS75/GoodThreshold/tree.4F.cut.bz2',
             #'75_edAm': {'75': '~/GENOMICUS75/GoodThreshold/tree.4F.cut.bz2'},
             '85_noed': '~/ws2/DUPLI_data85/tree.1.genomicus.bz2',
             '93_noed': '~/ws2/DUPLI_data93/tree.1.genomicus.bz2'}

params = []
for dataset in ('75_ed', '85_noed', '93_noed'):
    prottreefile = prottrees.get(dataset)
    if prottreefile is not None:
        params.append((dataset, prottreefile))


def load_each_file_once(input_data, loadfunc, *args, **kwargs):
    loaded = {}
    seen_files = {}  # Reverse search, to load each file only once.
    for dataset,fn in input_data.items():
        fn = op.realpath(op.expanduser(fn))
        fn_datasets = seen_files.setdefault(fn, set())
        if fn_datasets:
            loaded[dataset] = loaded[list(fn_datasets)[0]]
        else:
            loaded[dataset] = loadfunc(fn, *args, **kwargs)
        fn_datasets.add(dataset)
    return loaded


if __name__ == '__main__':
    logging.basicConfig()
    logger.setLevel(logging.INFO)

    parser = ap.ArgumentParser(description=__doc__)
    parser.add_argument('paramfile', nargs='?',
                        help=('tabular file with 3 columns: '
                              'the dataset name, the PhylTree path, '
                              'the ProteinTree forest path.\n'
                              'If none, will run with default values.'))
    parser.add_argument('-l', '--lineage-to', default='HomoPan')

    args = parser.parse_args()

    if args.paramfile:
        params = []
        outbasename = op.basename(op.splitext(args.paramfile))
        phyltreefiles = {}
        prottrees = {}
        with open(args.paramfile) as f:
            for line in f:
                if not line.lstrip(' ').startswith('#'):
                    fields = [field.strip(' ') for field in line.rstrip().split('\t')]
                    if len(fields) != 3:
                        print('Bad input file format: expected 3 tab-separated columns.',
                              file=stderr)
                        exit(2)
                    phyltreefiles[fields[0]] = phyltreefiles[fields[1]]
                    params.append(tuple(fields)[:2])
                
    logger.info('Will process:\n%s', '\n'.join(str(p) for p in params))
    
    setrecursionlimit(10000)

    phyltrees = load_each_file_once(phyltreefiles, myPhylTree.PhylogeneticTree)

    phyltree_ref = phyltrees[sorted(phyltrees.keys())[0]]
    other_phyltrees = [phyltrees[dataset] for dataset in sorted(phyltrees.keys())[1:]]
    other_ancestors = frozenset.intersection(*(p.listAncestr for p in other_phyltrees))
    ancestors = [link
                 for link in phyltree_ref.dicLinks[phyltree_ref.root][args.lineage_to]
                 if len(phyltree_ref.items[link]) != 1 and link in other_ancestors]

    logger.info('List of ancestors:\n%s', ' '.join(ancestors))

    # Interestingly, `partial` won't work with multiprocessing, if given args
    # are not picklable, but a closure will.
    #process = partial(get_robusts_from_prottree, ancestors=ancestors, phyltrees=phyltrees)
    def process(dataset, prottreefile):
        return get_robusts_from_prottree(dataset, prottreefile, ancestors, phyltrees)

    install_mp_handler(logger)
    with mp.Pool(ncores) as pool:
        all_datasets, all_counts = zip(*pool.starmap_async(process, params).get())

    print('Concatenating...', file=stderr)
    out = pd.concat(all_counts, axis=1, keys=all_datasets)
    print('Writing to file...', file=stderr)
    out.to_csv('%s_robust_counts.tsv' % outbasename, sep='\t')

