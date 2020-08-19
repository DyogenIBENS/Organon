#!/usr/bin/env python3

"""
Helper for batch runs of `genetree_drawer.py`:

Automatically prepare the data for drawing, from a single Ensembl genetree name.

Warnings: experimental. hard-coded paths and file structure.
"""

import sys
import os
import os.path as op


ENSEMBL_VERSION = 85


def check_extracted(genetrees, output):
    all_outputs = [output.format(genetree=gt) for gt in genetrees]
    new_genetrees = [gt for gt in genetrees if not op.exists(output.format(genetree=gt))]
    new_outputs = [output.format(genetree=gt) for gt in new_genetrees]
    print("The following genetrees are already extracted: %s" %
            (set(all_outputs) - set(new_outputs),), file=sys.stderr)
    return all_outputs, new_outputs, new_genetrees


def prepare(genetrees, ancestors, ensembl_version=ENSEMBL_VERSION, ncores=1,
            edited=True, subtrees_dir='subtrees_', dry_run=False):
    """
    Prepare genetree files (if needed) starting from raw Genomicus/Ensembl data:
      1. find them in the tree forest;
      2. reconcile them with the species tree;
      3. extract the subtrees starting at a given ancestor.
    """
    # 1. Find ancgene name from modern gene?
    # 2. Find the orthologs/ancestors in the given ancestors

    datadir = op.expanduser('~/ws2/DUPLI_data%d/alignments' % ensembl_version)
    assert op.exists(datadir)

    if edited:
        # Take gene tree from Genomicus
        treeforestfile = op.expanduser("~/GENOMICUS%d/"
                         "GoodThreshold/tree.4F.cut.bz2" % ensembl_version)
        withAncGenesNames = True
        field = 'family_name'
        output = op.join(datadir, '{genetree}', '{genetree}.nwk')
        fix_suffix = True
    else:
        # Take gene tree from Ensembl
        treeforestfile = op.expanduser("~/GENOMICUS%d/"
                         "tree.1.ensembl.bz2" % ensembl_version)
        withAncGenesNames = False
        field = 'tree_name'
        output = op.join(datadir, '{genetree}', '{genetree}_ensembl.nwk')
        fix_suffix = False

    all_outputs, new_outputs, new_genetrees = check_extracted(genetrees, output)
    
    # Extract genetrees
    if new_outputs:
        import ToolsDyogen.treeTools.ALL.extractMultipleGeneTree as xMulti
        xMulti.main(treeforestfile, new_genetrees, field=field, toNewick=True,
                    withAncSpeciesNames=True,
                    withAncGenesNames=withAncGenesNames,
                    output=output, mkdirs=True)
    else:
        print("No new genetrees to extract.", file=sys.stderr)

    # Create output dirs for prune2family
    gt_format = '{0}' if edited else '{0:.%d}' % len(genetrees[0])
    #print(gt_format)
    
    prune_outdir = op.join(datadir, gt_format, subtrees_dir)

    for gt in genetrees:
        p_outdir = prune_outdir.format(gt)
        print(p_outdir, end=' ')
        if not op.exists(p_outdir):
            print("make!")
            os.mkdir(p_outdir)
        else:
            print("exists.")

    # prune genetrees to the right ancestor
    import genchron.prune2family as prune2family

    return prune2family.parallel_save_subtrees(all_outputs, ancestors,
                                           ncores=ncores,
                                           outdir=prune_outdir,
                                           only_dup=False, one_leaf=True,
                                           fix_suffix=fix_suffix,
                                           dry_run=dry_run, ignore_errors=False,
                                           ensembl_version=ensembl_version)



