#!/usr/bin/env python3


from sys import stderr
import re
from glob import glob
import os.path as op


def get_genomicus_release(dirname):
    """try to find the number of the genomicus release from the directory"""
    m = re.search(r'/data([0-9]+)/$', dirname)
    if m:
        return int(m.group(1))
    else:
        return 0


def set_genomicus_paths(
    root = "/kingdoms/dyogen/workspace2/workspace2/alouis/GENOMICUS_SVN/",
    release=None):                     
    """ guess genomicus paths (ancGenes, tree, genesST) """

    if release:
        genomicus_dir = op.join(root, "data%d/" % release)
    else:
        print("Auto determine GENOMICUS data location", file=stderr)
        # automatically determine file location (most recent genomicus
        # version)
        genomicus_dirs = glob(op.join(root, "data*/"))
        genomicus_dir = max(genomicus_dirs, key=get_genomicus_release)
        release = get_genomicus_release(genomicus_dir)

    print("Using GENOMICUS release %s" % release, file=stderr)
        
    tree = op.join(genomicus_dir, "PhylTree.Ensembl.%d.conf" % release)
    genesST = op.join(genomicus_dir, "genes", "genesST.%s.list.bz2")
    if op.exists(op.join(genomicus_dir, "SplitGenesEdition")):
        ancGenes_dir = op.join(genomicus_dir, "SplitGenesEdition")
    elif op.exists(op.join(genomicus_dir, "goodThreshold")):
        ancGenes_dir = op.join(genomicus_dir, "goodThreshold")
    elif op.exists(op.join(genomicus_dir, "GoodThreshold")):
        ancGenes_dir = op.join(genomicus_dir, "GoodThreshold")
    else:
        print("ancGenes directory not found", file=stderr)
    ancGenes = op.join(ancGenes_dir, "ancGenes", "all",
                            "ancGenes.%s.list.bz2")
    return tree, genesST, ancGenes
