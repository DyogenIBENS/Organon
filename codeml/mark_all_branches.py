#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
"""Output all trees with one foreground branch marked ('#1') for codeml branch models.
If there are internal node names, they are used as suffix of the output file.
"""

import os.path
import ete3
import argparse


def mark_each(nwfile, format=1, outbase=None):
    """iterate each branch and save a new tree with the marked branch"""
    inbase, ext = os.path.splitext(nwfile)
    if not outbase: outbase = inbase

    tree = ete3.Tree(nwfile, format=format)
    node_nb = 0
    for node in tree.iter_descendants():
        oldname = node.name
        if oldname:
            node.name += ' #1'
            outfile = outbase + '-' + oldname + ext
        else:
            node.name = '#1'
            node_nb += 1
            outfile = outbase + '-%02d' % node_nb + ext

        tree.write(outfile=outfile, format=1)
        node.name = oldname


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('nwfile')
    parser.add_argument('-f', '--format', default=1, type=int,
                        help='newick subformat [%(default)s]')
    parser.add_argument('-o', '--outbase')
    
    args = parser.parse_args()
    mark_each(**vars(args))

