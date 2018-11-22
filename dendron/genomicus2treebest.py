#!/usr/bin/env python3


"""Convert LibsDyogen treeforest to **TreeBest** format (newick with special tags).
"""

from __future__ import print_function


from sys import stdin, stdout, setrecursionlimit
import argparse as ap
import LibsDyogen.myProteinTree as ProteinTree


if __name__ == '__main__':
    parser = ap.ArgumentParser(description=__doc__)
    parser.add_argument('forestfile', nargs='?', default=stdin)
    parser.add_argument('outfile', nargs='?', default=None,
                        help='If none: stdout.\nIf contains "{genetree}", each'\
                             ' individual genetree will be saved in the ' \
                             'corresponding file.')
    args = parser.parse_args()

    setrecursionlimit(20000)

    def close_outfile(outfile):
        pass

    indiv_files = False
    if args.outfile is None:
        outfile = stdout
    elif not '{genetree}' in args.outfile:
        outfile = open(args.outfile, 'w')
    else:
        indiv_files = True
        def get_outfile(tree):
            rootinfo = tree.info[tree.root]
            genetree = rootinfo.get('tree_name', rootinfo['family_name'])
            return open(args.outfile.format(genetree=genetree), 'w')

        for tree in ProteinTree.loadTree(args.forestfile):
            with get_outfile(tree) as outfile:
                tree.printNewick(outfile, withDist=True, withTags=True,
                                 withAncSpeciesNames=True, withAncGenesNames=True,
                                 withID=True)

    if not indiv_files:
        for tree in ProteinTree.loadTree(args.forestfile):
            tree.printNewick(outfile, withDist=True, withTags=True,
                             withAncSpeciesNames=True, withAncGenesNames=True,
                             withID=True)

        if args.outfile is not None:
            outfile.close()

