#!/usr/bin/env python3
# -*- coding: utf-8 -*-


"""Return only the ingroup tree (uses the NHX tag: is_outgroup)"""


from sys import stdin, stdout, stderr
import argparse as ap
import ete3
from dendro.parsers import read_multinewick
from codeml.subtrees_stats import find_ingroup_marked


def main():
    parser = ap.ArgumentParser(description=__doc__)
    parser.add_argument('intree', nargs='?', default=stdin,
                        type=ap.FileType('r'))
    parser.add_argument('outtree', nargs='?', default=stdout,
                        type=ap.FileType('w'))
    parser.add_argument('--outformat', '--of', type=int, default=0, 
                        help='[%(default)]')
    args = parser.parse_args()

    for newick in read_multinewick(args.intree):
        tree = ete3.Tree(newick, format=1)
        tree_features = set.union(*(n.features for n in tree.traverse()))
        tree_features -= set(('name', 'dist', 'support'))
        ingrouptree, outgroups = find_ingroup_marked(tree)
        outgroup_leaves = set().union(*(outg.get_leaf_names() for outg in outgroups))
        print('Outgroup leaves:' + ','.join(outgroup_leaves), file=stderr)
        args.outtree.write(ingrouptree.write(features=tree_features,
                                             format=args.outformat,
                                             format_root_node=True) + '\n')


if __name__=='__main__':
    main()
