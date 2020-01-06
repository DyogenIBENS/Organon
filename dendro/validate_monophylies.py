#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import sys
import ete3
import argparse as ap


def validate_monophylies(tree: ete3.TreeNode, clades: dict,
                         force_check=False):
    exit = 0
    for clade, leaves in clades.items():
        ismono, cladetype, badleaves = tree.check_monophyly(leaves, 'name')
        if ismono and not force_check:
            print('Clade %s: OK (%d leaves).' % (clade, len(leaves)))
        else:
            exit = 1
            print('Clade %s: NO (%s):' % (clade, cladetype))
            
            for node in tree.traverse('preorder'):
                if node.name.startswith(clade):
                    found_leaves = node.get_leaf_names()
                    print('  * Found node %s' % node.name)
                    print('    with extra leaves: ',
                          ' '.join(set(found_leaves).difference(leaves)))
                    print('    and missing leaves: ',
                          ' '.join(set(leaves).difference(found_leaves)))
                    break
            else:
                print("   * Not found:  '^%s.*'" % clade)

            mrca = tree.get_common_ancestor(leaves)
            print('  * MRCA is', mrca.name)
    return exit


def main(treefile: str, cladelistfile: str, force_check=False):
    tree = ete3.Tree(treefile, format=1)
    clades = {}
    with open(cladelistfile) as f:
        for line in f:
            clade, *leaves = line.split()
            assert clade not in clades, "Duplicate clade line %r" % clade
            clades[clade] = leaves
    return validate_monophylies(tree, clades, force_check)


if __name__ == '__main__':
    parser = ap.ArgumentParser(description=__doc__)
    parser.add_argument('treefile')
    parser.add_argument('cladelistfile')
    parser.add_argument('-f', '--force-check', action='store_true',
                        help='Display diff even when the Ete3 test passed.')
    args = parser.parse_args()
    sys.exit(main(**vars(args)))
