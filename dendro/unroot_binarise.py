#!/usr/bin/env python3
# -*- coding: utf-8 -*-


from sys import stdin, stdout
import argparse as ap
import ete3


def unroot_binarise(intreefile, outtreefile):
    tree = ete3.Tree(intreefile.read())
    tree.resolve_polytomies()
    tree.unroot()
    outtreefile.write(tree.write(format=5, quoted_node_names=False))



def main():
    parser = ap.ArgumentParser(description=__doc__)
    parser.add_argument('intreefile', nargs='?', type=ap.FileType('r'),
                        default=stdin)
    parser.add_argument('outtreefile', nargs='?', type=ap.FileType('w'),
                        default=stdout)
    args = parser.parse_args()
    unroot_binarise(**vars(args))


if __name__ == '__main__':
    main()
