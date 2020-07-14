#!/usr/bin/env python3
# -*- coding: utf-8 -*-


"""
Read clade conversions from file (like produced with `clade_match`).

If the converted name indicates a polyphyly (with a `+`), delete it.
The parent clade therefore becomes a polytomy.
"""


from sys import stdin
import argparse as ap
import ete3
from dendro.bates import iter_distleaves
import logging
logger = logging.getLogger(__name__)
logging.basicConfig()


def load_conversion(conversionfile):
    with open(conversionfile) as stream:
        conversion = {}
        for line in stream:
            if not line.startswith('#'):
                field1, field2, *extra = line.rstrip('\r\n').split('\t')
                if field1 and field2:
                    conversion[field1] = field2
    return conversion


def get_data(tree, nodedist):
    return [(ch, ch.dist) for ch in nodedist[0].children]


def delete_polyphylies(tree, conversion):
    for node in tree.traverse('postorder'):
        newname = conversion.get(node.name, node.name)
        if '+' in newname:
            # This clade is polyphyletic in the reference tree.
            # Transform the parent into a multifurcation.
            
            # Memorize the original distances to leaves, to check if it preserves lengths
            parent = node.up if node.up else node
            distleaves0 = sorted(iter_distleaves(parent, get_data),
                                 key=lambda dl: dl[0].name)

            # Propagate the deleted branch length to its children.
            for child in node.children:
                child.dist += node.dist
            node.delete(prevent_nondicotomic=False, preserve_branch_length=False)
            #TODO: check that newname.split('+') matches the node.children.
            logger.info("Delete node %r:%g (%s)", node.name, node.dist, newname)

            distleaves = sorted(iter_distleaves(parent, get_data),
                                key=lambda dl: dl[0].name)

            assert all((l0 == l) for (l0, dl0), (l, dl) in zip(distleaves0, distleaves)), \
                ('\n'.join('%15s\t%15s' % (l0.name, l.name)
                           for (l0,_), (l,_) in zip(distleaves0, distleaves)))
            assert all( ((dl0-dl) < 1e-5) for (l0, dl0), (l, dl) in zip(distleaves0, distleaves)),\
                ('Leaf distances were not preserved after deletion!!!\n'
                 'At node %r:%g\n' % (node.name, node.dist)
                 + '\n'.join('%15s:%g\t%15s:%g' % (l0.name, dl0, l.name, dl)
                             for (l0,dl0), (l,dl) in zip(distleaves0, distleaves)))


def main(conversionfile, treefile=None):
    conversion = load_conversion(conversionfile)

    if treefile is None:
        treefile = stdin.read()
    tree = ete3.Tree(treefile, format=1)
    
    delete_polyphylies(tree, conversion)

    print(tree.write(format=1, format_root_node=True))


if __name__ == '__main__':
    parser = ap.ArgumentParser(description=__doc__)
    parser.add_argument('conversionfile',
                        help='2 columns file as produced by `clade_match.py`.')
    parser.add_argument('treefile', nargs='?')
    
    args = parser.parse_args()
    main(**vars(args))

