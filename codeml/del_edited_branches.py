#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Read the tree forest, delete the longest branches stemming from edited nodes,
and output the new forest."""


import argparse
import LibsDyogen.myProteinTree as ProteinTree

#from glou_duphist import dfw_descendants_generalized


#def prottree_getchildren(tree, node):
#    return [child for child, _ in tree.data.get(node, [])]
DEL_COUNT = 0
DEL_LEAF_COUNT = 0

def filternodes(tree, node=None):
    global DEL_COUNT, DEL_LEAF_COUNT
    if node is None:
        node = tree.root
        DEL_COUNT = 0

    data = tree.data.get(node)

    if data:
        info = tree.info[node]

        children, dists = zip(*data)

        # If this node was edited:
        if info['Duplication'] == 3:
            # remove the longest branch (the edited one)
            #print('EDITED NODE: %s -> %s (%s)' % (node, children, dists))
            maxdist = max(dists)
            edited = dists.index(maxdist)
            children = list(children)
            edited_child = children.pop(edited)
            data.pop(edited)
            try:
                tree.data.pop(edited_child)
            except KeyError:
                # the child is a leaf.
                DEL_LEAF_COUNT += 1
                #pass
            finally:
                tree.info.pop(edited_child)
                DEL_COUNT += 1

            #print('NEW CHILDREN: %s (%s)\nNEW DATA: %s' % (children, dists, data))
            #children = [(child, dist) for child, dist in data if dist == maxdist]
        
        for child in children:
            filternodes(tree, node=child)


def main(treeforestfile, outfile):
    total_deleted = 0
    total_leaves_deleted = 0
    with open(outfile, 'w') as out:
        for tree in ProteinTree.loadTree(treeforestfile):
            filternodes(tree)
            total_deleted += DEL_COUNT
            tree.printTree(out)
            #break
    print("Deleted %d branches, of which %d leaves." % (total_deleted, total_leaves_deleted))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('treeforestfile')
    parser.add_argument('outfile')
    
    args = parser.parse_args()
    main(**vars(args))
