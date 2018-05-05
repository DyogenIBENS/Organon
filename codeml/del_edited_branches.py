#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Read the tree forest, delete the longest branches stemming from edited nodes,
and output the new forest."""


import argparse
import LibsDyogen.myProteinTree as ProteinTree

#from glou_duphist import dfw_descendants_generalized


#def prottree_getchildren(tree, node):
#    return [child for child, _ in tree.data.get(node, [])]

INFINITE_DIST = 10000


def lock_targets(nodedata, nodeinfo, tree):
    """Return which children to delete"""
    try:
        targets = {child: i for i,(child,dist) in enumerate(nodedata) if dist >= INFINITE_DIST}
    except ValueError as err:
        err.args += (nodedata,)
        raise

    # Gene splits
    targets.update({child: i for i,(child,_) in enumerate(nodedata) \
            if tree.info[child]['Duplication'] == 10)

    if nodeinfo['Duplication'] == 3:
        children, dists = zip(*nodedata)
        maxdist = max(dists)
        edited = dists.index(maxdist)
        children = list(children)
        edited_child = children.pop(edited)
        targets[edited_child] = edited
    #elif nodeinfo['Duplication'] == 10:
    #    targets[]

    return targets


def knock_targets(targets, tree, nodedata, nodeinfo):
    """Delete target nodes (i.e subtrees) from tree"""
    leaf_count = 0
    for target, target_index in sorted(tuple(targets.items()),
                                        key=lambda item: -item[1]):
        tree.info.pop(target)
        try:
            tree.data.pop(target)
        except KeyError as err:
            leaf_count += 1

        nodedata.pop(target_index)
    return leaf_count


def filterbranches(tree, node):
    del_count = 0
    del_leaf_count = 0

    #if node is None:
    #    node = tree.root

    data = tree.data.get(node)

    if data:
        info = tree.info[node]
        targets = lock_targets(data, info, tree)
        if targets:
            del_count      += len(targets)
            del_leaf_count += knock_targets(targets, tree, data, info)
        
        for child,_ in data:
            ch_dc, ch_dlc = filterbranches(tree, node=child)
            del_count      += ch_dc
            del_leaf_count += ch_dlc

    return del_count, del_leaf_count



def main(treeforestfile, outfile, dryrun=False):
    total_deleted = 0
    total_leaves_deleted = 0
    if dryrun:
        outfile = '/dev/null'
    with open(outfile, 'w') as out:
        for tree in ProteinTree.loadTree(treeforestfile):
            del_count, del_leaf_count = filterbranches(tree, node=tree.root)
            total_deleted        += del_count
            total_leaves_deleted += del_leaf_count
            tree.printTree(out)
            #break
    print("Deleted %d branches, of which %d leaves." % (total_deleted, total_leaves_deleted))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('treeforestfile')
    parser.add_argument('outfile')
    parser.add_argument('-n', '--dryrun', action='store_true',
                        help='Do not delete or output anything, just print '\
                             'counts.')
    
    args = parser.parse_args()
    main(**vars(args))
