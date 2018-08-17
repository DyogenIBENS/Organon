#!/usr/bin/env python3
# -*- coding: utf-8 -*-


from sys import stderr, setrecursionlimit
import LibsDyogen.myTools as myTools
import LibsDyogen.myFile as myFile
import LibsDyogen.myProteinTree as ProteinTree

from dendron.climb import dfw_descendants_generalized, dfw_pairs_generalized

setrecursionlimit(10000) # YOLO


def main(ensembltree, outputfile):
    get_ch = lambda tree, node: [x[0] for x in tree.data.get(node, [])]
    get_chd = lambda tree, nodedist: tree.data.get(nodedist[0], [])

    count_trees = 0
    count_treenodes = []
    count_splits = 0
    count_split_desc = 0
    with myFile.openFile(outputfile, 'w') as out:
        for tree in ProteinTree.loadTree(ensembltree):
            count_treenodes.append(0)
            for (node, dist), childrendists in dfw_descendants_generalized(tree,
                                    get_chd, queue=[(None, (tree.root, 0))]):
                count_treenodes[-1] += 1
                assert tree.info[node]['duplication'] != 10, \
                        "Unexpected. parent node is a split gene: %s: %s" % \
                                (node, tree.info[node])
                for child, chdist in childrendists:
                    if tree.info[child]['Duplication'] == 10:
                        # It's a gene split
                        # Recurse through all the descendants to remove them
                        count_splits += 1
                        for _, GS_descendant in reversed(list(
                                    dfw_pairs_generalized(tree, get_ch,
                                                          queue=[(None, child)],
                                                          include_root=True))):
                            tree.info.pop(GS_descendant)
                            tree.data.pop(GS_descendant)
                            count_split_desc += 1

                        tree.data[node].remove((child, chdist))

            tree.printTree(out)
    print("%d trees" % count_trees, file=stderr)
    print("treenodes:", " ".join(str(nn) for nn in count_treenodes, file=stderr)
    print("Splits: %d  Split descendants: %d" % (count_splits, count_split_desc), file=stderr)


if __name__ == '__main__':
    args = myTools.checkArgs(
                            [("ensembltree",myTools.File), ("outputfile",str)],
                            [],
                            __doc__)
    main(**args)

