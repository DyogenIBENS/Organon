#!/usr/bin/env python3


"""Edit chosen nodes from the tree: keep only their closest leaf."""

from sys import stdin, stdout, setrecursionlimit
import argparse as ap
from collections import defaultdict

from LibsDyogen import myProteinTree
from dendron.climber import iter_leaves, \
                            iter_distleaves, \
                            dfw_descendants_generalized, \
                            rev_dfw_descendants

import logging
logger = logging.getLogger(__name__)
#logger.setLevel(logging.INFO)
logging.basicConfig(level=logging.INFO, format="%(levelname)s:%(name)s:%(funcName)s:%(message)s")

MAXDIST = 10000


def get_data(tree, datum):
    return tree.data.get(datum[0], [])

def get_children(tree, child):
    return [x for x, _ in tree.data.get(child, [])]


def keep_closest_leaf(tree, node):
    leafdists = sorted(iter_distleaves(tree, node, get_data), key=lambda datum: datum[1])
    # This does not keep intermediate nodes.
    tree.data[node] = [leafdists[0]]
    # WARNING: the info of deleted nodes is kept, so don't rely on it!!!


def detach_toolongbranch(tree, node, maxdist=MAXDIST):
    n_detached = 0
    n_leaves_detached = []

    newdata = []
    for ch, dist in tree.data.get(node, []):
        if dist < maxdist:
            newdata.append((ch, dist))
        else:
            n_detached += 1
            n_leaves_detached.append(
                    len(list(iter_leaves(tree, get_children, [child]))))
    tree.data[node] = newdata
    return n_detached, n_leaves_detached


def select_any():
    pass


def edit_what_and_how(proteintrees, what_how, yield_unchanged=False):
    ### WORK in progress
    """selector: tell if a node should be edited.
    editor: choose an edition method (ex: update the data of a given node)."""
    to_delete = {}

    for tree in proteintrees:
        for node_dist, children_dists in dfw_descendants_generalized(tree,
                                                            get_data,
                                                            [(tree.root, 0)]):
            pass
            for reason_why, editor in what_how:
                to_delete[reason_why].update(editor(node_dist, children_dists, tree))


def edit_subtree(node_dist, tree, what_how):
    ### WORK in progress
    n_edits = {what: 0 for what in what_how}
    n_included = {what: 0 for what in what_how}

    for child, dist in tree.data.get(node_dist[0], []):
        # First retrieve the counts of the subtree
        child_n_edits, child_n_included = edit_subtree((child, dist), tree, what_how)

        for child_inc_what, child_inc_count in child_n_included.items():
            n_included[child_inc_what] += child_inc_count

        if editor((child, dist), tree):
            for child_edit_what, child_edit_count in child_n_edits.items():
                for what, editor in what_how.items():
                    n_included[child_what] += count

    return n_edits, n_included


def tree_edit_toolong(tree, maxdist=MAXDIST):
    initial_Nleaves = len(list(iter_leaves(tree, get_children)))

    counts_by_child = {}
    #detached_leaves = set()
    leafset = set()
    detached_ids = set()

    new_leaves = set()
    included_new_leaves = set()

    detached_subtrees = []
    #logger.debug("%d/%d\n%d/%d\n%d - %d = %d\n%d - %d = %d\nlen_iter_leaves = %d",
    #             len(tree.info), len(set(tree.info)),
    #             len(tree.data), len(set(tree.data)),
    #             len(tree.info), len(tree.data), len(tree.info) - len(tree.data),
    #             len(set(tree.info)),len(set(tree.data)),len(set(tree.info)-set(tree.data)),
    #             len(list(iter_leaves(tree, get_children))))
    #logger.debug([n for n,info in tree.info.items() if (n not in tree.data and not info.get('gene_name'))])
    for node_dist, children_dists in rev_dfw_descendants(tree, get_data,
                                                         include_leaves=True,
                                                         queue=[(tree.root, 0)]):
        # (detached, included, detached_leaves, leaves) #, detached_leafset_sizes
        node_counts = [0, 0, 0, 0, []]
        if not children_dists:
            # Add 1 to the leaf count.
            node_counts[3] += 1
            leafset.add(node_dist[0])

        else:
            newdata = []
            for child, dist in children_dists:
                child_counts = counts_by_child.pop(child)
                # Add the leaf count.
                node_counts[3] += child_counts[3]
                
                # Conditionally on this child being detached, update the counts.
                if dist >= maxdist or child in new_leaves:
                    if child in new_leaves:
                        new_leaves.remove(child)
                        included_new_leaves.add(child)

                    node_counts[0] += 1
                    node_counts[1] += sum(child_counts[:2])
                    node_counts[2] += child_counts[3]
                    #node_counts[3] += len(list(iter_leaves(tree, get_children, child)))
                    node_counts[4].append(child_counts[3])
                    detached_ids.add(child)

                    detached_data = {}
                    detached_info = {}
                    for (n,d),ndata in dfw_descendants_generalized(tree,
                                                    get_data,
                                                    include_leaves=True,
                                                    queue=[(child,0)]):
                        detached_info[n] = tree.info.pop(n)  #.pop
                        if ndata:
                            detached_data[n] = ndata

                    detached_subtrees.append(myProteinTree.ProteinTree(
                                                detached_data,
                                                detached_info,
                                                child))
                    ### TODO: add a 'tree_name' at the subtree.info[subtree.root]
                else:
                    node_counts[0] += child_counts[0]
                    node_counts[1] += child_counts[1]
                    node_counts[2] += child_counts[2]
                    node_counts[4].extend(child_counts[4])
                    newdata.append((child, dist))

            #if len(tree.data[node_dist[0]]) != len(newdata):
            #if any((x!=y) for x,y in zip(sorted(tree.data[node_dist[0]]),
            #                             sorted(newdata))):
                #import ipdb; ipdb.set_trace()
            if not newdata:
                logger.warning("All children detached at node %d from tree %d.",
                               node_dist[0], tree.root)
                new_leaves.add(node_dist[0])

            tree.data[node_dist[0]] = newdata

        counts_by_child[node_dist[0]] = node_counts
    
    assert len(counts_by_child) == 1
    root_counts = counts_by_child[tree.root]

    #DEBUG
    
    #import ipdb; ipdb.set_trace()
    #logger.debug("%d/%d\n%d/%d\n%d - %d = %d\n%d - %d = %d\nlen_iter_leaves = %d",
    #         len(tree.info), len(set(tree.info)),
    #         len(tree.data), len(set(tree.data)),
    #         len(tree.info), len(tree.data), len(tree.info) - len(tree.data),
    #         len(set(tree.info)),len(set(tree.data)),len(set(tree.info)-set(tree.data)),
    #         len(list(iter_leaves(tree, get_children))))
    logger.debug("Tree %d: %s,\ndeleted ids: %s",
                 tree.root, root_counts, ' '.join(str(i) for i in detached_ids))
    
    root_Nleaves = len(list(iter_leaves(tree, get_children, [tree.root])))
    assert len(leafset) == initial_Nleaves, "leafset %d != initial_Nleaves %d" %(len(leafset), initial_Nleaves)
    assert len(leafset) == root_Nleaves - len(new_leaves) + root_counts[2], \
        "Tree %d: leafset %d != %d root_Nleaves - %d new_leaves + %d detached_leaves" %\
            (tree.root, len(leafset), root_Nleaves, len(new_leaves), root_counts[2])
    
    assert initial_Nleaves == root_counts[3], "%d != %d, %s" %(initial_Nleaves,
                                                               root_counts[3],
                                                               root_counts)
    assert len(root_counts[4]) == root_counts[0]
    return root_counts, detached_subtrees


def edit_toolong(proteintrees, maxdist=MAXDIST):
    n_detached = 0
    n_included = 0
    n_leaves_detached = 0
    n_leaves_detached_distrib = defaultdict(int)

    for tree in proteintrees:
        root_counts, detached_subtrees = tree_edit_toolong(tree, maxdist)

        n_detached += root_counts[0]
        n_included += root_counts[1]
        n_leaves_detached += root_counts[2]
        for detached_leafset_size in root_counts[4]:
            n_leaves_detached_distrib[detached_leafset_size] += 1
        
        yield root_counts[0] > 0, tree
        for subtree in detached_subtrees:
            yield True, subtree

    logger.info("%d detached; +%d included; %d leaves_detached.",
                n_detached, n_included, n_leaves_detached)
        

def edit_toolong_flagged(flagged_proteintrees, maxdist=MAXDIST):
    n_detached = 0
    n_included = 0
    n_leaves_detached = 0
    n_leaves_detached_distrib = defaultdict(int)

    for flag, tree in flagged_proteintrees:
        root_counts, detached_subtrees = tree_edit_toolong(tree, maxdist)

        n_detached += root_counts[0]
        n_included += root_counts[1]
        n_leaves_detached += root_counts[2]
        for detached_leafset_size in root_counts[4]:
            n_leaves_detached_distrib[detached_leafset_size] += 1
        
        yield flag, root_counts[0] > 0, tree
        for subtree in detached_subtrees:
            #TODO: actually accurately flag if subtree has an inner change
            yield flag, True, subtree
    
    logger.info("\n  %9d detached;\n"
                "  +%8d included;\n"
                "  %9d leaves_detached.",
                n_detached, n_included, n_leaves_detached)



def edit_from_selection(proteintrees, badnodes):
    """Keep only one child, the one leading to the closest leaves."""
    
    n_edits = 0
    n_included = 0  # Number of nodes not edited, because already in a larger
                    # edited subtree
    n_morethan2 = 0  # Number of nodes with >2 leaves

    for tree in proteintrees:
        tree_badnodes = badnodes.intersection(tree.data)
        if not tree_badnodes:
            logger.debug("No nodes to remove.")
            yield False, tree
            continue

        if tree.root in tree_badnodes:
            logger.error("Root node %d listed as a node to remove!!!", tree.root)
            raise NotImplementedError("Root node %d marked for removal: "
                                      "don't know what to do." % tree.root)
        # First the bad nodes must be sorted according to their rootwardness.
        # (in order to edit only the most basal in case of inclusion, and
        # ignore the included ones)
        badnode_leaves = [(badnode,
                           set(iter_leaves(tree, get_children, [badnode])))
                          for badnode in tree_badnodes]

        edited_leafsets = []
        # Indices in the above list where the leafset size decreases
        larger_size_i = 0
        # Current size
        size = len(badnode_leaves[0][1])

        for badnode, leaves in sorted(badnode_leaves, key=lambda x: len(x[1]), reverse=True):

            # Check leaf number: then only check intersection with strictly larger sets.
            new_size = len(leaves)
            if new_size < size:
                #logger.debug("Size decrease!")
                larger_size_i = len(edited_leafsets)
                size = new_size

            # Edit only if no ancestral node has already been edited.
            if not any((leaves & edited_l) for edited_l in edited_leafsets[:larger_size_i]):
                edited_leafsets.append(leaves)
                keep_closest_leaf(tree, badnode)
                n_edits += 1
            else:
                n_included += 1
        
        assert size == 2
        n_morethan2 += larger_size_i

        badnodes.difference_update(tree_badnodes)
        logger.debug("Tree: %d: %d edited nodes; %d implicitely edited; %d with >2 leaves",
                     tree.root, n_edits, n_included, n_morethan2)
        yield True, tree

    if badnodes:
        logger.warning("%d nodes not found: %s",
                       len(badnodes), ' '.join(str(n) for n in badnodes))
    logger.info("\n  %9d edited nodes\n"
                "  +%8d implicitely edited\n"
                "  %9d with >2 leaves",
                n_edits, n_included, n_morethan2)


def or_combine_flagged_iterable(flagged_iterable, flagging_iterator, *args, **kwargs):
    # It is still looping twice instead of once (not efficient).
    flags1, elems1 = zip(*flagged_iterable)
    for flag1, (flag2, elem2) in zip(flags1, flagging_iterator(elems1, *args, **kwargs)):
        yield flag1 | flag2, elem2


def main(badnodelistfile, forestfile, badnode_col=0, maxdist=MAXDIST,
         print_unchanged=True, dryrun=False):
    with open(badnodelistfile) as f:
        header_line = next(f)
        badnodes = set(int(line.rstrip().split()[0]) for line in f)

    logger.info('%d nodes to remove.', len(badnodes))

    proteintrees = myProteinTree.loadTree(forestfile)
    #for has_changed, tree in edit_from_selection(badnodes, proteintrees):
    #for has_changed, tree in or_combine_flagged_iterable(
    #                                edit_from_selection(proteintrees, badnodes),
    #                                edit_toolong,
    #                                maxdist=maxdist):
    n_unprinted = 0

    if dryrun:
        def output(tree, flag1, flag2):
            nonlocal n_unprinted
            n_unprinted += 1
            return int(flag1 | flag2)
    elif print_unchanged:
        def output(tree, flag1, flag2):
            tree.printTree(stdout)
            return int(flag1 | flag2)
    else:
        def output(tree, flag1, flag2):
            if flag1 | flag2:
                tree.printTree(stdout)
                return 1
            else:
                nonlocal n_unprinted
                n_unprinted += 1
                return 0

    n_edited_trees = 0
    n_edited_trees_fromsel = 0
    n_edited_trees_toolong = 0

    for change1, change2, tree in edit_toolong_flagged(
                                    edit_from_selection(proteintrees, badnodes),
                                    maxdist=maxdist):
        n_edited_trees_fromsel += int(change1)
        n_edited_trees_toolong += int(change2)
        n_edited_trees += output(tree, change1, change2)

    logger.info('\n%9d edited trees. %d unprinted trees.\n'
                ' -> %9d from node selection,\n'
                ' -> %9d because of too long branches.\n',
                n_edited_trees, n_unprinted,
                n_edited_trees_fromsel, n_edited_trees_toolong)


if __name__ == '__main__':
    setrecursionlimit(10000)
    parser = ap.ArgumentParser(description=__doc__)
    parser.add_argument('badnodelistfile')
    parser.add_argument('forestfile', nargs='?', default=stdin)
    parser.add_argument('-m', '--maxdist', type=float, default=MAXDIST,
                        help='[%(default)s]')
    parser.add_argument('-d', '--debug', action='store_true')
    parser.add_argument('-n', '--dryrun', action='store_true')
    parser.add_argument('-c', '--changed-only', action='store_false',
                        dest='print_unchanged',
                        help='Output only tree that were changed.')

    dictargs = vars(parser.parse_args())

    if dictargs.pop('debug'):
        logger.setLevel(logging.DEBUG)

    main(**dictargs)

