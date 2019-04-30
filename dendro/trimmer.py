#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Cut tree branches in various ways"""
# Anc. gr. name: tom


import dendro.bates as dclimb
import logging
logger = logging.getLogger(__name__)


def reset_data(tree, node, data):
    node.children = []
    for child, dist in data:
        node.add_child(child=child, dist=dist)


def thin(tree, node, leaves,
         get_data=lambda tree, data: [(ch, ch.dist) for ch in data.children],
         get_name=lambda tree, data: data[0].name
         ):
    """Keep only the intermediate nodes corresponding to the most recent common ancestors
    of the given leaves."""
    raise NotImplementedError
    
    #for leafpath in dclimb.iter_leaf_paths(tree, get_data, [(node,0)]):
    #    if

    newdata = []
    for child, chdist in get_data(tree, node, 0):
        unwanted_child = thin(tree, child, leaves, get_data, get_name)
    set_data(tree, node, data)
    return unwanted


def thin_ete3(root, keptleaves):
    """
    Thin down a tree to retain the minimal topology keeping the most recent
    ancestors of the keptleaves.

    Return the latest common ancestor of the found keptleaves if the root
    supports some keptleaves, else return None.

    This edits the structure (copying might be needed).
    """
    if root.is_leaf():
        return root if root in keptleaves else None
    for child in root.get_children():
        if thin_ete3(child, keptleaves) is None:
            root.remove_child(child)
    if not root.children:
        return None
    if len(root.children) == 1:
        newroot = root.children[0]
        newroot.dist += root.dist  # Transfer the branch length to the *child*
        root.delete(prevent_nondicotomic=False, preserve_branch_length=False)
        root = newroot

    return root


def thin_prottree(tree, root, rootdist, keptleaves):
    """
    Thin down a tree to retain the minimal topology keeping the most recent
    ancestors of the keptleaves.

    Return the latest common ancestor of the found keptleaves if the root
    supports some keptleaves, else return None.

    This edits the structure (copying might be needed).
    """
    if root not in tree.data:
        if root in keptleaves:
            return root, rootdist
        else:
            del tree.info[root]
            return None, None
    
    newnodedata = []
    for child,dist in tree.data[root]:
        thinned_child, thinned_dist = thin_prottree(tree, child, dist, keptleaves)
        if thinned_child is not None:
            newnodedata.append((thinned_child, thinned_dist))
    
    if not newnodedata:
        del tree.data[root]
        del tree.info[root]
        return None, None

    if len(newnodedata) == 1:
        del tree.data[root]
        del tree.info[root]
        root, newdist = newnodedata[0]
        rootdist += newdist  # Transfer the branch length to the *child*
    else:
        tree.data[root] = newnodedata

    return root, rootdist


def fuse_single_child_nodes_ete3(tree, copy=True):
    if copy:
        tree = tree.copy()
    count = 0
    for node in tree.iter_descendants(strategy='postorder'):
        if len(node.children) == 1:
            child, = node.children
            child.dist += node.dist
            # This does not preserve the order of the children,
            # if the parent had >=2 children.
            #node.delete(prevent_nondicotomic=False,
            #            preserve_branch_length=False)
            node_i = node.up.children.index(node)
            node.up.children[node_i] = child
            count += 1

    while len(tree.children) == 1:
        tree = tree.children[0].detach()
        count += 1

    logger.debug('Fused %d nodes', count)
    return tree
