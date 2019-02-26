#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import dendron.climber as dclimb

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


