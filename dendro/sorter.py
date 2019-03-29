#!/usr/bin/env python3
# -*- coding: utf-8 -*-


"""Module containing some algorithms to reorder branches of rooted tree topologies"""


from dendron.climber import dfw_descendants_generalized, rev_dfw_descendants


def children_sort(tree, get_children=None, attribute='name'):
    """Rotate sister branches according to the node attribute"""


def leaf_sort(tree, root, get_children, assign_children=None, get_attribute=None,
              reverse=False):
    """Rotate sister branches according to the leaf attribute (e.g. the name).

    Sort **inplace**.

    Example of `get_attribute` to get the 'gene_name' of a leaf in a
        Dyogen ProteinTree:

        >>> get_children = lambda tree, node_and_dist: tree.data.get(node_and_dist[0], [])
        >>> get_attribute = lambda tree, node_and_dist: tree.info[node_and_dist[0]].get('gene_name')
    """

    if get_attribute is None:
        def get_attribute(tree, node):
            return node

    # It might not be necessary to reassign the children. Default is "do nothing".
    if assign_children is None:
        def assign_children(tree, node, children):
            pass

    previous_attributes = {}
    children_get_attribute = lambda child: previous_attributes[child]

    for node, children in rev_dfw_descendants(tree, get_children,
                                              include_leaves=True,
                                              queue=[root]):
        if not children:
            # This is a leaf. Store the attribute information
            previous_attributes[node] = get_attribute(tree, node)

        else:
            children.sort(key=children_get_attribute, reverse=reverse)
            assign_children(tree, node, children)
            previous_attributes[node] = min(previous_attributes.pop(ch) for ch in children)


def ladderize_bylevels(tree, root, get_children=None, short_on_top=True):
    """Pivot nodes to reorder the branches in a visually nice manner (branches
    with the most numerous nodes (on *one* lineage) at the bottom)"""
    if not get_children:
        get_children = lambda tree, node: tree.get(node, [])

    dfw = dfw_descendants_generalized(tree, get_children, include_leaves=False,
                                      queue=[root])
    cumul_dist = {}
    key_ladder = lambda child: cumul_dist.setdefault(child, 0)
    for parent, children in reversed(list(dfw)):
        #print parent, children
        children.sort(key=key_ladder, reverse=short_on_top)
        cumul_dist[parent] = max(cumul_dist[ch] for ch in children) + 1
        tree[parent] = children


def ladderize(tree, root, get_children=None, light_on_top=True, 
              assign=None):
    """Pivot nodes to reorder the branches in a visually nice manner (nodes
    with the most numerous leaves at the bottom).
    
    An example get_children function for the LibsDyogen.myPhylTree class:
    get_children = lambda ph,node: ph.items.get(node, [])
    
    An example assign function:
    assign = lambda ph, node, children: ph.items.update({node: children})
    """

    if not get_children:
        get_children = lambda tree, node: tree.get(node, [])
    if not assign:
        assign = lambda tree, node, children: tree.update({node: children})

    dfw = dfw_descendants_generalized(tree, get_children, include_leaves=False,
                                      queue=[root])
    cumul_mass = {}
    key_ladder = lambda child: cumul_mass.setdefault(child, 1)

    for parent, children in reversed(list(dfw)):
        #print parent, children
        #children.sort(key=key_ladder, reverse=light_on_top)
        assign(tree, parent, sorted(children, key=key_ladder, reverse=light_on_top))
        cumul_mass[parent] = sum(cumul_mass[ch] for ch in children)

