#!/usr/bin/env python3
# -*- coding: utf-8 -*-


"""Module containing some algorithms to reorder branches of rooted tree topologies"""


from .climber import dfw_descendants_generalized


def rev_dfw_descendants(*args, **kwargs):
    return reversed(list(dfw_descendants_generalized(*args, **kwargs)))


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

### TODO
def ladderize():
    raise NotImplementedError
