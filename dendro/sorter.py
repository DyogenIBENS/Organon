#!/usr/bin/env python3
# -*- coding: utf-8 -*-


"""Module containing some algorithms to reorder branches of rooted tree topologies"""


from dendro.bates import dfw_descendants_generalized, \
                         rev_dfw_descendants, \
                         iter_leaves, \
                         bfw_descendants_generalized
import logging
logger = logging.getLogger(__name__)


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


def ladderize(tree, root, get_children=None, heavy_on_top=False, 
              assign=None):
    """Pivot nodes to reorder the branches in a visually nice manner (nodes
    with the most numerous leaves at the bottom).
    
    An example get_children function for the LibsDyogen.myPhylTree class:
    get_children = lambda ph,node: ph.items.get(node, [])
    
    An example assign function:
    assign = lambda ph, node, children: ph.items.update({node: children})
    """
    #TODO: with odd numbers of children: the lightest in the middle.
    
    if not get_children:
        get_children = lambda tree, node: tree.get(node, [])
    if not assign:
        assign = lambda tree, node, children: tree.update({node: children})

    dfw = dfw_descendants_generalized(tree, get_children, include_leaves=False,
                                      queue=[root])
    cumul_mass = {}
    key_ladder = lambda child: cumul_mass.setdefault(child, 1)

    for parent, children in reversed(list(dfw)):
        #logger.debug(parent, children)
        #children.sort(key=key_ladder, reverse=heavy_on_top)
        assign(tree, parent, sorted(children, key=key_ladder, reverse=heavy_on_top))
        cumul_mass[parent] = sum(cumul_mass[ch] for ch in children)
    return cumul_mass


def pyramid_old(tree, root, get_children=None, assign=None):
    """Spread internal nodes outwards as they are more recent.
    Convenient for pretty plotting the tree with older ancestral nodes
    towards the center.
    """
    if not get_children:
        get_children = lambda tree, node: tree.get(node, [])
    if not assign:
        assign = lambda tree, node, children: tree.update({node: children})

    basal_nodes = get_children(tree, root)
    while len(basal_nodes)==1:
        basal_nodes = get_children(tree, basal_nodes[0])

    #sizes = {node: len(list(iter_leaves(tree, get_children, [node]))) for node in basal_nodes}
    
    sorted_base = []
    
    # Put outward by alternating insert and append:
    for node,_ in sorted(sizes.items(), key=lambda item: item[1]):
        sorted_base.reverse()  # So that the heaviest is always down.
        sorted_base.append(node)

    n_base = len(sorted_base)

    for n in sorted_base[:n_base//2]:
        ladderize(tree, n, get_children, heavy_on_top=True, assign=assign)

    for n in sorted_base[n_base//2:]:
        ladderize(tree, n, get_children, assign=assign)

    assign(tree, root, sorted_base)


def pyramid(tree, root, get_children=None, assign=None):
    """Spread internal nodes outwards as they are more recent.
    Convenient for pretty plotting the tree with older ancestral nodes
    towards the center.
    """
    if not get_children:
        get_children = lambda tree, node: tree.get(node, [])
    if not assign:
        assign = lambda tree, node, children: tree.update({node: children})

    sizes = ladderize(tree, root, get_children, assign=assign)

    pivot, split = find_pivot(tree, root, get_children, assign, sizes)
    logger.debug('Pyramid pivot at %r', pivot)
    node = root
    outgroup = []
    while node != pivot:  # Fill the outgroups down to the pivot node.
        children = get_children(tree, node)
        outgroup.extend(children if len(children)==1 else children[:-1])
        node = children[-1]
    outgroup.extend(get_children(tree, node)[:split])

    logger.debug('%d outgroup nodes: %s', len(outgroup), outgroup)
    for n in outgroup:
        ladderize(tree, n, get_children, heavy_on_top=True, assign=assign)
    # No need to re-ladderize the pivot, it was done when computing sizes.


# for each node in levelorder:
#     1. put the smallest child one the same side than the outgroup.
#     2a. break if the outgroup becomes bigger than the biggest child.
#     2b. else, continue only to the biggest child.
def find_pivot(tree, root, get_children=None, assign=None, sizes=None):
    if not get_children:
        get_children = lambda tree, node: tree.get(node, [])
    if not assign:
        assign = lambda tree, node, children: tree.update({node: children})

    if sizes is None:
        def skip_assign(*a):
            pass
        sizes = ladderize(tree, root, get_children, assign=skip_assign)
    nleaves = sizes[root]

    outgroup_size = 0
    queue = [root]
    pivot, psplit = root, 1
    while queue:
        node = queue.pop(0)
        children = get_children(tree, node)
        sorted_children, chsizes = zip(*sorted(((ch, sizes[ch]) for ch in children),
                                               key=lambda t: t[1]))
        #Check to remove:
        assert all(isinstance(s, int) for s in chsizes) and all((ch in children) for ch in sorted_children)

        if any(sorted_ch != ch for sorted_ch, ch in zip(sorted_children, children)):
            # Children order changed: must assign new:
            assign(tree, node, list(sorted_children))
        
        # Find the splitting point (works for polytomy):
        assert outgroup_size < sum(chsizes)
        logger.debug('%r outgroup: %d  VS  other: %d' % (node, outgroup_size, sum(chsizes)))
        for split in range(1, len(children)):
            logger.debug('- split %d: %d VS %d' % (split, outgroup_size + sum(chsizes[:split]), sum(chsizes[split:])))
            if outgroup_size + sum(chsizes[:split]) == sum(chsizes[split:]):
                return node, split
            if outgroup_size + sum(chsizes[:split]) > sum(chsizes[split:]):
                # Choose the split that minimizes the delta:
                previous_delta = sum(chsizes[(split-1):]) - (outgroup_size + sum(chsizes[:(split-1)]))  # >0
                delta = outgroup_size + sum(chsizes[:split]) - sum(chsizes[split:])  # = - (previous_delta - 2*chsizes[0])
                assert previous_delta > 0
                assert delta > 0
                logger.debug('previous delta: %d ; next delta: %d' % (previous_delta, delta))
                if previous_delta <= delta:  # Select split with mininum delta
                    return (pivot,psplit) if (split==1) else (node, split-1)
                else:
                    return node,split

                #return (pivot,psplit) if (split==1) else (node, split)
        
        # There was no good split yet (outgroup is still too small):
        outgroup_size += sum(chsizes[:-1])
        queue.append(sorted_children[-1])
        pivot, psplit = node, split

    # Should not be returned.
    #return pivot, psplit

    # Rewrite:
    #pivot, psplit = root, 1  # Store previous stage
    #children = get_children(tree, root)
    #sorted_children, chsizes = zip(*sorted(((ch, sizes[ch]) for ch in children),
    #                                     key=lambda t: t[1]))
    #outgroup_size = chsizes[0]

    #if any(sorted_ch != ch for sorted_ch, ch in zip(sorted_children, children)):
    #    # Children order changed: must assign new:
    #    assign(tree, node, list(sorted_children))

    #while outgroup_size < nleaves - outgroup_size:
    #    split += 1
    #    if split == len(children):
    #        pivot =  sorted_children[-1]
    #        children = get_children(tree, pivot)
    #        sorted_children, chsizes = zip(*sorted(((ch, sizes[ch]) for ch in children),
    #                                             key=lambda t: t[1]))

    #        if any(sorted_ch != ch for sorted_ch, ch in zip(sorted_children, children)):
    #            # Children order changed: must assign new:
    #            assign(tree, node, list(sorted_children))
    #        split = 1

    #    outgroup_size += chsizes[split]





    #level_deltas = [{}, {}]  # Max possible delta.
    #stored_deltas = {root: nleaves}
    #for node, children in bfw_descendants_generalized(tree, get_children,
    #                                                  queue=[root]):
    #    if node == root or len(children)==1:
    #        continue
    #    alt_deltas = []
    #    for child in children:
    #        chsize = sizes[child]
    #        alt_deltas.append(abs(chsize - (nleaves - chsize)))

    #    stored_delta[node] = level_deltas[0].pop(node)

    #    if not level_deltas[0]:
    #        # We finished one level.
    #        level_deltas.pop(0)      # The current level replaces the parent level.
    #        if min(stored_deltas.values()) < min(level_deltas[0].values()):  # The children didn't improve.
    #            break
    #        level_deltas.append({})  # Init the next level.


