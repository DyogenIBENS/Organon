#!/usr/bin/env python3
# -*- coding: utf-8 -*-


"""Working with trees stored in dataframes (using Pandas)"""


from math import isnan
import pandas as pd
import ete3
from dendro.bates import dfw_pairs_generalized
from dendro.any import ete3 as ete3_methods, \
                       myPhylTree as phyltree_methods, \
                       myProteinTree as prottree_methods
import logging
logger = logging.getLogger(__name__)


# ~~> dendro.bates.framed
# FIXME: too slow.
def roll_rootwards_indices(df, parent_column='parent', type_column=None, node_column=None,
                           strategy='levelorder'): #tree_column=None, 
    """Given a dataframe with rows representing tree nodes,
    Iterate row indices in a postorder manner:

    yield (parent_i, tuple(children_is))

    If the parent is the root, yield None. [Not implemented: relies on the parent being NaN]
    """
    if node_column is not None:
        df = df.set_index(node_column)
        if df.index.has_duplicates:
            logger.error('Should not index forest with duplicated node names.')

    sister_groups = df.groupby(parent_column, sort=False)
    sister_items = [(parent, g.index) for parent,g in sister_groups]
    logger.debug('%d node groups.', sister_groups.ngroups)
    # Find leaves:
    if type_column is None:
        # this could be done using:
        leaves = set(df.index.difference(sister_groups.groups))
    else:
        leaves = set(df.index[df[type_column] == 'leaf'])
    # Better with list(leaves) ? (order accelerating the search)

    if strategy == 'levelorder':
        yield from levelorder_rootwards(df, sister_groups.groups, leaves)
    elif strategy == 'postorder':
        yield from postorder_rootwards(df, sister_items, leaves)
    else:
        raise ValueError('Invalid `strategy` ["levelorder"/"postorder"]')


def get_sorted_sister_groups(df, parent_column):
    # This is going to be awfully slow but let's try
    # groupby(sort=False) does the same.
    pass


def levelorder_rootwards(df, sister_groups, leaves):
    # suboptimal
    # Worst case scenario: with dataframe of N rows:
    # one single caterpillar tree, and the cherry is the last encountered sister_group.
    # Will iterate N*(N-1)/2 times
    max_iter = df.shape[0]
    max_iter *= (max_iter - 1)/2.
    i = 0
    # You only want to fetch those parents whose subtrees have been visited.
    # Therefore, visited nodes will be set as leaves as we go rootwards.
    while sister_groups and i <= max_iter:
        logger.debug('rootwards: iteration #%d (%d groups, %d leaves remaining)',
                     i, len(sister_groups), len(leaves))
        for parent, children in list(sister_groups.items()):
            i += 1
            if children.difference(leaves).empty:
                # Both children are current leaves: proceed rootward.
                yield parent, children
                # This may work, but unsure of the side effects.
                del sister_groups[parent]
                leaves.difference_update(children)
                leaves.add(parent)
                #break  if not levelorder, we can break here.
    if i == max_iter:
        logger.error('Reached maximum iter. Tree is probably wrong.')
    else:
        logger.debug('Finished in %d iterations.', i)


def postorder_rootwards(df, sister_items, leaves, is_sorted=True):
    # suboptimal
    # Worst case scenario: with dataframe of N rows:
    # one single caterpillar tree, and the cherry is the last encountered sister_item.
    # Will iterate N*(N-1)/2 times
    max_iter = df.shape[0]
    max_iter *= (max_iter - 1)/2.
    i = 0
    while sister_items and i <= max_iter:
        #logger.debug('rootwards: iteration #%d (%d groups, %d leaves remaining)',
        #             i, len(sister_items), len(leaves))
        parent, children = sister_items.pop(0)
        i += 1
        if children.difference(leaves).empty:  # If sorted, always true.
            # Both children are current leaves: proceed rootward.
            yield parent, children
            # This may work, but unsure of the side effects.
            leaves.difference_update(children)
            leaves.add(parent)
            #break  if not levelorder, we can break here.
        else:
            #sister_items.insert(1, (parent, children))
            # Track the next item having `parent` in the children
            if is_sorted:
                raise RuntimeError('Expected sorted nodes (postorder).')
            pass
            
    if i == max_iter:
        logger.error('Reached maximum iter. Tree is probably wrong.')
    else:
        logger.debug('Finished in %d iterations.', i)


def roll_leafwards_indices(df, parent_column='parent', root_value=None,
                           type_column=None, node_column=None,
                           include_leaves=False): #tree_column=None, 
    """Given a dataframe with rows representing tree nodes,
    Iterate row indices in a postorder manner:

    yield (parent_i, tuple(children_is))

    If the parent is the root, yield None.
    """
    if node_column is not None:
        df = df.set_index(node_column)
        if df.index.has_duplicates:
            logger.error('Should not index forest with duplicated node names.')

    sister_grouped = df.groupby(parent_column)

    # Find roots:
    if root_value is None:
        roots = df.index[df[parent_column].isin((root_value,))].tolist()
    else:
        roots = sister_grouped.get_group(root_value).index.tolist()

    while roots:
        node = roots.pop()  # .pop(0) for a levelorder.
        try:
            children = sister_grouped.get_group(node).index
        except KeyError:
            # node is a leaf. skip.
            if include_leaves:
                children = pd.Index([])
            else:
                continue
        yield node, children
        roots.extend(children)


# ~~> dendro.converters.framed
def to_parentdata(tree, get_items, root_item, get_label=None):
    """Encodes a tree into a DataFrame with index 'node' and columns ['parent', 'dist']
    """
    if get_label is None:
        def get_label(tree, node): return node

    parents, nodes, dists = [], [], []
    for (parent,_), (child, dist) in dfw_pairs_generalized(tree,
                                              get_items,
                                              queue=[((None, 0), root_item)]):
        parents.append(None if parent is None else get_label(tree, parent))
        nodes.append(get_label(tree, child))
        dists.append(dist)

    return pd.DataFrame({'parent': parents, 'dist': dists}, index=nodes)


def ete3_to_parentdata(tree):
    return to_parentdata(tree, ete3_methods.get_items, (tree, tree.dist),
                         ete3_methods.get_label)


def parentdata_to_ete3(df, dist_column='dist', root_value=None): #, parent_column='parent'
    roots = []
    trees = {}
    #if isinstance(root_value, float) and isnan(root_value):
    #    def is_root(node): return isnan(df.loc[node, parent_column])
    #else:
    #    def is_root(node): return df.loc[node, parent_column] == root_value

    if dist_column is None:
        def get_dist(node): return 1
    else:
        def get_dist(node):
            dist = df.loc[node, dist_column]
            # NaN creates an error with ete3.Tree.show()
            return 0 if isnan(dist) else dist

    for nodename, children in roll_leafwards_indices(df, root_value=root_value):
        #root, root_children = next(iter_leafwards)
        #if is_root(nodename):
        try:
            node = trees.pop(nodename)
        except KeyError:
            node = ete3.TreeNode(name=nodename, dist=get_dist(nodename))
            roots.append(node)

        for ch in children:
            if ch in trees:
                logger.error('Node name already used: %r (override).', ch)

            trees[ch] = node.add_child(name=ch, dist=get_dist(ch))
            logger.debug('node %r %r -> child %r dist[%r]=%s', nodename, node,
                         ch, dist_column, trees[ch].dist)

    return roots


def get_topo_time(df, **kwargs):
    """Set arbitrary branch lengths: a parent node is at distance 1 of the
    closest child, and all leaves are at age 0."""
    topo = pd.DataFrame(columns=['topo_age', 'topo_brlen'],
                        index=df.index)
    is_leaf = df['type'] == 'leaf'
    topo.loc[is_leaf, 'topo_age'] = 0
    for parent, children in roll_rootwards_indices(df, type_column='type', **kwargs):
        topo.loc[parent, 'topo_age'] = topo.loc[children, 'topo_age'].max() + 1
        topo.loc[children, 'topo_brlen'] = topo.loc[parent, 'topo_age'] - topo.loc[children, 'topo_age']
        logger.debug('parent %r age=%s, children %s ages=%s', parent,
                     topo.loc[parent, 'topo_age'], children,
                     topo.loc[children, 'topo_age'].tolist())

    return topo


