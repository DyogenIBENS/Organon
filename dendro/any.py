
"""Module for wrapping common functions required by other dendro modules such as 
dendro.bates.

Defines the get_item, get_children, get_root functions for instances of
ete3.Tree, Bio.Phylo.Newick.Tree, LibsDyogen.myPhylTree, LibsDyogen.myProteinTree.
"""

import warnings


class default(object):
    # Not meant to be used directly
    
    @staticmethod
    def get_dist(tree, node):
        return node.dist
    
    @staticmethod
    def get_root(tree):
        return tree.root  # Not ete3 but more often implemented like this.
    
    @staticmethod
    def get_label(node):
        return node

    @staticmethod
    def get_children(tree, node):
        return node.children  # ete3

    @classmethod
    def get_items(cls, tree, nodedist):
        return [(child, cls.get_dist(tree, child))
                for child in cls.get_children(tree, nodedist[0])]


class nodebased(default):
    @staticmethod
    def get_label(node):
        return node.name


class itembased(default):
    """ For trees based on collections of items (child, dist).
    
    Simply overwrite the `get_items` method.

    Not meant to be used directly.
    """
    @classmethod
    def get_children(cls, tree, node):
        return [child for child, dist in cls.get_items(node, [])]


class ete3(nodebased):
    
    @staticmethod
    def get_root(tree):
        return tree.get_tree_root()
    
    @staticmethod
    def copy_children(tree, node):
        return tree.get_children()

    @staticmethod
    def set_items(tree, node, new_items):
        for child, dist in new_items:
            node.add_child(child, dist=dist)

    @staticmethod
    def set_children(tree, node, new_children):
        for child in new_children:
            node.add_child(child)


class BioPhylo(nodebased):
    
    @staticmethod
    def get_dist(tree, node):
        return node.branch_length
    
    @staticmethod
    def get_children(tree, node):
        return node.clades


class myProteinTree(itembased):
    
    @staticmethod
    def get_dist(tree, node):
        warnings.warn('Unefficient strategy for myProteinTree, use get_items instead')
        for items in tree.data.values():
            for child, dist in items:
                if child == node:
                    return dist
    
    @staticmethod
    def get_items(tree, nodedist):
        return tree.data.get(nodedist[0], [])


class myPhylTree(itembased):

    @staticmethod
    def get_dist(tree, node):
        if node == tree.root:
            return getattr(tree, rootdist, None)
        else:
            return tree.parent.get(node, (None, None))[1]
    
    @staticmethod
    def get_items(tree, nodedist):
        return tree.items.get(nodedist[0], [])
