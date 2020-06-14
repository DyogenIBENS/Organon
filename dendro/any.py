
"""Module for wrapping common functions required by other dendro modules such as 
dendro.bates.

Defines the get_item, get_children, get_root functions for instances of
ete3.Tree, Bio.Phylo.Newick.Tree, LibsDyogen.myPhylTree, LibsDyogen.myProteinTree.
"""

import warnings


class TreeMethod(object):
    # Not meant to be used directly
    
    @staticmethod
    def get_dist(tree, node):
        return node.dist
    
    @staticmethod
    def get_root(tree):
        return tree.root  # Not ete3 but more often implemented like this.
    
    @staticmethod
    def get_label(tree, node):
        return node

    @staticmethod
    def get_children(tree, node):
        return node.children  # ete3

    @classmethod
    def get_items(cls, tree, nodedist):
        return [(child, cls.get_dist(tree, child))
                for child in cls.get_children(tree, nodedist[0])]


class nodebased(TreeMethod):
    @staticmethod
    def get_label(tree, node):
        return node.name


class itembased(TreeMethod):
    """ For trees based on collections of items (child, dist).
    
    Simply overwrite the `get_items` method.

    Not meant to be used directly.
    """
    @classmethod
    def get_children(cls, tree, node):
        return [child for child, _ in cls.get_items(tree, (node, None))]


class ete3(nodebased):
    
    @staticmethod
    def get_root(tree):
        return tree.get_tree_root()
    
    @staticmethod
    def copy_children(tree, node):
        # This should in general be favored over 'get_children', to avoid unexpected behavior.
        return tree.get_children()

    @staticmethod
    def set_items(tree, node, new_items):
        for child, dist in new_items:
            node.add_child(child, dist=dist)

    @staticmethod
    def set_children(tree, node, new_children):
        # If new_children is the same list object than node.children, weird behavior expected.
        # therefore, just re-init a new list for the children:
        node.children = []
        for child in new_children:
            node.add_child(child)

    @staticmethod
    def print_newick(tree, stream=None, root=None, format=1, format_root_node=True, **kwargs):
        if root is not None:
            tree = tree&root
        print(tree.write(outfile=None, format=format,
                         format_root_node=format_root_node, **kwargs), file=stream)

    @classmethod
    def set_dist(cls, tree, node, dist):
        node.dist = dist

class BioPhylo(nodebased):
    
    @staticmethod
    def get_dist(tree, node):
        return node.branch_length
    
    @staticmethod
    def get_children(tree, node):
        return node.clades

    @staticmethod
    def set_children(tree, node, children):
        node.clades = children


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

    @staticmethod
    def set_items(tree, nodedist, items):
        tree.data[nodedist[0]] = items

    @staticmethod
    def get_label(tree, node):
        return tree.info[node].get('family_name')


class myPhylTree(itembased):

    @staticmethod
    def get_dist(tree, node):
        if node == tree.root:
            return getattr(tree, 'rootlength', None)
        else:
            return tree.parent.get(node, (None, None))[1]
    
    @staticmethod
    def get_items(tree, nodedist):
        return tree.items.get(nodedist[0], [])

    @staticmethod
    def set_items(tree, nodedist, items):
        tree.items[nodedist[0]] = items

    @staticmethod
    def print_newick(tree, stream=None, root=None, commonnames=True, symbols=True, **kwargs):
        tree.printNewick(stream, root=root, commonnames=commonnames, symbols=symbols, **kwargs)

    @classmethod
    def set_dist(cls, tree, node, dist):
        if node == cls.get_root(tree):
            tree.rootlength = dist
        else:
            parent = tree.parents[node]
            cls.set_items(tree, (parent, None),
                    [(child, (dist if child==node else chdist)) for
                        child,chdist in cls.get_items(tree, (parent, None))])


def skip_set_children(tree, node, children):
    pass


class BioNexus(TreeMethod):
    @staticmethod
    def get_root(tree):
        return tree.node(tree.root)

    @staticmethod
    def get_label(tree, node):
        return node.data.taxon

    @staticmethod
    def get_children(tree, node):
        """Take an instance of Bio.Nexus.Trees.Tree and a Bio.Nexus.Nodes.Node"""
        return [tree.node(child_id) for child_id in node.get_succ()]

    @staticmethod
    def get_dist(tree, node):
        return node.data.branchlength


methodchoice = {'phyltree': myPhylTree,
                'proteintree': myProteinTree,
                'ete3': ete3,
                'biophylo': BioPhylo,
                'bionexus': BioNexus}
