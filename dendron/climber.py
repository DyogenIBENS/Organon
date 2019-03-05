#!/usr/bin/env python3


#import queue
import logging
#logging.basicConfig(format='%(levelname)s:%(funcName)s:%(message)s')
logger = logging.getLogger(__name__)


def bfw_pairs(phyltree, queue=None, closest_first=True):
    """ breadth First Walk
    iterate over each branch (yields a pair : ancestor,child),
    in the order of the child's age"""

    if queue is None:  # initialize
        queue = [(None, phyltree.root)]

    if len(queue) == 0:  # terminate
        raise StopIteration

    # iterate
    parent, node = queue.pop(0)
    if parent is not None:  # do not yield root
        yield parent, node

    children = phyltree.items.get(node, [])
    sort_coeff = 1 if closest_first else -1
    children.sort(key=lambda x: sort_coeff * x[1])

    for child, _ in children:
        # append pair (parent, child)
        queue.append((node, child))

    for nextnode, nextchild in bfw_pairs(phyltree, queue, closest_first):
        yield nextnode, nextchild


def dfw_pairs(phyltree, queue=None, closest_first=False, nosinglechild=False):
    """Depth-first walk"""
    if queue is None:  # initialize
        queue = [(None, phyltree.root)]

    if len(queue) == 0:  # terminate
        raise StopIteration

    # iterate
    parent, node = queue.pop()  # != bfw
    if parent is not None:
        yield parent, node

    children = phyltree.items.get(node, [])
    if nosinglechild:
        # Skip children that have only one child (update to next child)
        for i, (ch, dist) in enumerate(children):
            nextch = phyltree.items.get(ch, [])
            while len(nextch) == 1:
                #print(".", end='')
                dist += nextch[0][1]
                children[i] = (nextch[0][0], dist)
                nextch = phyltree.items.get(nextch[0][0], [])

    sort_coeff = 1 if closest_first else -1
    children.sort(key=lambda x: sort_coeff * x[1])

    for child, _ in children:
        # append pair (parent, child)
        queue.append((node, child))

    for nextnode, nextchild in dfw_pairs(phyltree, queue, closest_first,
                                         nosinglechild):
        yield nextnode, nextchild


def dfw_pairs_generalized(tree, get_children, queue=None, include_root=False):
    """Depth-first walk"""
    if queue is None:  # initialize
        # try to find the `root` attribute of the tree
        try:
            queue = [(None, tree.root)]
        except AttributeError as e:
            logger.error("No root attribute found. Please initialize "
                          "the queue with the root: queue=[(None, root)]")
            queue = []  # will raise Stop iteration

    if len(queue) == 0:  # terminate
        raise StopIteration

    # iterate
    parent, node = queue.pop()  # != bfw
    if parent is not None or include_root:
        yield parent, node

    children = get_children(tree, node)

    # reverse the append so that the iteration follows the order of the children
    for child in reversed(children):
        # append pair (parent, child)
        queue.append((node, child))

    for nextnode, nextchild in dfw_pairs_generalized(tree, get_children, queue,
                                                     include_root):
        #print("yield %s, %s" % (nextnode, nextchild))
        yield nextnode, nextchild


def bfw_descendants(phyltree, include_leaves=False, queue=None):
    """
    Breadth First Walk

    Yields pair (parent, list of direct descendants)

    Arguments:
    - queue: list of lists of descendants to evaluate
    - include_leaves: whether to return the pair (species, empty list)

    Note: it is useful to use it in reverse (climb down the tree by gathering
    branches)

    Example : recreate the dictionary phyltree.species:

    # initialize leaves
    redo_species = {sp: frozenset([sp]) for sp in phyltree.species}
    for anc1, anc2list in reversed(list(bfw_descendants)):
        redo_species[anc1] = frozenset()
        for anc2 in anc2list:
            # we know anc2 is in the keys
            redo_species[anc1] |= redo_species[anc2]
    """

    if queue is None:  # initialize
        queue = [phyltree.root]

    if len(queue) == 0:  # terminate
        raise StopIteration

    # iterate
    parent = queue.pop(0)

    children = phyltree.items.get(parent, [])
    if include_leaves or children:
        children.sort(key=lambda x: x[1])
        descendants = [ch[0] for ch in children]

        queue.extend(descendants)

        yield parent, descendants

    for nextparent, nextdescendants in bfw_descendants(phyltree,include_leaves,
                                                       queue):
        yield nextparent, nextdescendants


def bfw_descendants_generalized(tree, get_children, include_leaves=False,
                                queue=None):
    """
    Breadth First Walk on any kind of tree, using `get_children` to recurse.

    Yields pair (parent, list of direct descendants)

    Arguments:

    - tree
    - get_children: function taking (tree, node) and returning the list of
                    children.
    - queue: list of lists of descendants to evaluate
    - include_leaves: whether to return the pair (species, empty list)

    Note: it is useful to use it in reverse (climb down the tree by gathering
    branches)
    """

    if queue is None:  # initialize
        # try to find the `root` attribute of the tree
        try:
            queue = [tree.root]
        except AttributeError as e:
            logger.error("No root attribute found. Please initialize "
                          "the queue with the root: queue=[root]")
            queue = []  # will raise Stop iteration

    if len(queue) == 0:  # terminate
        raise StopIteration

    # iterate
    parent = queue.pop(0)

    children = list(get_children(tree, parent))
    if include_leaves or children:
        queue.extend(children)
        #queue.append(descendants)

        yield parent, children

    for nextparent, nextchildren in bfw_descendants_generalized(tree,
                                        get_children, include_leaves, queue):
        yield nextparent, nextchildren


def dfw_descendants(phyltree, closest_first=False,
                    include_leaves=False, queue=None,
                    nosinglechild=False):     # trees
    """Depth-first walk. yields node and its children list"""
    if queue is None:  # initialize
        queue = [phyltree.root]

    if len(queue) == 0:  # terminate
        raise StopIteration

    # iterate
    parent = queue.pop()  # != bfw

    children = phyltree.items.get(parent, [])

    if nosinglechild:
        # Skip children that have only one child (update to next child)
        for i, (ch, dist) in enumerate(children):
            nextch = phyltree.items.get(ch, [])
            try:
                while len(nextch) == 1:
                    #print(".", end='')
                    dist += nextch[0][1]
                    children[i] = (nextch[0][0], dist)
                    nextch = phyltree.items.get(nextch[0][0], [])
            except:
                raise
            #finally:
            #    print(parent, children, nextch, file=sys.stderr)

    if include_leaves or children:
        sort_coeff = 1 if closest_first else -1
        children.sort(key=lambda x: sort_coeff * x[1])
        descendants = [ch[0] for ch in children]

        queue.extend(descendants)

        yield parent, descendants

    for nextnode, nextdescendants in dfw_descendants(phyltree, closest_first,
                                                     include_leaves, queue,
                                                     nosinglechild):
        yield nextnode, nextdescendants


def dfw_descendants_generalized(tree, get_children, include_leaves=False,
                                queue=None, copy=True):         # trees
    """Depth-first walk. yields node and its children list
    You can give any tree structure as long as you provide the appropriate
    get_children function

    Examples:
      * for a LibsDyogen.Phyltree instance:
        lambda tree, node: [x[0] for x in tree.items.get(node, [])]

      * LibsDyogen.ProteinTree instance:
        lambda tree, node: [x[0] for x in tree.data.get(node, [])]

      * ProteinTree, but yielding (node, [(descendant, dist), ...]):
        lambda tree, datum: tree.data.get(datum[0], [])
        
        but also give `queue = [(tree.root, 0)]`.
      
      * nested dictionary:
        lambda tree, node: tree.get(node, {}).keys()

      * ete3.Tree:
        lambda tree, node: node.children

      * Bio.Phylo.Tree:
        lambda tree, node: node.clades
    """
    if queue is None:  # initialize
        # try to find the `root` attribute of the tree
        try:
            queue = [tree.root]
        except AttributeError as e:
            logger.error("No root attribute found. Please initialize "
                          "the queue with the root: queue=[root]")
            queue = []  # will raise Stop iteration

    if len(queue) == 0:  # terminate
        raise StopIteration

    # iterate
    parent = queue.pop()  # != bfw

    descendants = get_children(tree, parent)
    if copy:
        descendants = list(descendants)  # copying is important!
    if include_leaves or descendants:
        yield parent, descendants
        queue.extend(reversed(descendants))

    for nextnode, nextdescendants in \
            dfw_descendants_generalized(tree, get_children,include_leaves, queue, copy):
        yield nextnode, nextdescendants


def dfw_lineage_generalized(tree, get_children, queue=None):
    """Unfinished"""
    raise NotImplementedError
    ### TODO
    if queue is None:  # initialize
        # try to find the `root` attribute of the tree
        try:
            queue = [[tree.root]]
        except AttributeError as e:
            logger.error("No root attribute found. Please initialize "
                          "the queue with the root: queue=[root]")
            queue = []  # will raise Stop iteration

    if len(queue) == 0:  # terminate
        raise StopIteration

    # iterate
    prev_lineage ## TODO = queue.pop()  # != bfw

    descendants = list(get_children(tree, parent)) # copying is important!
    if include_leaves or descendants:
        queue.extend(reversed(descendants))
        yield parent, descendants

    for nextnode, nextdescendants in dfw_lineage_generalized(tree,
                                        get_children,include_leaves, queue):
        yield nextnode, nextdescendants

#def iter_all_descending_paths(tree, get_children, root):
#    """Iterate over all possible descending paths (list of directly connected
#    nodes, from all ancestors to all of its descendants, leaves and internals).
#    
#    Will yield the terminal paths first."""
#    straight_lines = [[root]]
#    for node, descendants in dfw_descendants_generalized(tree, get_children,
#                                                         include_leaves=True,
#                                                         queue=[root]):
#        if descendants:
#        #try:
#            straight_lines[-1].append(descendants.pop())
#        else:
#        #except IndexError:
#            # We reached a leaf. Let's yield all paths on the corresponding
#            # straight line.
#            yield straight_lines[-1]
#            #N = len(straight_lines[0])
#            #for start, end in combinations(range(N), 2):
#            #    yield straight_lines[0][start:end]
#
#            straight_lines.pop()
#
#        for child in reversed(descendants):
#            straight_lines.append([node, child])


def rev_dfw_descendants(*args, **kwargs):
    return reversed(list(dfw_descendants_generalized(*args, **kwargs)))


def rev_bfw_descendants(*args, **kwargs):
    return reversed(list(bfw_descendants_generalized(*args, **kwargs)))


def iter_all_ancestor_descendant_pairs(tree, get_children, root=None):
    startqueue = None if root is None else [(None, root)]

    for _, startnode in dfw_pairs_generalized(tree, get_children, startqueue,
                                                 include_root=True):
        nextstartqueue = [(None, startnode)]
        for _, endnode in dfw_pairs_generalized(tree, get_children, nextstartqueue):
            yield startnode, endnode


def iter_all_paths_fromroot(tree, get_children, path):
    """path should be initialized by [root]"""
    assert isinstance(path, list)
    for child in get_children(tree, path[-1]):
        nextpath = path + [child]
        yield nextpath
        for nnextpath in iter_all_paths_fromroot(tree, get_children, nextpath):
            yield nnextpath


def iter_leaf_paths(tree, get_children, path):
    assert isinstance(path, list)
    children = get_children(tree, path[-1])
    if not children:
        yield path
    else:
        for child in children:
            nextpath = path + [child]
            for leafpath in iter_leaf_paths(tree, get_children, path=nextpath):
                yield leafpath


def iter_all_paths(tree, get_children, root):
    for rootpath in iter_all_paths_fromroot(tree, get_children, [root]):
        # Yield in a specific order: start from the included chunks to the larger
        # chunks.
        for i in range(len(rootpath)-1, 0, -1):
            yield rootpath[(i-1):]


def iter_leaves(tree, get_children, queue=None):
    for node, children in dfw_descendants_generalized(tree, get_children,
                                                      queue=queue,
                                                      include_leaves=True):
        if not children:
            yield node


def iter_distleaves(tree, get_data, root=None, root_dist=0):
    """Iterate over the pairs (leaf, distance from the root).
    param: `get_data` takes a tuple (tree, node) and return a
               tuple (child, dist) for each child of the given node (as a list)."""
    if root is None:
        try:
            root = tree.root
        except AttributeError:
            root = tree
    for leafpath in iter_leaf_paths(tree, get_data, [(root, root_dist)]):
        leaf = leafpath[-1][0]
        leafdist = sum(dist for node, dist in leafpath[1:])
        yield leaf, leafdist

