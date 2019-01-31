
"""Module to convert Python tree objects between different classes.

(ete3, LibsDyogen, scipy.hclust.linkage...)"""

import numpy as np
from climber import rev_dfw_descendants, dfw_descendants_generalized
import logging
logger = logging.getLogger(__name__)


def get_data(tree, nodedist):
    return tree.data.get(nodedist[0], [])


def ete3_to_linkage(tree):
    ### PB: distances are identical on each side of any fork.
    Z = []
    node_is = {}
    leaf_i = 0
    for node_dist, children_dists in \
                                rev_dfw_descendants(tree,
                                                get_data,
                                                include_leaves=True,
                                                queue=[(tree.root, tree.dist)]):
        if not children_dists:
            node_is[node_dist[0]] = leaf_i
            leaf_i += 1
            continue

        children_is = [node_is.pop(ch) for ch,_ in children_dists]
        Z.append([])

        
def ProtTree_to_ete3(prottree):
    import ete3
    tree = ete3.Tree(name=prottree.root)
    current_nodes = [tree]
    while current_nodes:
        current = current_nodes.pop()
        current.add_features(**prottree.info[current.name])
        for child, dist in prottree.data.get(current.name, []):
            current_nodes.append(current.add_child(name=child, dist=dist))
    return tree

def PhylTree_to_ete3(phyltree):
    return phyltree.to_ete3()

def ete3_to_ete3(tree):
    return tree


converterchoice = {'PhylTree': {'Ete3': PhylTree_to_ete3},
                   'ProtTree': {'Ete3': ProtTree_to_ete3},
                   'Ete3': {'Ete3': ete3_to_ete3}}
for from_, todict in list(converterchoice.items()):
    for to_, func in list(todict.items()):
        todict[to_.lower()] = func
    converterchoice[from_.lower()] = todict
