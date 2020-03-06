
"""Module to convert Python tree objects between different classes.

(ete3, LibsDyogen, scipy.hclust.linkage...)"""


from .bates import rev_dfw_descendants


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
    tree = ete3.Tree(name=prottree.root, dist=0)
    current_nodes = [tree]
    while current_nodes:
        current = current_nodes.pop()
        current.add_features(**prottree.info[current.name])
        for child, dist in prottree.data.get(current.name, []):
            current_nodes.append(current.add_child(name=child, dist=dist))
    return tree


def PhylTree_to_ete3(phyltree, nosinglechild=False):
    # TODO: do not import ete3 here, just try to use it and raise error
    import ete3
    tree = ete3.Tree(name=phyltree.root, dist=getattr(phyltree, 'rootlength', 0))
    tree.add_features(treename=getattr(phyltree, 'name', ''))
    current_nodes = [tree]
    while current_nodes:
        current = current_nodes.pop()
        current.add_features(age=phyltree.ages.get(current.name),
                             commonNames=phyltree.commonNames.get(current.name),
                             fileName=phyltree.fileName.get(current.name),
                             indBranch=phyltree.indBranches.get(current.name),
                             indName=phyltree.indNames.get(current.name),
                             )
                             #officialName=phyltree.officialName.get(current.name),
        childlist = phyltree.items.get(current.name, [])
        if not childlist:
            current.add_features(Esp2X=(current.name in phyltree.lstEsp2X),
                                 Esp6X=(current.name in phyltree.lstEsp6X),
                                 EspFull=(current.name in phyltree.lstEspFull))
            # dicGenes
            # dicGenomes
        while nosinglechild and (len(childlist) == 1):
            current.name = childlist[0][0]
            current.dist += childlist[0][1]
            childlist = phyltree.items.get(current.name, [])

        for child, dist in childlist:
            current_nodes.append(current.add_child(name=child, dist=dist))
    return tree


def PhylTree_draw(phyltree, images=None):
    """Draw tree in a gui window using Ete3"""
    import ete3
    ptree = PhylTree_to_ete3(phyltree)
    ts = ete3.TreeStyle()
    ts.scale = 2
    ts.show_branch_length = True
    nameface = ete3.faces.AttrFace("name", fsize=12)
    def mylayout(node):
        if node.is_leaf():
            # add species image
            pass
        else:
            ete3.faces.add_face_to_node(nameface, node, column=0,
                                        position='branch-top')
    ptree.show(layout=mylayout, tree_style=ts, name='Species tree')


def ete3_to_ete3(tree):
    return tree



converterchoice = {'PhylTree': {'Ete3': PhylTree_to_ete3},
                   'ProtTree': {'Ete3': ProtTree_to_ete3},
                   'Ete3': {'Ete3': ete3_to_ete3}}
for from_, todict in list(converterchoice.items()):
    for to_, func in list(todict.items()):
        todict[to_.lower()] = func
    converterchoice[from_.lower()] = todict
