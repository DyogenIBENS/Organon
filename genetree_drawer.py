#!/usr/bin/env python3

import re
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import ete3

import LibsDyogen.myPhylTree as PhylTree
import LibsDyogen.myProteinTree as ProteinTree

from glou_duphist import bfw_descendants_generalized, ladderize
from codeml.select_leaves_from_specieslist import convert_gene2species


ANCGENE2SP = re.compile(r'([A-Z][A-Za-z_.-]+)ENS')

TESTTREE = "/users/ldog/glouvel/ws2/DUPLI_data85/alignments/ENSGT00810000125388/subtrees2/RodentiaENSGT00810000125388.A.a.a.c.a.b.nwk"

### Matplotlib graphical parameters ###
grey10 = '#1a1a1a'
mpl.rcParams['text.color'] = grey10
mpl.rcParams['axes.edgecolor'] = grey10
mpl.rcParams['axes.labelcolor'] = grey10
mpl.rcParams['axes.spines.top'] = False
mpl.rcParams['axes.spines.right'] = False
mpl.rcParams['xtick.color'] = grey10
mpl.rcParams['ytick.color'] = grey10
mpl.rcParams['grid.color'] = grey10
mpl.rcParams['patch.edgecolor'] = grey10
mpl.rcParams['boxplot.flierprops.markeredgecolor'] = grey10
mpl.rcParams['boxplot.capprops.color'] = grey10
mpl.rcParams['legend.facecolor'] = '#777777'
mpl.rcParams['legend.framealpha'] = 0.2
#mpl.rcParams['legend.edgecolor'] = grey10
mpl.rcParams['savefig.facecolor'] = 'none'
mpl.rcParams['savefig.frameon'] = False  #background frame transparent
#mpl.rcParams['savefig.transparent'] = True # all background transparent
                                            # (including ggplot2 style)
#mpl.style.use('ggplot')

def get_taxon(node):
    """from a gene name in my newick gene trees, find the taxon:
        either:
            - node is a leaf (e.g ENSMUSG00...)
            - node is internal (e.g Mus.musculusENSGT...)"""
    if node.is_leaf():
        taxon = convert_gene2species(node.name)
    else:
        try:
            taxon = ANCGENE2SP.match(node.name).group(1).replace('.', ' ')
        except AttributeError:
            raise RuntimeError("Can not match species name in %r" % node.name)
    return taxon


def get_taxa_set(etegenetree):
    taxa_set = set()
    for node in etegenetree.traverse():
        taxon = get_taxon(node)
        taxa_set.add(taxon)
    return taxa_set


def bfw_phylsubtree(phyltree, taxa):
    """Return an iterator that will progress on the tree from leaves to root.
    
    It also rotates forks such that the displayed tree will show branches
    that are most distant from the root (in number of nodes) on one side.

    Taxa returned by the iterator are space-separated.

    Use breadth-first iteration.
    """
    root, subtree = phyltree.getSubTree(taxa)

    # reorder branches in a visually nice manner:
    ladderize(subtree, root)
    get_children = lambda tree, node: tree.get(node, [])
    bfw = bfw_descendants_generalized(subtree, get_children,
                                      include_leaves=False, queue=[root])
    return reversed(list(bfw))


def iter_species_coords(phyltree, taxa):
    """Assign a pair x,y of coordinates for each node in the tree.
    Yield (parent name, parent_xy, child name, child_xy).
    
    Just the topology. Branch lengths are ignored.
    """
    coords = {}
    y0 = 0
    for parent, children in bfw_phylsubtree(phyltree, taxa):
        children_xs = []
        children_ys = []
        for child in children:
            child_xy = coords.get(child)
            if not child_xy:
                x = 0
                y = y0
                y0 -= 1
                coords[child] = (x, y)
            else:
                x, y = child_xy
            children_xs.append(x)
            children_ys.append(y)

        # Along the X axis: move one step to the left
        # Along the Y axis: take the middle of the children Y coordinates.
        parent_x = min(children_xs) - 1
        parent_y = (min(children_ys) + max(children_ys)) / 2
        coords[parent] = (parent_x, parent_y)

        for child in children:
            yield parent, coords[parent], child, coords[child]


class GenetreeDrawer(object):
    """Draw a gene tree inside a species tree"""

    ensembl_version = 85
    phyltreefile = "/users/ldog/glouvel/ws_alouis/GENOMICUS_SVN/data{0}/" \
                   "PhylTree.Ensembl.{0}.conf"
    
    def __init__(self):
        self.phyltree = None
        pass

    def load_phyltree(self, phyltreefile=None, ensembl_version=None):
        phyltreefile = phyltreefile if phyltreefile else self.phyltreefile
        ensembl_version = ensembl_version if ensembl_version else \
                          self.ensembl_version
        self.phyltree = PhylTree.PhylogeneticTree(phyltreefile.format(
                                                            ensembl_version))

    def load_reconciled_genetree(self, filename, format=1):
        """Load gene tree with all species nodes present.
        newick format with internal node labelling"""
        self.genetree = ete3.Tree(filename, format=format)
        # add only the meaningful taxa (not those with one child and age = 0)
        alldescendants = self.phyltree.allDescendants[get_taxon(self.genetree)]
        self.taxa = set()
        for taxon in alldescendants:
            if not (len(self.phyltree.items.get(taxon, [])) == 1 and
                    self.phyltree.ages[taxon] == 0):
                self.taxa.add(taxon)


    def draw_species_tree(self, equal_angles=False, branch_width=0.8):
        """Init figure + draw branches of the species tree.
        
        branch_width: proportion of vertical space between two branches taken
                      by the branch polygon."""
        # TODO: implement equal_angles
        self.species_coords   = {}
        self.species_branches = {} 

        self.fig, ax0 = plt.subplots(frameon=False) # set figsize later
        #ax0 = self.fig.add_axes([0.1,0.1,0.9,0.9]) #, adjustable='box-forced')
        #ax0 = self.fig.add_axes([0.1, 0.1, 0.9, 0.9])
        ax0.axis('off')

        ymin = 0
        for parent, (px, py), child, (cx, cy) in iter_species_coords(self.phyltree,
                                                                     self.taxa):
            self.species_branches[child] = (parent, cx - px, cy - py)

            self.species_coords[child] = (cx, cy)
            if cy < ymin:
                ymin = cy
            coords = np.array([(px, py),
                               (cx, cy),
                               (cx, cy - branch_width),
                               (px, py - branch_width)])
            ax0.add_patch(patches.Polygon(coords, 
                                          facecolor='#e5e5e5',
                                          edgecolor='#e5e5e5'))
            ax0.text(cx, cy, child, ha='center', fontsize=10,
                     fontstyle='italic', family='serif')

        self.species_coords[parent] = (px, py)
        ax0.text(px, py, parent, ha='center', fontsize=10, fontstyle='italic',
                 family='serif')
        
        ax0.set_xlim(px, 1)
        ax0.set_ylim(ymin - 1, 1)
        self.ax0 = ax0


    def set_gene_coords(self):
        self.gene_coords = {} # {species: [gene list]}
        #self.dup_branchings = [] # (genename, x, x_child1, x_child2, y_child1, y_child2)
        #self.spe_branchings = []
        self.branchings = {} # {genename: [species, x, y, dup/spe, (children)], ...]
        
        # temporary tree structure holding nodes that are inbetween 2 speciations
        interspecies_trees = {}

        # x and y are relative to the branch considered. (between 0 and 1)

        for node in self.genetree.traverse('postorder'):
            taxon = get_taxon(node)
            taxon_gene_coords = self.gene_coords.setdefault(taxon, [])

            if node.is_leaf():
                taxon_gene_coords.append(node.name)
                interspecies_trees[node.name] = {'taxon': taxon,
                                                 'ndup': 0,
                                                 'x': 1,
                                                 'y': len(taxon_gene_coords) - 1}
            else:
                children_taxa = set(interspecies_trees[ch.name]['taxon'] for ch
                                    in node.children)
                if len(children_taxa) == 1 and len(node.children) > 1:
                    # It is a duplication
                    children_ys   = []
                    children_ndup = []
                    for ch in node.children:
                        children_ys.append(interspecies_trees[ch.name]['y'])
                        children_ndup.append(interspecies_trees[ch.name]['ndup'])

                    node_y    = (min(children_ys) + max(children_ys)) / 2
                    node_ndup = max(children_ndup) + 1

                    interspecies_trees[node.name] = {
                            'taxon': taxon,
                            'ndup': node_ndup,
                            'x': 1, # modified later
                            'y': node_y}
                else:
                    # It is a speciation
                    taxon_gene_coords.append(node.name)

                    node_y = len(taxon_gene_coords) - 1
                    interspecies_trees[node.name] = {
                            'taxon': taxon,
                            'ndup': 0,
                            'x': 0,
                            'y': node_y}
                    self.branchings[node.name] = [taxon, 0, node_y, 'spe',
                                                  [ch.name for ch in
                                                      node.children]]

                    for ch in node.children:
                        n_steps = interspecies_trees[ch.name]['ndup'] + 1

                        delta_x = 1 / n_steps

                        # climb up this subtree until next speciation
                        nextnodes = [ch]
                        while nextnodes:
                            nextnode = nextnodes.pop(0)
                            nextnode_ndup = interspecies_trees[nextnode.name]['ndup']
                            if nextnode_ndup > 0:
                                nextnode_x = 1 - nextnode_ndup * delta_x
                                nextnode_y = interspecies_trees[nextnode.name]['y']
                                nextnode_taxon = interspecies_trees[nextnode.name]['taxon']
                                nextnode_ch = [nnch.name for nnch in nextnode.children]

                                self.branchings[nextnode.name] = [
                                                                nextnode_taxon,
                                                                nextnode_x,
                                                                nextnode_y,
                                                                'dup',
                                                                nextnode_ch]
                                nextnodes.extend(nextnode.children)

                            interspecies_trees.pop(nextnode.name)
        print(interspecies_trees)

                                                                
    def draw_gene_tree(self, branch_width=0.8):
        self.real_gene_coords = {}

        for node in self.genetree.traverse('postorder'):
            genename = node.name
            print(genename)
            
            if node.is_leaf():
                species = convert_gene2species(genename)
                real_x, real_y = self.species_coords[species]
                pos_list = self.gene_coords[species]
                pos = pos_list.index(genename) + 1
                real_y -= branch_width * pos/(len(pos_list) + 1)
            else:
                species, rel_x, rel_y, event, children = self.branchings[genename]

                children_ys = [self.real_gene_coords[ch][1] for ch in children]

                if event == 'dup':
                    parent_sp, Dx, Dy = self.species_branches[species]
                    parent_x, parent_y = self.species_coords[parent_sp]
                    real_x = parent_x + Dx * rel_x
                    real_y = parent_y + Dy * rel_x - branch_width * rel_y
                    nodecolor = 'red'

                    delta_y = (children_ys[0] - children_ys[-1]) / 2

                    for ch, delta in zip(children, [1, -1]):
                        ch_real_x, ch_real_y = self.real_gene_coords[ch]
                        self.ax0.plot((real_x, ch_real_x),
                                      (real_y + delta * delta_y, ch_real_y),
                                      color=grey10)
                    
                    self.ax0.plot((real_x, real_x),
                                  (real_y - delta_y, real_y + delta_y),
                                  color=grey10)

                else:
                    real_x, real_y = self.species_coords[species]
                    pos_list = self.gene_coords[species]
                    pos = pos_list.index(genename) + 1
                    real_y -= branch_width * pos/(len(pos_list) + 1)
                    nodecolor = 'blue'
                    children_xs = [self.real_gene_coords[ch][0] for ch in children]

                    for ch_real_x, ch_real_y in zip(children_xs, children_ys):
                        self.ax0.plot((real_x, ch_real_x), (real_y, ch_real_y),
                                      color=grey10)
                
                self.ax0.plot((real_x,), (real_y,), 'o', color=nodecolor)
                self.ax0.text(real_x, real_y, species, fontsize='xx-small',
                                color=nodecolor)

            self.real_gene_coords[genename] = (real_x, real_y)



if __name__ == '__main__':
    gd = GenetreeDrawer()
    gd.load_phyltree()
    gd.load_reconciled_genetree(TESTTREE)
    gd.draw_species_tree()
    gd.set_gene_coords()
    gd.draw_gene_tree()

