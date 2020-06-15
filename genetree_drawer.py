#!/usr/bin/env python3


"""
Draw a gene tree inside a species tree.

"""
# USAGE:
# 
#     ./genetree_drawer.py <genetreefile>
# 
# ARGUMENTS:
#   - outputfile:   pdf file, or '-'. If '-', will use Qt to display the figure.
#   - genetreefile: must be a genetree (nwk format with internal nodes labelling)
#                   reconciled with species tree.

import sys
import os
import os.path as op
import re
try:
    import argparse_custom as argparse # Less verbose help message
except ImportError:
    import argparse
from itertools import zip_longest

import numpy as np
import matplotlib as mpl
#mpl.use('Agg')
if __name__=='__main__': mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import patches
from matplotlib import lines
from matplotlib import patheffects
#from matplotlib.transforms import Affine2D
from matplotlib.backends.backend_pdf import PdfPages, FigureCanvas
from matplotlib.path import Path
MOVETO, CURVE3, LINETO = Path.MOVETO, Path.CURVE3, Path.LINETO
import ete3

import LibsDyogen.myPhylTree as PhylTree
#import LibsDyogen.myProteinTree as ProteinTree

from dendro.bates import rev_dfw_descendants
from dendro.sorter import ladderize, pyramid
from dendro.reconciled import get_taxon, get_taxon_treebest, infer_gene_event

import logging
logger = logging.getLogger(__name__)


ENSEMBL_VERSION = 85
PHYLTREEFILE = "/users/ldog/glouvel/GENOMICUS{0}/PhylTree.Ensembl.{0}.conf"
ANCGENE2SP = re.compile(r'([A-Z][A-Za-z_.-]+)ENS') # NOT USED
FIGSIZE = re.compile(r'(?<=:)(?:(\d+(?:\.\d+)?)x(\d+(?:\.\d+)?)|a[0-4])$', re.I)
PAPERSIZE = {'a4': (8.27, 11.7),
             'a3': (11.7, 16.54),
             'a2': (16.54, 23.4),
             'a1': (23.4, 33.1),
             'a0': (33.1, 46.8)}
fontsizes=['xx-small', 'x-small', 'small', 'medium', 'large', 'x-large', 'xx-large']
basefontsize='medium'  #'x-small'

### Matplotlib graphical parameters ###
grey10 = '#1a1a1a'
mpl.rcParams['figure.figsize'] = (8.27, 11.7) # a4
mpl.rcParams['lines.linewidth'] = 0.8
mpl.rcParams['lines.markersize'] = 3
mpl.rcParams['text.color'] = grey10
mpl.rcParams['axes.edgecolor'] = grey10
mpl.rcParams['axes.labelcolor'] = grey10
mpl.rcParams['axes.spines.top'] = False
mpl.rcParams['axes.spines.right'] = False
mpl.rcParams['xtick.color'] = grey10
mpl.rcParams['ytick.color'] = grey10
mpl.rcParams['grid.color'] = grey10
mpl.rcParams['patch.edgecolor'] = grey10
mpl.rcParams['patch.linewidth'] = 0.8
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

def get_common_name(phyltree, latin_name):
    """Get the common species name from the latin one"""
    for name in phyltree.commonNames[latin_name]:
        if name != latin_name and \
                not isinstance(name, int) and \
                not '_' in name:
            return name
    return None



def walk_phylsubtree(phyltree, taxa):
    """Return an iterator that will progress on the tree from leaves to root.
    
    It also rotates forks such that the displayed tree will show branches
    that are most distant from the root (in number of nodes) on one side.

    Taxa returned by the iterator are space-separated.
    """
    if "root" in taxa:
        logger.debug('"root" is in taxa')
        root = "root"
        taxa.remove("root")
        oldroot, subtree = phyltree.getSubTree(taxa)
        subtree[root] = [oldroot]
    else:
        root, subtree = phyltree.getSubTree(taxa)

    # reorder branches in a visually nice manner:
    get_children = lambda tree, node: tree.get(node, [])
    ladderize(subtree, root, get_children, heavy_on_top=True)
    #pyramid(subtree, oldroot, get_children)
    return rev_dfw_descendants(subtree, get_children,
                              include_leaves=False, queue=[root])


def iter_species_coords(phyltree, taxa, angle_style=0, ages=False):
    """Assign a pair x,y of coordinates for each node in the tree.
    Yield (parent name, parent_xy, child name, child_xy).
    
    - angle_style: 0 to 5, see parser help
    - ages:
        True: use real species ages
        'log': use log10(1 + species age)
        else: just the topology, real branch lengths are ignored.

    X-axis is left-to-right; older nodes on the left, leaves at 0;
    Y-axis is bottom-to-top; leaves encountered first are at the top (top at 0). 
        internal nodes encountered first are at the bottom .....
    """
    # Store the graph coords (x,y) of a node, and its weight w (number of leaves)
    coords = {}
    y0 = 0
    for parent, children in walk_phylsubtree(phyltree, taxa):
        children_xs = []
        children_ys = []
        children_ws = []
        for child in children:
            child_xyw = coords.get(child)
            if not child_xyw:
                x = 0 # x = phyltree.ages[child]
                y = y0
                y0 -= 1
                w = 1
                coords[child] = (x, y, w)
            else:
                x, y, w = child_xyw
            children_xs.append(x)
            children_ys.append(y)
            children_ws.append(w)

        parent_w = sum(children_ws)
        # Along the X axis: move one step to the left
        if ages:
            if parent == 'root':
                parent_x = max(children_xs) * 1.05
            else:
                parent_x = - phyltree.ages[parent]
            if ages == 'log':
                parent_x == - np.log10(1 + parent_x)
        else:
            parent_x = min(children_xs) - 1

        if angle_style == 1:
            # average inversely-weighted by the number of descendant leaves
            #parent_y = sum(cy/cw for cy,cw in zip(children_ys, children_ws))
            #parent_y /= sum(1/cw for cw in children_ws)
            
            # Almost the same result
            # Reverse weights to make the parent closer to the lightest node.
            parent_y = sum(cy*cw for cy,cw in zip(children_ys, reversed(children_ws)))
            parent_y /= parent_w
        elif angle_style == 2:
            # Along the Y axis: take the average of the children Y coordinates.
            parent_y = sum(children_ys) / len(children_ys)
        elif angle_style == 4:
            parent_y = sum(cy*cw for cy,cw in zip(children_ys, children_ws))
            parent_y /= parent_w

        else:
            # TODO: dx < 0 ?
            dy = max(children_ys) - min(children_ys)
            if dy == 0: dy = 1 # when only one child
            dx = max(children_xs) - min(children_xs)

            if angle_style == 0:
                step = dy / (2+dx) if len(children)>1 else 0
                #parent_y = max(children_ys) - step

                # The y-step is biased. Add the smallest step to the furthest
                # child along x.
                # (because the tree is ladderized, it's the bottom one).
                parent_y = min(children_ys) + step
            elif angle_style == 3 and dy>dx:
                # branches are all at 45 degrees
                #step from the x-closest node. dy and dx are positive.
                step = (dy - dx)/2.
                #step = dy - np.sqrt(dx*dy) if dx else dy/2.
                parent_x = min(children_xs) - step
                closest_x = children_xs.index(min(children_xs))
                parent_y = children_ys[closest_x]
                if children_ys[closest_x] < max(children_ys):  # child is lower.
                    parent_y += step
                else:  # child is higher.
                    parent_y -= step
                    #assert min(parent_x - cx for cx in children_xs) \
                    #        == min(parent_y - cy for cy in children_ys)
                stepf = max(children_xs) - parent_x
                logger.debug(
                    "%s-%s: dx=%g, dy=%g, step=%g, stepf=%g; Xs=%s, Ys=%s, closest_x=%d, parent_y=%g",
                    parent, '.'.join(children), dx, dy, step, stepf, children_xs, children_ys, closest_x, parent_y)
                #logger.debug("closest_x=%d")
                assert len(children_ys)==1 or (min(children_ys) <= parent_y <= max(children_ys))
                #assert stepf + step == dy, "%g+%g != %g" % (stepf, step, dy)
                #assert stepf - step == dx, "%g-%g != %g" % (stepf, step, dx)
                #assert not dx or step/dy == step/(step+dx), \
                #        "step=%g, dy=%g, dx=%g (%s-%s) Xs=%s, Ys=%s" % (
                #        step, dy, dx, parent, child, children_xs, children_ys)
            elif angle_style == 5 or dy <= dx:
                # equal angles
                step = dy / (2+dx)
                parent_y = max(children_ys) - step
            else:
                raise ValueError("Invalid angle_style value: %r" % angle_style)

        coords[parent] = (parent_x, parent_y, parent_w)

        #for child in children:
        #    yield parent, coords[parent], child, coords[child]
        yield parent, coords[parent], [(child, coords[child]) for child in children]



class GenetreeDrawer(object):
    """Draw a gene tree inside a species tree"""

    ensembl_version = ENSEMBL_VERSION
    phyltreefile = "/users/ldog/glouvel/ws_alouis/GENOMICUS_SVN/data{0}/" \
                   "PhylTree.Ensembl.{0}.conf"
    
    def __init__(self, phyltreefile=None, ensembl_version=None,
                 colorize_clades=None, commonname=False, latinname=False,
                 angle_style=0, ages=False, internal=None, treebest=False, show_cov=False, debug=False):
        """Options:
            - colorize_clades: grouping of species names to colorize
            - commonname: display the common english name of the species
            - latinname: display the scientific name of the species
            - treebest: whether input reconciled gene trees are formatted with
                        TreeBest tags (B=bootstrap score, S=species,
                        D="Y"/"N" for duplication)
            - show_cov: add a grey gradient to distinguish species genomes with
                        good (>6X), medium (6X) or bad (2X) coverage.
        """
        self.debug = debug
        if ensembl_version: self.ensembl_version = ensembl_version
        if phyltreefile: self.phyltreefile = phyltreefile
        self.get_taxon = get_taxon_treebest if treebest else get_taxon
        self.commonname = commonname
        self.latinname = latinname
        self.angle_style = angle_style
        self.ages = ages
        self.drawn_count = 0

        self.phyltree = PhylTree.PhylogeneticTree(self.phyltreefile.format(
                                                        self.ensembl_version))
        self.internal = set() if internal is None else internal

        # add legend elements for coverage information
        self.show_cov = show_cov
        if show_cov:
            self.legend_coverage = [
                lines.Line2D([], [], alpha=1, linestyle=' ', label='high coverage'),
                lines.Line2D([], [], alpha=0.6, linestyle=' ', label='6X coverage'),
                lines.Line2D([], [], alpha=0.3, linestyle=' ', label='2X coverage')
               ]

        # add legend elements for the clade information
        self.colorize_species  = {}
        self.legend_clades = []
        if colorize_clades is not None:
            cmap = plt.get_cmap('Set3', len(colorize_clades))
            for i, clade in enumerate(colorize_clades):
                self.legend_clades.append(patches.Patch(color=cmap(i, alpha=0.7),
                                                         label=clade))
                self.colorize_species.update({sp: cmap(i, alpha=0.7) for sp in
                                              self.phyltree.species[clade]})
            self.clade_cmap = cmap

        self.ancgene2sp = re.compile('(' + 'root|'
                        + '|'.join(re.escape(s) for s in
                                    list(self.phyltree.listSpecies) +
                                    sorted(self.phyltree.listAncestr,
                                           key=lambda a: len(a),
                                           reverse=True)).replace(r' ', r'.')
                        + ')(.*)$')
        
        self.taxa = set(self.phyltree.allNames)
        # Cleared and refilled if a genetree is given
        self.max_age = self.phyltree.ages[self.phyltree.root]

    def load_reconciled_genetree(self, filename, format=1, genetreename=None):
        """Load gene tree with all species nodes present.
        newick format with internal node labelling"""

        if genetreename:
            self.genetreename = genetreename
        elif filename.startswith('/dev/fd/'):
            self.genetreename = filename
        elif filename == '-':
            filename = '/dev/stdin'
            self.genetreename = 'stdin'
        else:
            self.genetreename = op.splitext(op.basename(filename))[0]

        #self.genetree = ete3.Tree(filename, format=format)
        with open(filename) as gt_input:
            # Allow reading multiple trees from a single input file.
            genetree_texts = [gt for gt in gt_input.read().split(';')
                              if gt.rstrip()]
            try:
                genetrees = [ete3.Tree(gt_txt+';', format=format)
                             for gt_txt in genetree_texts]
            except ete3.parser.newick.NewickError as err:
                err.args += tuple(genetree_texts)
                raise

        assert genetrees, "Input file contains no genetree."

        self.genetrees = []
        roots = set()
        alldescendants = set()
        self.taxa = set()

        for genetree in genetrees:
            ### TODO: actually, ladderizing is useless, it should sort speciation
            ###       nodes according to the average y position of next spe node.

            root = self.get_taxon(genetree, self.ancgene2sp, self.ensembl_version)
            alldescendants |= self.phyltree.allDescendants[root]

            # Add the branch leading to the current root, if **duplications** in
            # this branch, so that the new root is a speciation.
            lower_root_node = self.phyltree.parent.get(root)
            if lower_root_node:
                lower_root = lower_root_node.name

                # This check is not especially necessary.
                while len(self.phyltree.items.get(lower_root, [])) == 1 \
                        and self.phyltree.ages[lower_root] == 0:
                    lower_root = self.phyltree.parent[lower_root].name
            else:
                lower_root = "root"

            roots.add(lower_root)
            #alldescendants.add(lower_root)
            self.taxa.add(lower_root)

            #rerooted_genetree = ete3.TreeNode(name=lower_root) # Might break
            lower_root = root  # Put the root in the middle of the branch.
            rerooted_genetree = ete3.TreeNode(name=genetree.name.replace(root,
                                                                    lower_root))

            # Only needed for the TreeBest format
            rerooted_genetree.add_feature('S', lower_root)
            # See the meaning of the 'P' ("paralogy") and 'A' ("asymmetrical")
            # in dendro/reconciledtree2paralogytree.py
            rerooted_genetree.add_feature('P', getattr(genetree, 'P', 0))
            rerooted_genetree.add_child(child=genetree)
            self.genetrees.append(rerooted_genetree)

        #for genetree in genetrees:
        #    rerooted_genetree.add_child(child=genetree)
        #self.genetree = rerooted_genetree
        #self.genetrees = [reroot_genetree...]

        # add only the meaningful taxa (not those with one child and age = 0)
        for taxon in alldescendants:
            if not (len(self.phyltree.items.get(taxon, [])) == 1 and
                    self.phyltree.ages[taxon] == 0):
                self.taxa.add(taxon)

        #if len(roots) == 1:
        #    root, = roots
        #else:
        
        #print("Roots: " + ", ".join(roots))
        if self.ages:
            lca = self.phyltree.lastCommonAncestor(list(roots.difference(("root",))))
            self.max_age = self.phyltree.ages[lca]
            print('roots: %s; max_age: %g' % (roots, self.max_age))

        root = "root" if "root" in roots else self.phyltree.lastCommonAncestor(list(roots))
        self.taxa.add(root)
        

    def draw_species_tree(self, branch_width=0.8):
        """Init figure + draw branches of the species tree.
        
        branch_width: proportion of vertical space between two branches taken
                      by the branch polygon.
        angle_style: See the parser help.
        rootwardgrowth: widening factor of the branches towards the root, per x unit.
        internal: restricted list/set of internal nodes to print, or "all".
        """
        #TODO: implement the colorize_clades as branch color
        self.species_coords   = {}
        self.species_branches = {}  # Back link to the parent, and branch length.
        rootwardgrowth = 3./self.max_age if self.ages else 0.1 
        #rootwardgrowth = 0

        self.fig, ax0 = plt.subplots() #frameon=False) # set figsize later
        #ax0 = self.fig.add_axes([0.1,0.1,0.9,0.9]) #, adjustable='box-forced')
        #ax0 = self.fig.add_axes([0.1, 0.1, 0.9, 0.9])
        if not self.debug: ax0.axis('off')
        # TODO: if ages: keep x axis with xlabel "age"

        # I need a width from data coordinates to linewidth units:
        ymin = 0
        any_show_cov = False
        for parent, (px, py, _), children_coords in \
              iter_species_coords(self.phyltree, self.taxa, self.angle_style, self.ages):
            py += branch_width*abs(px)*rootwardgrowth*0.5
            pwidth = branch_width*(1+rootwardgrowth*abs(px))
            # Sort branches by their angle (branches going down first)
            sorted_children_coords = sorted(children_coords, key=lambda item: (item[1][1]-py) / (item[1][0]-px))
            #sorted_children_coords = list(reversed(children_coords))
            for i, (child, (cx, cy, _)) in enumerate(sorted_children_coords):
                cy += branch_width*abs(cx)*rootwardgrowth*0.5
                cwidth = branch_width*(1+rootwardgrowth*abs(cx))
                logger.debug('%-16s %-16s (%g, %g) -> (%g, %g); widths: %g -> %g', parent, child,
                             px, py, cx, cy, pwidth, cwidth)
                self.species_branches[child] = (parent, cx - px, cy - py, cwidth)

                self.species_coords[child] = (cx, cy, cwidth)
                if cy < ymin:  # update ymin.
                    ymin = cy
                coords = np.array([(px, py),
                                   (cx, cy),
                                   (cx, cy - cwidth),
                                   (px, py - pwidth)])
                    
                if i>0:
                    # If this not the first child (i.e. drawn over the others), update the polygon coords,
                    # so that shadows fall correctly at the fork.
                    # Compute the intersection of the bottom line going to child i, 
                    # and the top line going to child i-1
                    orig_yi = -pwidth  # Origin is at py
                    slope_i = (cy-cwidth - (py-pwidth)) / (cx-px)

                    child2, (cx2, cy2, _) = sorted_children_coords[i-1]
                    cy2 += branch_width*abs(cx2)*rootwardgrowth*0.5

                    #ax0.plot([px, cx2], [py, cy2], '-.', [px, cx], [py-pwidth, cy-cwidth], '-')
                    orig_y2 = 0
                    slope_2 = (cy2 - py)/(cx2-px)

                    #assert slope_2 < slope_i
                    x_inter = (orig_yi - orig_y2) / (slope_2 - slope_i) 

                    y_inter = slope_i*x_inter + orig_yi
                    #y_inter2 = slope_2*x_inter + orig_y2
                    #assert np.isclose(y_inter, y_inter2), '%g != %g' %(y_inter, y_inter2)

                    logger.debug('Child list: [%s]', ' '.join('%s(%g)'%(ch,achy) for ch,(_,achy,_), in sorted_children_coords))
                    logger.debug('cy2=%g VS species_coords[%s] = ... %g ...', cy2, child2, self.species_coords[child2][1])
                    logger.debug('intersect %s: %g + %g*x with %s: %g + %g*x at x=%g',
                                 child, orig_yi, slope_i, child2, orig_y2, slope_2, x_inter)
                    #assert x_inter > 0
                    #assert x_inter < cx - px
                    
                    #assert py-pwidth < y_inter < cy-cwidth, "%g %g %g" % (py-pwidth, y_inter, cy-cwidth)
                    #assert py > y_inter > cy2, "%g %g %g" % (py, y_inter, cy2)
                    coords[3,:] = [px + x_inter, py + y_inter]
                    #ax0.plot([px, px+x_inter, cx2], [py, py+line_2_start_y+line_2_slope*x_inter, cy2], '-')
                    #ax0.plot([px, px+x_inter, cx], [py-pwidth, py+line_i_start_y+line_i_slope*x_inter, cy-cwidth], '-')
                #else:
                #    logger.debug('No shade update needed at %s', child)
                # Shadow
                if coords[2,0] > coords[3,0]:
                    ax0.fill_between(coords[2:,0], coords[2:,1], coords[2:,1]-0.03, color='k', alpha=0.2,
                             zorder=-1)
                ax0.add_patch(patches.Polygon(coords, #facecolor='#e5e5e5', edgecolor='none')) #color='#808080'))
                                              facecolor='#80808033', edgecolor='none'))
                                              #edgecolor='#80808033'))
                              #path_effects=[patheffects.withSimplePatchShadow((1,-1))]))
                              #path_effects=[patheffects.Normal(),
                              #      patheffects.Normal(),
                              #      patheffects.Normal(),
                              #      patheffects.SimpleLineShadow((0.1, 0.1))]))
                              #path_effects=[patheffects.Normal(),
                              #      patheffects.SimplePatchShadow((0.1, 0.1)),
                              #      patheffects.Normal(),
                              #      patheffects.SimplePatchShadow((0.1, 0.1))]))

                # data-specific coloring (low-coverage)
                alpha = 1
                if self.show_cov and child in getattr(self.phyltree, "lstEsp6X", ()):
                    alpha = 0.6
                    any_show_cov = True
                elif self.show_cov and child in getattr(self.phyltree, "lstEsp2X", ()):
                    alpha = 0.3
                    any_show_cov = True

                bgcolor = self.colorize_species.get(child, '#ffffff00')

                text_xy = (cx, cy)
                ha = 'right'
                va = 'top'
                if cx >= 0:
                    # It is a species, add the common name
                    if not self.latinname:
                        common_name = get_common_name(self.phyltree, child)
                        if self.commonname:
                            child = common_name
                        elif common_name:
                            child += ', %s' % common_name
                    ha = 'left'
                elif cy-py >= 0:  # Branch going up
                    va = 'bottom'
                    alpha = 0.75
                else:  # Branch going down
                    text_xy = (cx, cy-cwidth)
                    alpha = 0.75

                if 'all' in self.internal or child in self.internal or cx >= 0:
                    ax0.annotate(child, text_xy, (2, 0), textcoords='offset points',
                                 ha=ha, va=va, fontsize=basefontsize,
                                 fontstyle='italic', family='serif', alpha=alpha,
                                 bbox={'facecolor': bgcolor, 'pad': 0, 'edgecolor': 'none'})
                    #logger.debug()

        # include root.
        self.species_coords[parent] = (px, py, pwidth)
        ax0.text(px, py, parent, ha='right', fontsize=basefontsize,
                 fontstyle='italic', family='serif')
        
        # Add legend in case of coverage information
        if self.show_cov and any_show_cov:
            legend_cov = ax0.legend(handles=self.legend_coverage,
                                    loc="upper left",
                                    bbox_to_anchor=(0, 1),
                                    title="Sequencing coverage",
                                    prop={'size': basefontsize, 'family': 'serif',
                                          'style': 'italic'},
                                    handlelength=0,
                                    facecolor='inherit')
            for label, line in zip(legend_cov.get_texts(), legend_cov.get_lines()):
                label.set_alpha(line.get_alpha())
                #label.set_style('italic')
                #label.set_family('serif')

        # Add legend in case of colorized_clades
        if self.legend_clades:
            ax0.legend(handles=self.legend_clades,
                       loc="upper left",
                       bbox_to_anchor=(0, 0.9),
                       title="Clades",
                       prop={'size': basefontsize, 'family': 'serif', 'style': 'italic'},
                       facecolor='inherit')
            if self.show_cov and any_show_cov: ax0.add_artist(legend_cov)

        if self.debug:
            ax0.autoscale_view()
        else:
            ax0.set_xlim(px, 0)
            ax0.set_ylim(ymin - 1, 1)
        self.ax0 = ax0
        return self.fig


    def minimize_gene_branch_crossing(self):
        """Reduce the number of gene branches crossing at speciation nodes, 
        by rotating duplication forks."""

        _, phylsubtree = self.phyltree.getSubTree(self.taxa)

        # Store for each node the number of species deletion going down or up.
        deletion_count = {} # +1 for deletion in high taxon, -1 in low taxon
        get_taxon_y = lambda taxon: self.species_coords[taxon][1]
        cached_taxa = {}

        def orient_deletion(children_taxa, expected_children_taxa):
            """Tell if the deleted gene branch belong to a species branch going
            up (+1) or down (-1)"""
            assert expected_children_taxa, children_taxa
            expected_children_taxa_ys = sorted(expected_children_taxa, key=get_taxon_y)
            low_taxon = expected_children_taxa_ys[0]
            high_taxon = expected_children_taxa_ys[-1]

            if low_taxon not in children_taxa:
                return -1
            elif high_taxon not in children_taxa:
                return 1
            else:
                #raise ValueError(children_taxa, expected_children_taxa,
                #                 expected_children_taxa_ys)
                return 0

        def iter_missing_branches(child_taxon, expected_children_taxa):
            """When intermediate speciation nodes are missing between a taxon and
            descendant taxa, list the implicit branches (a pair of a taxon and its child
            taxon)."""
            logger.info("  Missing between %s and %s: %s",
                        child_taxon, expected_children_taxa,
                        child_taxon not in expected_children_taxa)
            while child_taxon not in expected_children_taxa:
                tmp_taxon = self.phyltree.parent[child_taxon].name
                logger.info("  - get intermediate parent of %s: %s",
                            child_taxon, tmp_taxon)
                while tmp_taxon not in phylsubtree:
                    try:
                        tmp_taxon = self.phyltree.parent[tmp_taxon].name
                    except KeyError as err:
                        err.args += ('from child %s' % child_taxon,)
                        logger.error('%s %s', err, phylsubtree)
                        raise

                yield tmp_taxon, child_taxon
                
                child_taxon = tmp_taxon

        for node in (n for genetree in self.genetrees \
                                      for n in genetree.traverse('postorder')):
            taxon = self.get_taxon(node, self.ancgene2sp, self.ensembl_version)
            cached_taxa[node] = taxon
            children_taxa = set(cached_taxa[ch] for ch in node.children)
            event = infer_gene_event(node, taxon, children_taxa)

            deletion_count[node] = sum(deletion_count[ch] for ch in node.children)

            if event == 'spe':

                expected_children_taxa = set(phylsubtree[taxon])
                
                assert taxon not in children_taxa, \
                    "This should be a duplication: %s -> %s" % (taxon, children_taxa)
                if expected_children_taxa != children_taxa:
                    logger.warning('Expected (spe): %s ; got: %s (from %s)',
                                   expected_children_taxa, children_taxa, taxon)

                    tmp_children_taxa = []
                
                    # if intermediate speciation nodes are missing, count the
                    # implied deletions.
                    for child_taxon in children_taxa:

                        # It works when child_taxon is **descendant** of expected_taxa.
                        for tmp_taxon, child_taxon in iter_missing_branches(child_taxon,
                                                          expected_children_taxa):

                            tmp_expected_children_taxa = phylsubtree[tmp_taxon]

                            deletion_count[node] += orient_deletion(set((tmp_taxon,)), tmp_expected_children_taxa)
                            logger.info('+1 deletion')
                        tmp_children_taxa.append(child_taxon)

                    deletion_count[node] += orient_deletion(children_taxa, expected_children_taxa)

            elif event == 'dup':
                if len(children_taxa) > 1:
                    for child in node.children:
                        child_taxon = cached_taxa[child]
                        #is_transfer = (getattr(node, 'T', None) is not None
                        #            and getattr(child, 'T', None) is not None)
                        is_transfer = ((getattr(node, 'T', None) in ('1', '0'))
                                    and (getattr(child, 'T', None) in ('-1', '0')))
                        logger.info('node.T = %r; child.T = %r => is_transfer = %s',
                                    getattr(node, 'T', None),
                                    getattr(child, 'T', None), is_transfer)
                        if child_taxon != taxon and not is_transfer:
                            logger.warning('Expected (dup): %s ; '
                                           'got: %s (from %s)',
                                           taxon, children_taxa, taxon)
                            for tmp_taxon, _ in iter_missing_branches(child_taxon,
                                                        set((taxon,))):
                                tmp_expected_children_taxa = phylsubtree[tmp_taxon]
                                deletion_count[child] += orient_deletion(set((tmp_taxon,)), tmp_expected_children_taxa)

                node.children = sorted(node.children, key=deletion_count.get)

        self.genetrees = sorted(self.genetrees, key=deletion_count.get)


    def set_gene_coords(self, asymmetric=False, colorize_descent=None):
        """- asymmetric: whether to draw *asymmetric* divisions in the gene
                        tree. If True, the tree nodes descending from a
                        duplication must contain a "A=1" tag when they are
                        children from the main branch.
        """
        ### TODO: enable repeted calls to set multiple gene trees
        ### TODO: propagate deleted lineage : draw invisible lines, just to
        ###       consume vertical space.
        
        self.gene_coords = {} # {species: [gene list]}
        #self.dup_branchings = [] # (nodeid/genename, x, x_child1, x_child2, y_child1, y_child2)
        #self.spe_branchings = []
        self.branchings = {} # {nodeid: [species, x, y, dup/spe,
                             #           (children), {node features}], ...]
        
        # *temporary* tree structure holding nodes that are inbetween 2 speciations
        interspecies_trees = {}

        # x and y are relative to the branch considered. (between 0 and 1)

        self.minimize_gene_branch_crossing()

        descent_i = 0
        # Node names triggering special color for all descendant lines.
        self.colorize_descent = []  # Sorted in the order they were found.
        if not colorize_descent: colorize_descent = set()

        for nodeid, node in enumerate(n for genetree in self.genetrees \
                                      for n in genetree.traverse('postorder')):
            # set up a *unique* node ID, to avoid relying on node.name.
            node.id = nodeid
            #if node.is_root(): print('At Root:', node.name,node.id)

            taxon = self.get_taxon(node, self.ancgene2sp, self.ensembl_version)
            taxon_gene_coords = self.gene_coords.setdefault(taxon, [])

            children_taxa = set(interspecies_trees[ch.id]['taxon'] for ch
                                in node.children)
            event = infer_gene_event(node, taxon, children_taxa)
            #FIXME: take a decision on whether to copy those features or just use the ete3.TreeNode object.
            node_features = {ft: getattr(node, ft) for ft in node.features}
            if node.name in colorize_descent:
                node_features['C'] = descent_i  # colorize the *leading* branch
                for descendant in node.iter_descendants():
                    try:
                        # Not working if that is a dup.
                        self.branchings[descendant.id][5]['C'] = descent_i
                    except KeyError:
                        #print('%s (id=%d) ... %s leaf=%s' % (node.name, nodeid, descendant.name, descendant.is_leaf()), file=sys.stderr)
                        #print('branchings:', sorted(self.branchings), 'interspecies_trees', sorted(interspecies_trees), file=sys.stderr)
                        interspecies_trees[descendant.id]['features']['C'] = descent_i
                descent_i += 1
                self.colorize_descent.append(node.name)
            elif colorize_descent:
                node_features['C'] = -1  # falls under the colormap range -> black

            if event == 'leaf':
                taxon_gene_coords.append(nodeid)
                node_y = len(taxon_gene_coords) - 1 # This list is being
                                                    # extended later on.
                interspecies_trees[nodeid] = {'taxon': taxon,
                                              'ndup': 0,
                                              'x': 1,
                                              'y': node_y,
                                              'asymmetry': int(getattr(node, 'A', 0)),
                                              'features': node_features}
                self.branchings[nodeid] = [taxon, 0, node_y, 'leaf', [], node_features]
            elif event == 'dup':

                # It is a **duplication**
                # At this stage, branchings[nodeid] is unset, it will only be set
                # on the second pass (preorder from the ancestral speciation).
                children_ys   = []
                children_ndup = []
                for ch in node.children:
                    children_ys.append(interspecies_trees[ch.id]['y'])
                    children_ndup.append(interspecies_trees[ch.id]['ndup'])
                if not colorize_descent and not node.is_root():
                    # We added a fake root node in load_reconciled_genetree but it's recognized as a dup.
                    # Colorize by duplication fork, as initially implemented.
                    for ch in node.children:
                        interspecies_trees[ch.id]['features']['C'] = descent_i
                        logger.debug('Child of dup: set color #%d at %s' % (descent_i, ch.name))
                        # NOTE: very implicit behavior: if child is a speciation, 
                        # its branchings dict will not be updated. However the 
                        # node_features in interspecies_trees and branchings is
                        # the same object, so this will be updated in both places.
                    descent_i += 1
                        
                # TODO: tilt according to the parent speciation node.
                # should fit in the triangle defined by parent and descendants.
                #node_y    = (min(children_ys) + max(children_ys)) / 2
                
                node_ndup = max(children_ndup) + 1
                if not asymmetric:
                    # weighted average of the children positions.
                    node_y    = sum(cy*(cd+1) for cy,cd in \
                                               zip(children_ys, children_ndup))
                    node_y   /= (sum(children_ndup) + len(children_ndup))
                else:
                    main_branches = [ch for ch in node.children
                                     if not interspecies_trees[ch.id]['asymmetry'] # int(['features'].get('A', 0))
                                    ] or node.children
                    node_y = interspecies_trees[main_branches[0].id]['y']

                interspecies_trees[nodeid] = {
                        'taxon': taxon,
                        'ndup': node_ndup,
                        'x': 1, # modified later
                        'y': node_y,
                        'asymmetry': int(getattr(node, 'A', 0)),
                        'features': node_features}
            if event == 'spe' or node.is_root():
            #else:
                # It is a **speciation**
                if event == 'spe':
                    taxon_gene_coords.append(nodeid)
                    node_y = len(taxon_gene_coords) - 1
                    node_ndup = 0
                    node_x = 0

                    interspecies_trees[nodeid] = {
                            'taxon': taxon,
                            'ndup': node_ndup,#0, 
                            'x': node_x,
                            'y': node_y,
                            'asymmetry': int(getattr(node, 'A', 0)),
                            'features': node_features}
                else:
                    # If it is the root, place it in the middle of the branch.
                    node_x = 1 - node_ndup/(node_ndup+1)
                    interspecies_trees[nodeid]['x'] = node_x
                    #interspecies_trees[nodeid]['features']['C'] = -1  # Black root edge, always.

                self.branchings[nodeid] = [taxon, node_x, node_y, event,
                                           [ch.id for ch in node.children],
                                           node_features]
                #if event == 'dup':
                #    continue

                #print('Speciation %r' % taxon, self.branchings[nodeid])
                
                for ch in node.children:
                    n_steps = interspecies_trees[ch.id]['ndup'] + 1
                    if node.is_root(): n_steps += 1

                    delta_x = 1 / n_steps

                    # climb up this subtree until next speciation
                    # (toward present)
                    nextnodes = [ch]
                    while nextnodes:
                        nextnode = nextnodes.pop(0)
                        try:
                            nextnode_ndup = interspecies_trees[nextnode.id]['ndup']
                        except KeyError as err:
                            raise KeyError('Non unique node name: %r', nextnode.id)

                        if nextnode_ndup > 0:
                            nextnode_x = 1 - nextnode_ndup * delta_x
                            nextnode_y = interspecies_trees[nextnode.id]['y']
                            nextnode_taxon = interspecies_trees[nextnode.id]['taxon']
                            nextnode_ch = [nnch.id for nnch in nextnode.children]
                            nextnode_features = interspecies_trees[nextnode.id]['features']

                            self.branchings[nextnode.id] = [
                                                            nextnode_taxon,
                                                            nextnode_x,
                                                            nextnode_y,
                                                            'dup',
                                                            nextnode_ch,
                                                            nextnode_features]
                            nextnodes.extend(nextnode.children)

                        interspecies_trees.pop(nextnode.id)

        #assert not interspecies_trees, \
        #        "Temporarily stored nodes have not been processed: %s" % interspecies_trees
        print("Root:", interspecies_trees)


    def draw_gene_tree(self, extratitle='', genenames=False, tags="",
                       fork_style='curved'):
        self.drawn_count += 1
        print(' --- Drawing genetree %d ---' % self.drawn_count)
        # Duplicate the species tree axis to separate the plotting
        self.ax1 = self.ax0.twinx()
        self.ax1.set_ylim(self.ax0.get_ylim())
        if not self.debug:
            self.ax1.axis('off')
        title = self.genetreename
        if extratitle: title += ' -- %s' % extratitle
        self.ax1.set_title(title)
        
        self.real_gene_coords = {}

        # Colors for edges
        #cmap = plt.get_cmap('Dark2', 10)
        #colorcycle = mpl.rcParams['axes.prop_cycle'].by_key()['color']
        # Based on Seaborn 'bright' palette, without the dark grey color, extended with some colors from 'colorblind'
        # dark:['#001c7f', '#b1400d', '#12711c', '#8c0800', '#591e71', '#592f0d', '#a23582', '#006374', '#b8850a', # too dark
        # muted:['#4878d0', '#ee854a', '#6acc64', '#d65f5f', '#956cb4', '#8c613c', '#dc7ec0', '#82c6e2', '#d5bb67', # Ok.
        # default:['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#17becf', '#bcbd22',
        cmap = mpl.colors.ListedColormap(
            ['#023eff', '#ff7c00', '#1ac938', '#e8000b', '#8b2be2', '#9f4800', '#f14cc1', '#00d7ff', '#ffc400',
             '#fbafe4', '#56b4e9', '#ece133', '#029e73'], name='custom')
        cmap.set_under('k')  # Better grey0.1 ?
        #color_index = 0

        # Gene node names to be displayed (dup, leaf, spe)
        genenames = "leaf,dup,spe" if genenames in (True,"all") else genenames
        genenames = set(genenames.split(',')) if genenames else set()
        # node tags to be displayed
        tags = set(tags.split(",")) if tags else set()

        if fork_style == "curved":
            fork_instructions = [MOVETO, CURVE3, CURVE3]
        elif fork_style == "square":
            fork_instructions = [MOVETO, LINETO, LINETO]
        else:
            raise ValueError("fork_style = %r (Should be in %s)" % \
                             (fork_style, ("curved", "squared")))
        
        for node in (n for genetree in self.genetrees \
                     for n in genetree.traverse('postorder')):
            nodeid = node.id
            
            species, rel_x, rel_y, event, children, features = self.branchings[nodeid]
            pos_list = self.gene_coords[species]
            nranks = len(pos_list) + 1
            children_real_coords = [self.real_gene_coords[ch] for ch in children]
            children_features = [self.branchings[ch][5] for ch in children]
            
            #print(nodeid, event, species, children)
            if event == 'dup':

                # coordinates of the **species** branch vector
                # Dx = x (taxon) - x (ancestor) >0
                parent_sp, Dx, Dy, width = self.species_branches[species]
                parent_x, parent_y, parent_width = self.species_coords[parent_sp]
                # When branch width is constant, it's something like that:
                #branch_transform = Affine2D().translate(parent_x, parent_y-width).scale(abs(Dx), abs(Dy)).skew(0, np.arctan2(Dy, abs(Dx)))
                #real_x, real_y = branch_transform.transform_point((rel_x, rel_y))
                # Branch width at the x position of the node
                branch_width = parent_width * (1-rel_x) + width * rel_x
                real_x = parent_x + Dx * rel_x
                real_y = parent_y + Dy * rel_x - branch_width * (rel_y+1)/nranks
                
                nodecolor = 'black' if node.is_root() else 'red'

                nch = len(children)
                children_rel_ys = [self.branchings[ch][2] for ch in children]

                assert all((ch_rel_y >= 0) for ch_rel_y in children_rel_ys)
                if any((children_rel_ys[i+1] - children_rel_ys[i] < 0) for i in range(nch-1)): 
                    logger.error("Children's relative Y not sorted! (%s: %s)",
                                  event, node.name)
                elif any((children_rel_ys[i+1] - children_rel_ys[i] == 0) for i in range(nch-1)):
                    logger.error("Some children's relative Y are identical! (%s: %s)",
                                  event, node.name)

                for ch, ch_rel_y, (ch_real_x, ch_real_y), ch_ft in \
                        zip(children, children_rel_ys, children_real_coords,
                            children_features):

                    # Draw thicker line when it represents a paralogy
                    linewidth = 4 if int(ch_ft.get('P', 0)) else 1
                    
                    delta_y = (rel_y - ch_rel_y)/nranks
                    fork_coords = [(ch_real_x, ch_real_y),
                                   (real_x, real_y + delta_y),
                                   (real_x, real_y)]
                    color_cat = ch_ft.get('C', -1)
                    branch_color = cmap(color_cat) if color_cat<0 else cmap(color_cat % cmap.N)
                    u_fork_finger = patches.PathPatch(
                                    Path(fork_coords, fork_instructions),
                                    fill=False,
                                    edgecolor=('k' if int(ch_ft.get('A', -1))==0
                                               else branch_color),
                                    alpha=0.5,
                                    linewidth=linewidth,
                                    joinstyle='round')  # don't see a change
                    self.ax1.add_patch(u_fork_finger)

            else:  # event == 'spe' or 'leaf'
                real_x, real_y, branch_width = self.species_coords[species]
                pos = pos_list.index(nodeid) + 1
                real_y -= branch_width * pos/nranks
                nodecolor = 'blue'
                #nodecolor = 'none'

                for ch, (ch_real_x, ch_real_y), ch_ft in \
                        zip(children, children_real_coords, children_features):
                    linewidth = 4 if int(ch_ft.get('P', 0)) else 1
                    color_cat = ch_ft.get('C', -1)
                    v_fork_finger = patches.PathPatch(
                                             Path([(real_x, real_y),
                                                   (ch_real_x, ch_real_y)],
                                                  [MOVETO, LINETO]),
                                             edgecolor=(cmap(color_cat) if color_cat<0 else cmap(color_cat % cmap.N)),
                                             alpha=0.5,
                                             linewidth=linewidth,
                                             joinstyle='round')
                    self.ax1.add_patch(v_fork_finger)
                                  #(':' if node.is_root() else '-'),
            
            #if event != 'leaf':
            if event == 'dup' and not node.is_root():  # This line is still here because historically it was.
                # Add a *dot*
                self.ax1.plot((real_x,), (real_y,), '.',
                              #('<' if node.is_root() else '.'),
                              color=nodecolor, alpha=0.5)
                              #, markeredgewidth=linewidth)
                #self.ax1.text(real_x, real_y, species, fontsize='xxx-small',
                #              color=nodecolor)
            if event in genenames:
                node_annot = [] if '-name' in tags else [node.name]
                for ft in tags.difference(('-name',)):
                    try:
                        node_annot.append("%s=%s" % (ft, features[ft]))
                    except KeyError:
                        logger.debug("Attribute %r not found in node %s",
                                     ft, ' '.join(node_annot))
                        pass

                self.ax1.text(real_x, real_y, ' '.join(node_annot), alpha=0.5,
                              fontsize=fontsizes[max(0, fontsizes.index(basefontsize)-1)],
                              ha=('left' if event=='leaf' else 'right'),
                              va='top')

            self.real_gene_coords[nodeid] = (real_x, real_y)
            #self.fig.draw(self.fig.canvas.get_renderer())
            #plt.show()
            #input('Press Enter to continue.')
            ### TODO: add onpick action: display node name
        if self.colorize_descent:
            self.ax1.legend(handles=[lines.Line2D([],  [], color=cmap(i % cmap.N))
                            for i,_ in enumerate(self.colorize_descent)],
                       labels=self.colorize_descent,
                       loc="lower left",
                       title="Descending from:",
                       prop={'size': basefontsize, 'family': 'serif', 'style': 'italic'},
                       facecolor='#e5e5e5', framealpha=1)  # Same as species branch.

    def draw(self, genetree, extratitle='',
             genenames=False, tags="", asymmetric=False,
             colorize_descent=None):
             #fork_style="curved"):
        """Once phyltree is loaded, perform all drawing steps."""
        self.load_reconciled_genetree(genetree)
        self.draw_species_tree()
        self.set_gene_coords(asymmetric, colorize_descent)
        self.draw_gene_tree(extratitle, genenames=genenames, tags=tags,
                            fork_style=("square" if asymmetric else "curved"))
        self.fig.tight_layout()
        return self.fig


###
### Part of the module that tries to automatically prepare the data for drawing,
### from a single ensembl genetree name.
###

def check_extracted(genetrees, output):
    all_outputs = [output.format(genetree=gt) for gt in genetrees]
    new_genetrees = [gt for gt in genetrees if not op.exists(output.format(genetree=gt))]
    new_outputs = [output.format(genetree=gt) for gt in new_genetrees]
    print("The following genetrees are already extracted: %s" %
            (set(all_outputs) - set(new_outputs),), file=sys.stderr)
    return all_outputs, new_outputs, new_genetrees


def prepare(genetrees, ancestors, ensembl_version=ENSEMBL_VERSION, ncores=1,
            edited=True, subtrees_dir='subtrees_', dry_run=False):
    """
    Prepare genetree files (if needed) starting from raw Genomicus/Ensembl data:
      1. find them in the tree forest;
      2. reconcile them with the species tree;
      3. extract the subtrees starting at a given ancestor.
    """
    # 1. Find ancgene name from modern gene?
    # 2. Find the orthologs/ancestors in the given ancestors

    datadir = '/users/ldog/glouvel/ws2/DUPLI_data%d/alignments' % ensembl_version
    assert op.exists(datadir)

    if edited:
        # Take gene tree from Genomicus
        treeforestfile = "/users/ldog/glouvel/ws_alouis/GENOMICUS_SVN/data%d/"\
                         "GoodThreshold/tree.4F.cut.bz2" % ensembl_version
        withAncGenesNames = True
        field = 'family_name'
        output = op.join(datadir, '{genetree}', '{genetree}.nwk')
        fix_suffix = True
    else:
        # Take gene tree from Ensembl
        treeforestfile = "/users/ldog/glouvel/ws_alouis/GENOMICUS_SVN/data%d/"\
                         "tree.1.ensembl.bz2" % ensembl_version
        withAncGenesNames = False
        field = 'tree_name'
        output = op.join(datadir, '{genetree}', '{genetree}_ensembl.nwk')
        fix_suffix = False

    all_outputs, new_outputs, new_genetrees = check_extracted(genetrees, output)
    
    # Extract genetrees
    if new_outputs:
        import ToolsDyogen.treeTools.ALL.extractMultipleGeneTree as xMulti
        xMulti.main(treeforestfile, new_genetrees, field=field, toNewick=True,
                    withAncSpeciesNames=True,
                    withAncGenesNames=withAncGenesNames,
                    output=output, mkdirs=True)
    else:
        print("No new genetrees to extract.", file=sys.stderr)

    # Create output dirs for prune2family
    gt_format = '{0}' if edited else '{0:.%d}' % len(genetrees[0])
    #print(gt_format)
    
    prune_outdir = op.join(datadir, gt_format, subtrees_dir)

    for gt in genetrees:
        p_outdir = prune_outdir.format(gt)
        print(p_outdir, end=' ')
        if not op.exists(p_outdir):
            print("make!")
            os.mkdir(p_outdir)
        else:
            print("exists.")

    # prune genetrees to the right ancestor
    import codeml.prune2family as prune2family

    return prune2family.parallel_save_subtrees(all_outputs, ancestors,
                                           ncores=ncores,
                                           outdir=prune_outdir,
                                           only_dup=False, one_leaf=True,
                                           fix_suffix=fix_suffix,
                                           dry_run=dry_run, ignore_errors=False,
                                           ensembl_version=ensembl_version)



#TESTTREE = "/users/ldog/glouvel/ws2/DUPLI_data85/alignments/ENSGT00810000125388/subtrees2/RodentiaENSGT00810000125388.A.a.a.c.a.b.nwk"
#TESTTREE = "/users/ldog/glouvel/ws2/DUPLI_data85/alignments/ENSGT00850000132243/subtrees2/SimiiformesENSGT00850000132243.b.q.b.a.a.a.b.nwk"
TESTTREE = "/users/ldog/glouvel/ws2/DUPLI_data85/alignments/ENSGT00850000132243/subtrees2/SimiiformesENSGT00850000132243.b.q.b.b.a.b.b.a.b.c.a.a.a.nwk"


def run(genetrees, gene_params, outfile, genenames=False, tags="", asymmetric=False,
        colorize_descent=None, **kwargs):

    figsize = PAPERSIZE['a4']
    match_figsize = FIGSIZE.search(outfile)
    sizestr = 'a4'
    if match_figsize:
        sizestr = match_figsize.group().lower()
        outfile = outfile[:match_figsize.start()-1]
        #print(match_figsize.group())
        try:
            figsize = PAPERSIZE[sizestr]
        except KeyError:
            figsize = tuple(float(x) for x in match_figsize.groups())
    logger.debug('Will set output figsize to: %d x %d inches.', *figsize)
    gd = GenetreeDrawer(**kwargs)
    display = lambda: plt.show() # Display function for shell or notebook usage
    if __name__=='__main__' and outfile == '-':
        try:
            plt.switch_backend('Qt5Agg')
        except ImportError:
            try:
                plt.switch_backend('Qt4Agg')
            except ImportError:
                plt.switch_backend('TkAgg')
        #mpl.use('Qt4Agg')
        #from importlib import reload; reload(plt)
    elif outfile.endswith('.pdf'):
        pdf = PdfPages(outfile)
        display = lambda: (pdf.savefig(bbox_inches='tight',
                           papertype=(sizestr if sizestr in PAPERSIZE else 'a4')),
                           plt.close())
    else:
        assert len(genetrees) <= 1, "multipage output only supported for pdf"
        display = lambda: (plt.savefig(outfile, bbox_inches='tight'),
                           plt.close())

    for (genetree, gene_kwargs) in zip_longest(genetrees, gene_params, fillvalue={}):
        # Use global gene kwargs as default:
        gene_kwargs = {'genenames': genenames, 'tags': tags,
                       'asymmetric': asymmetric, 'colorize_descent': colorize_descent,
                       **gene_kwargs}
        try:
            genetree, *extratitles = genetree.split(',')
            extratitle = ', '.join(extratitles)
            print('INPUT FILE:', genetree, '(%s)' % extratitle)
            gd.draw(genetree, extratitle, **gene_kwargs)
            gd.fig.set_size_inches(figsize)
            display()
        except BaseException as err:
            if outfile.endswith('.pdf'):
                pdf.close()
            raise

    if not genetrees:
        gd.draw_species_tree().set_size_inches(figsize)
        display()

    if outfile.endswith('.pdf'):
        pdf.close()
    #plt.show()
    return gd


if __name__ == '__main__':

    logging.basicConfig(format='%(levelname)s:%(name)s:l.%(lineno)d:%(funcName)s:%(message)s')

    def commasplit(arg):
        return arg.split(',')

    parser = argparse.ArgumentParser(description=__doc__,
                        #formatter_class=CustomHelpFormatter)
                        formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('outfile',
                        help=("pdf file, or '-'. If '-', will use Qt to "
                              "display the figure. You can control the figsize"
                              " by suffixing with for example \":6x8\" (width "
                              "6, height 8 inches)."))
    parser.add_argument('genetrees', nargs='*', default=[],
        help=("must be a genetree (nwk format with internal nodes labelling) "
            "reconciled with species tree, or a genetree formatted like "
            "`TreeBest` output. '-' means standard input."))
    #parser.add_argument_group("Common arguments")
    parser.add_argument('--fromfile', action='store_true',
                        help='take genetree paths and description from a file')
    parser.add_argument('-p', '--phyltreefile', default=PHYLTREEFILE,
                        help='species tree in phyltree or newick format [%(default)s]')
    parser.add_argument('-e', '--ensembl-version', type=int,
                        default=ENSEMBL_VERSION)
    parser.add_argument('-t', '--treebest', action='store_true',
                        help='The input genetree is a treebest output')
    
    parser.add_argument('-a', '--angle-style', type=int, choices=range(6),
                        default=0,
                        help=(
                              "0: parent node positioned at x-1 and parent y is"
                              " closer to the child with furthest x.\n" # branch angles are equal
                              "1: parent node y: inversed weighted-average (less"
                              " crowded node gets more weight);\n"
                              "2: parent node y: mean of children nodes;\n"
                              "3: branches always at 45 degrees;\n"
                              "4: parent node y: weighted average of children y;\n"
                              "5: equal angles.\n[%(default)s]"
                              ))
    parser.add_argument('-C', '--commonname', action='store_true', 
                        help='Species common names only')
    parser.add_argument('-L', '--latinname', action='store_true', 
                        help='Species scientific names only')
    parser.add_argument('-i', '--internal', default='',
                        help='Internal node names to display or "all" (comma-separated).')
    
    parser.add_argument('-c', '--colorize-clade',
                        dest='colorize_clades', metavar='CLADES', default=None,
                        type=commasplit,
                        help='Set a specific label color to species in these '\
                             'clades (comma-separated)')
    parser.add_argument('-s', '--show-cov', action='store_true',
                        help='Show genome coverage information (grey shading)')
    parser.add_argument('-A', '--ages', action='store_true',
                        help='Place species nodes at their real age.')
    parser.add_argument('-d', '--debug', action='store_true',
                        help='More verbose output and graphical hints (axes tick values)')
    g_pars = parser.add_argument_group("Gene tree control",
                        "Can be superseeded by key=value pairs in the input file (--fromfile)")
    g_pars.add_argument('-g', '--genenames',
                        help='Display gene names: \n' \
                             '- comma-sep list of "leaf", "dup", "spe";\n'
                             '- "all" (identical to "leaf,dup,spe").')
    g_pars.add_argument('-T', '--tags', default="",
                        help="Additional node tags to write")
    g_pars.add_argument('-k', '--asymmetric', action='store_true',
                        help='Draw *asymmetric* gene duplications: one branch'\
                             ' stays the main branch, and other are children.')
    g_pars.add_argument('-l', '--colorize-descent', default=None, metavar='NODE NAMES',
                        type=commasplit,
                        help='Set a specific color to all gene lines descending'\
                             'from the given node (comma-separated)')
    gene_param_converter = {'colorize_descent': commasplit,
                            'asymmetric': lambda v: (True if v.lower()=='true' else
                                                      False if v.lower()=='false' else v)
                            }
                        
    #parser.add_argument('-m', '--multiple-pdfs', action='store_true',
    #                    help='output one pdf file per genetree. [NOT implemented]')
    args = parser.parse_args()
    dictargs = vars(args)
    #if not dictargs.get('genetrees'):
    #    dictargs['genetrees'] = [TESTTREE]
    genetrees = dictargs.pop('genetrees')
    gene_params = []
    if dictargs.pop('fromfile'):
        genetreelistfile = genetrees.pop()
        if genetrees:
            print('Argument error: only one input file allowed with --fromfile')
            sys.exit(2)
        with (sys.stdin if genetreelistfile=='-' else open(genetreelistfile)) as stream:
            for line in stream:
                if not line.startswith('#'):
                    genetree_str, *gene_kwarg_str = line.rstrip().split('\t')
                    genetrees.append(genetree_str)
                    gene_kwargs = {}
                    for keyval in gene_kwarg_str:
                        key, val = keyval.split('=')
                        try:
                            gene_kwargs[key] = gene_param_converter[key](val)
                        except KeyError:
                            gene_kwargs[key] = val
                    gene_params.append(gene_kwargs)
                    
            
    # TODO: add into run()
    #ANCGENE2SP = re.compile(r'([A-Z][A-Za-z_.-]+)%s' % ancgene_regex)

    logger.setLevel(logging.DEBUG if args.debug else logging.INFO)

    gd = run(genetrees, gene_params, **dictargs)

