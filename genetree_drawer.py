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
import os.path
import re
import warnings
try:
    import argparse_custom as argparse # Less verbose help message
except ImportError:
    import argparse

import numpy as np
import matplotlib as mpl
#mpl.use('Agg')
if __name__=='__main__': mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.lines as lines
from matplotlib.backends.backend_pdf import PdfPages, FigureCanvas
from matplotlib.path import Path
MOVETO, CURVE3, LINETO = Path.MOVETO, Path.CURVE3, Path.LINETO
import ete3

import LibsDyogen.myPhylTree as PhylTree
#import LibsDyogen.myProteinTree as ProteinTree

from glou_duphist import dfw_descendants_generalized, ladderize
from codeml.select_leaves_from_specieslist import convert_gene2species
from prot2gene import convert_prot2species
from seqtools.specify import load_conversion


ENSEMBL_VERSION = 85
PHYLTREEFILE = "/users/ldog/glouvel/GENOMICUS{0}/PhylTree.Ensembl.{0}.conf"
ANCGENE2SP = re.compile(r'([A-Z][A-Za-z_.-]+)ENS') # NOT USED

try:
    ucsc_conv_filename = "~glouvel/ws2/UCSC_genome_releases_full.tsv"
    UCSC_CONVERSION = load_conversion(ucsc_conv_filename)
except FileNotFoundError:
    warnings.warn("UCSC species name conversion file %r not found" % ucsc_conv_filename)
    UCSC_CONVERSION = {}

### Matplotlib graphical parameters ###
grey10 = '#1a1a1a'
mpl.rcParams['figure.figsize'] = [11.7, 8.27] # a4
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


def get_taxon(node, ancgene2sp, ensembl_version=ENSEMBL_VERSION):
    """from a gene name in my newick gene trees, find the taxon:
        either:
            - node is a leaf (e.g ENSMUSG00...)
            - node is internal (e.g Mus.musculusENSGT...)"""
    if node.is_leaf():
        try:
            taxon = convert_gene2species(node.name, ensembl_version)
        except RuntimeError:
            try:
                taxon = convert_prot2species(node.name, ensembl_version)
            except KeyError:
                taxon = UCSC_CONVERSION[re.match('[A-Za-z0-9]+', node.name).group()]
    else:
        try:
            taxon = ancgene2sp.match(node.name).group(1).replace('.', ' ')
        except AttributeError:
            raise RuntimeError("Can not match species name in %r" % node.name)
    return taxon


def get_taxon_treebest(node, *args):
    """get the taxon of a gene node, using a treebest output tree.
    
    *args are not used, they are here for compatibility with `get_taxon`"""
    try:
        return node.S.replace('.', ' ')
    except AttributeError:
        print(node.name)
        print(node)
        raise


def walk_phylsubtree(phyltree, taxa):
    """Return an iterator that will progress on the tree from leaves to root.
    
    It also rotates forks such that the displayed tree will show branches
    that are most distant from the root (in number of nodes) on one side.

    Taxa returned by the iterator are space-separated.
    """
    if "root" in taxa:
        root = "root"
        taxa.remove("root")
        oldroot, subtree = phyltree.getSubTree(taxa)
        subtree[root] = [oldroot]
    else:
        root, subtree = phyltree.getSubTree(taxa)

    # reorder branches in a visually nice manner:
    get_children = lambda tree, node: tree.get(node, [])
    ladderize(subtree, root, get_children, light_on_top=False)
    dfw = dfw_descendants_generalized(subtree, get_children,
                                      include_leaves=False, queue=[root])
    return reversed(list(dfw))


def iter_species_coords(phyltree, taxa, angle_style=0, ages=False):
    """Assign a pair x,y of coordinates for each node in the tree.
    Yield (parent name, parent_xy, child name, child_xy).
    
    - angle_style: 0 to 5, see parser help
    - ages:
        True: use real species ages
        'log': use log10(1 + species age)
        else: just the topology, real branch lengths are ignored.
    """
    # Store the graph coords (x,y) of a node, + its weight w (number of leaves)
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
                parent_x = phyltree.ages[parent]
            if ages == 'log':
                parent_x == np.log10(1 + parent_x)
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
            elif angle_style == 3:
                # branches are all at 45 degrees
                step = (dy - dx) / 2
                parent_x = min(children_xs) - step
                parent_y = max(children_ys) - step
            elif angle_style == 5:
                # equal angles
                step = dy / (2+dx)
                parent_y = max(children_ys) - step
            else:
                raise ValueError("Invalid angle_style value: %r" % angle_style)

        coords[parent] = (parent_x, parent_y, parent_w)

        for child in children:
            yield parent, coords[parent], child, coords[child]


def infer_gene_event(node, taxon, children_taxa):
    """Tell whether a gene tree node (ete3 format) is a leaf, a speciation or
    a duplication, using the taxon information.
    
    Children taxa must be a set (because len() is used)"""
    #if node.is_leaf():
    if not children_taxa:
        return 'leaf'
    elif (getattr(node, 'D', None) == 'Y') or taxon in children_taxa:
       #(len(children_taxa) == 1 and taxon in children_taxa): #\
                    #and len(node.children) > 1 \
       #or \
       #(len(node.children) == 1 and taxon in children_taxa):
        # Should rather check if taxon in children_taxa.
        if getattr(node, 'D', None) == 'N':
            warnings.warn("The node %r -> %s "\
                  "is marked as a speciation but was infered as a duplication"\
                  " after taxon information." % \
                  (node.name, [ch.name for ch in node.children]))
        if len(children_taxa) > 1:
            warnings.warn("The node %r -> %s "\
                  "is a duplication + a speciation. Not " \
                  "truly reconciled tree." % \
                  (node.name, [ch.name for ch in node.children]))
        elif len(node.children) == 1:
            warnings.warn("The node %r -> %s "\
                  " is a duplication/speciation with one " \
                  "descendant." % \
                  (node.name, [ch.name for ch in node.children]))
        return 'dup'
    else:
        ### FIXME
        if len(children_taxa)==1:
            msg = "The node %r -> %s" % \
                  (node.name, [ch.name for ch in node.children])
            if len(node.children) > 1:
                msg += " is a duplication + a speciation. Not " \
                       "truly reconciled tree."
            #else:
            #    msg += " is a duplication/speciation with one " \
            #           "descendant."
                #if not node.is_root(): event = 'dup'
                # TODO: check this in the duplication block above.
                warnings.warn(msg)

        return 'spe'


class GenetreeDrawer(object):
    """Draw a gene tree inside a species tree"""

    ensembl_version = ENSEMBL_VERSION
    phyltreefile = "/users/ldog/glouvel/ws_alouis/GENOMICUS_SVN/data{0}/" \
                   "PhylTree.Ensembl.{0}.conf"
    
    def __init__(self, phyltreefile=None, ensembl_version=None,
                 colorize_clades=None, commonname=False, latinname=False,
                 treebest=False, show_cov=False):
        if ensembl_version: self.ensembl_version = ensembl_version
        if phyltreefile: self.phyltreefile = phyltreefile
        self.get_taxon = get_taxon_treebest if treebest else get_taxon
        self.commonname = commonname
        self.latinname = latinname

        self.phyltree = PhylTree.PhylogeneticTree(self.phyltreefile.format(
                                                        self.ensembl_version))
        
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
                self.legend_clades.append(patches.Patch(color=cmap(i),
                                                         label=clade))
                self.colorize_species.update(**{sp: cmap(i) for sp in 
                                                self.phyltree.species[clade]})

        self.ancgene2sp = re.compile(r'(' + r'root|'
                        + r'|'.join(list(self.phyltree.listSpecies) +
                                    sorted(self.phyltree.listAncestr,
                                           key=lambda a: len(a),
                                           reverse=True)).replace(' ','\.')
                        + r')(.*)$')
        
        self.taxa = set(self.phyltree.allNames)
        # Cleared and refilled if a genetree is given

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
            self.genetreename = os.path.splitext(os.path.basename(filename))[0]

        #self.genetree = ete3.Tree(filename, format=format)
        with open(filename) as gt_input:
            # Allow reading multiple trees from a single input file.
            genetrees = [ete3.Tree(gt+';', format=format) \
                            for gt in gt_input.read().split(';') \
                            if gt.rstrip()]
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
        root = "root" if "root" in roots else self.phyltree.lastCommonAncestor(list(roots))
        self.taxa.add(root)
        
    def draw_species_tree(self, figsize=None, angle_style=0, branch_width=0.8,
                          colorize_clades=None, ages=False):
        """Init figure + draw branches of the species tree.
        
        branch_width: proportion of vertical space between two branches taken
                      by the branch polygon.
        angle_style: See the parser help."""
        self.species_coords   = {}
        self.species_branches = {}

        if figsize is None: figsize = (8, 8) # 20x20cm
        self.fig, ax0 = plt.subplots(figsize=figsize) #frameon=False) # set figsize later
        #ax0 = self.fig.add_axes([0.1,0.1,0.9,0.9]) #, adjustable='box-forced')
        #ax0 = self.fig.add_axes([0.1, 0.1, 0.9, 0.9])
        ax0.axis('off')
        # TODO: if ages: keep x axis with xlabel "age"

        ymin = 0
        any_show_cov = False
        for parent, (px, py, _), child, (cx, cy, _) in \
              iter_species_coords(self.phyltree, self.taxa, angle_style, ages):

            print(parent, child)
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
            ha = 'center' if cx < 0 else 'left'
            va = 'bottom' if cx < 0 else 'top'
            if cx == 0:
                cx += 0.1 # Add some padding
            # data-specific coloring (low-coverage)
            alpha = 1
            if self.show_cov and child in getattr(self.phyltree, "lstEsp6X", ()):
                alpha = 0.6
                any_show_cov = True
            elif self.show_cov and child in getattr(self.phyltree, "lstEsp2X", ()):
                alpha = 0.3
                any_show_cov = True

            bgcolor = self.colorize_species.get(child, '#ffffff00')

            if cx >= 0:
                # It is a species, add the common name
                if not self.latinname:
                    common_name = get_common_name(self.phyltree, child)
                    if self.commonname:
                        child = common_name
                    elif common_name:
                        child += ', %s' % common_name

            ax0.text(cx, cy, child, ha=ha, va=va, fontsize='x-small',
                     fontstyle='italic', family='serif', alpha=alpha,
                     backgroundcolor=bgcolor)

        # include root.
        self.species_coords[parent] = (px, py)
        ax0.text(px, py, parent, ha='right', fontsize='x-small',
                 fontstyle='italic', family='serif')
        
        # Add legend in case of coverage information
        if self.show_cov and any_show_cov:
            legend_cov = ax0.legend(handles=self.legend_coverage,
                                    loc="upper left",
                                    bbox_to_anchor=(0, 1),
                                    title="Sequencing coverage",
                                    prop={'size': 'x-small', 'family': 'serif',
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
                       prop={'size': 'x-small', 'family': 'serif', 'style': 'italic'},
                       facecolor='inherit')
            if self.show_cov and any_show_cov: ax0.add_artist(legend_cov)

        ax0.set_xlim(px, 1)
        ax0.set_ylim(ymin - 1, 1)
        self.ax0 = ax0


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

        def iter_missing_branches(child_taxon, expected_children_taxa):
            """When intermediate speciation nodes are missing between a taxon and
            descendant taxa, list the implicit branches (a pair of a taxon and its child
            taxon)."""
            print("  Missing between %s and %s: %s" % \
                    (child_taxon, expected_children_taxa,
                     child_taxon not in expected_children_taxa),
                    file=sys.stderr)
            while child_taxon not in expected_children_taxa:
                tmp_taxon = self.phyltree.parent[child_taxon].name
                print("  - get intermediate parent of %s: %s" % \
                        (child_taxon, tmp_taxon), file=sys.stderr)
                while tmp_taxon not in phylsubtree:
                    try:
                        tmp_taxon = self.phyltree.parent[tmp_taxon].name
                    except KeyError as err:
                        print(phylsubtree, file=sys.stderr)
                        err.args += (child_taxon,)
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
                    print('Expected (spe): %s ; got: %s (from %s)' % \
                            (expected_children_taxa, children_taxa, taxon),
                            file=sys.stderr)

                    tmp_children_taxa = []
                
                    # if intermediate speciation nodes are missing, count the
                    # implied deletions.
                    for child_taxon in children_taxa:

                        # It works when child_taxon is *descendant* of expected_taxa.
                        for tmp_taxon, child_taxon in iter_missing_branches(child_taxon,
                                                          expected_children_taxa):

                            tmp_expected_children_taxa = phylsubtree[tmp_taxon]

                            deletion_count[node] += orient_deletion(set((tmp_taxon,)), tmp_expected_children_taxa)
                            print('+1 deletion', file=sys.stderr)
                        tmp_children_taxa.append(child_taxon)

                    deletion_count[node] += orient_deletion(children_taxa, expected_children_taxa)

            elif event == 'dup':
                if len(children_taxa) > 1:
                    for child in node.children:
                        child_taxon = cached_taxa[child]
                        if child_taxon != taxon:
                            print('Expected (dup): %s ; got: %s (from %s)' % \
                                    (taxon, children_taxa, taxon),
                                    file=sys.stderr)
                            for tmp_taxon, _ in iter_missing_branches(child_taxon,
                                                        set((taxon,))):
                                tmp_expected_children_taxa = phylsubtree[tmp_taxon]
                                deletion_count[child] += orient_deletion(set((tmp_taxon,)), tmp_expected_children_taxa)

                node.children = sorted(node.children, key=deletion_count.get)

        self.genetrees = sorted(self.genetrees, key=deletion_count.get)


    def set_gene_coords(self):
        ### TODO: enable repeted calls to set multiple gene trees
        ### TODO: propagate deleted lineage : draw invisible lines, just to
        ###       consume vertical space.
        ### TODO: genes whose children got deleted in one lineage should be
        ###       placed away from this lineage (avoid line crossing)
        
        self.gene_coords = {} # {species: [gene list]}
        #self.dup_branchings = [] # (nodeid/genename, x, x_child1, x_child2, y_child1, y_child2)
        #self.spe_branchings = []
        self.branchings = {} # {nodeid: [species, x, y, dup/spe,
                             #           (children), {node features}], ...]
        
        # *temporary* tree structure holding nodes that are inbetween 2 speciations
        interspecies_trees = {}

        # x and y are relative to the branch considered. (between 0 and 1)

        self.minimize_gene_branch_crossing()

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
            node_features = {ft: getattr(node, ft) for ft in node.features}

            if event == 'leaf':
                taxon_gene_coords.append(nodeid)
                node_y = len(taxon_gene_coords) - 1 # This list is being
                                                    # extended later on.
                interspecies_trees[nodeid] = {'taxon': taxon,
                                               'ndup': 0,
                                               'x': 1,
                                               'y': node_y}
                self.branchings[nodeid] = [taxon, 0, node_y, 'leaf', [], node_features]
            elif event == 'dup':

                # It is a **duplication**
                children_ys   = []
                children_ndup = []
                for ch in node.children:
                    children_ys.append(interspecies_trees[ch.id]['y'])
                    children_ndup.append(interspecies_trees[ch.id]['ndup'])

                # TODO: tilt according to the parent speciation node.
                # should fit in the triangle defined by parent and descendants.
                #node_y    = (min(children_ys) + max(children_ys)) / 2
                
                # weighted average of the children positions.
                node_y    = sum(cy*(cd+1) for cy,cd in \
                                           zip(children_ys, children_ndup))
                node_y   /= (sum(children_ndup) + len(children_ndup))
                node_ndup = max(children_ndup) + 1

                interspecies_trees[nodeid] = {
                        'taxon': taxon,
                        'ndup': node_ndup,
                        'x': 1, # modified later
                        'y': node_y}
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
                            'y': node_y}
                else:
                    # If it is the root, place it in the middle of the branch.
                    node_x = 1 - node_ndup/(node_ndup+1)
                    interspecies_trees[nodeid]['x'] = node_x

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

                            self.branchings[nextnode.id] = [
                                                            nextnode_taxon,
                                                            nextnode_x,
                                                            nextnode_y,
                                                            'dup',
                                                            nextnode_ch,
                                                            node_features]
                            nextnodes.extend(nextnode.children)

                        interspecies_trees.pop(nextnode.id)

        #assert not interspecies_trees, \
        #        "Temporarily stored nodes have not been processed: %s" % interspecies_trees
        print("Root:", interspecies_trees)


    def draw_gene_tree(self, extratitle='', genenames=False, branch_width=0.8):
        print(' --- Drawing genetree ---')
        # Duplicate the species tree axis to separate the plotting
        self.ax1 = self.ax0.twinx()
        self.ax1.set_ylim(self.ax0.get_ylim())
        self.ax1.axis('off')
        title = self.genetreename
        if extratitle: title += ' -- %s' % extratitle
        self.ax1.set_title(title)
        
        self.real_gene_coords = {}

        cmap = plt.get_cmap('Dark2', 10) # TODO: as many as duplications

        color_index = 0

        # Gene node names to be displayed (dup, leaf, spe)
        genenames=set(genenames.split(',')) if genenames else set()

        for node in (n for genetree in self.genetrees \
                     for n in genetree.traverse('postorder')):
            nodeid = node.id
            
            species, rel_x, rel_y, event, children, _ = self.branchings[nodeid]
            pos_list = self.gene_coords[species]
            nranks = len(pos_list) + 1
            children_real_coords = [self.real_gene_coords[ch] for ch in children]
            children_features = [self.branchings[ch][5] for ch in children]
            
            #print(nodeid, event, species, children)
            if event == 'dup':

                # coordinates of the **species** branch vector
                # Dx = x (taxon) - x (ancestor)
                parent_sp, Dx, Dy = self.species_branches[species]
                parent_x, parent_y = self.species_coords[parent_sp]
                real_x = parent_x + Dx * rel_x
                real_y = parent_y + Dy * rel_x - branch_width * (rel_y+1)/nranks
                
                if not node.is_root():
                    nodecolor = 'red'
                    branches_color = cmap(color_index)
                    color_index = (color_index + 1) % cmap.N
                else:
                    nodecolor='black'
                    branches_color = 'black'

                nch = len(children)
                children_rel_ys = [self.branchings[ch][2] for ch in children]

                assert all((ch_rel_y >= 0) for ch_rel_y in children_rel_ys)
                if any((children_rel_ys[i+1] - children_rel_ys[i] < 0) for i in range(nch-1)): 
                    warnings.warn("Children's relative Y not sorted! (%s: %s)" % (event, node.name))
                elif any((children_rel_ys[i+1] - children_rel_ys[i] == 0) for i in range(nch-1)):
                    warnings.warn("Some children's relative Y are identical! (%s: %s)" % (event, node.name))

                for ch, ch_rel_y, (ch_real_x, ch_real_y), ch_ft in \
                        zip(children, children_rel_ys, children_real_coords,
                            children_features):

                    # Draw thicker line when it represents a paralogy
                    linewidth = 4 if ch_ft.get('P') == 'True' else 1
                    
                    delta_y = (rel_y - ch_rel_y)/nranks * branch_width
                    fork_coords = [(ch_real_x, ch_real_y),
                                   (real_x, real_y + delta_y),
                                   (real_x, real_y)]
                    u_fork_finger = patches.PathPatch(
                                    Path(fork_coords, [MOVETO, CURVE3, CURVE3]),
                                    fill=False,
                                    edgecolor=branches_color,
                                    alpha=0.5,
                                    linewidth=linewidth,
                                    joinstyle='round')  # don't see a change
                    self.ax1.add_patch(u_fork_finger)

            else:  # event == 'spe' or 'leaf'
                real_x, real_y = self.species_coords[species]
                pos = pos_list.index(nodeid) + 1
                real_y -= branch_width * pos/nranks
                nodecolor = 'blue'
                #nodecolor = 'none'

                for ch, (ch_real_x, ch_real_y), ch_ft in \
                        zip(children, children_real_coords, children_features):
                    linewidth = 4 if ch_ft.get('P') == 'True' else 1
                    v_fork_finger = patches.PathPatch(
                                             Path([(real_x, real_y),
                                                   (ch_real_x, ch_real_y)],
                                                  [MOVETO, LINETO]),
                                             edgecolor='black',
                                             alpha=0.5,
                                             linewidth=linewidth,
                                             joinstyle='round')
                    self.ax1.add_patch(v_fork_finger)
                                  #(':' if node.is_root() else '-'),
            
            #if event != 'leaf':
            if event == 'dup':  # This line is still here because historically it was.
                # Add a *dot*
                self.ax1.plot((real_x,), (real_y,), '.', color=nodecolor, alpha=0.5)
                #self.ax1.text(real_x, real_y, species, fontsize='xxx-small',
                #              color=nodecolor)
            if event in genenames:
                self.ax1.text(real_x, real_y, node.name, alpha=0.5,
                              fontsize='xx-small',
                              ha=('left' if event=='leaf' else 'right'),
                              va='top')

            self.real_gene_coords[nodeid] = (real_x, real_y)
            #self.fig.draw(self.fig.canvas.get_renderer())
            #plt.show()
            #input('Press Enter to continue.')
            ### TODO: add onpick action: display node name

    def draw(self, genetree, extratitle='', angle_style=0, ages=False,
             figsize=None, genenames=False):
        """Once phyltree is loaded, perform all drawing steps."""
        self.load_reconciled_genetree(genetree)
        self.draw_species_tree(figsize=figsize, angle_style=angle_style, ages=ages)
        self.set_gene_coords()
        self.draw_gene_tree(extratitle, genenames=genenames)


def check_extracted(genetrees, output):
    all_outputs = [output.format(genetree=gt) for gt in genetrees]
    new_genetrees = [gt for gt in genetrees if not os.path.exists(output.format(genetree=gt))]
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
    assert os.path.exists(datadir)

    if edited:
        # Take gene tree from Genomicus
        treeforestfile = "/users/ldog/glouvel/ws_alouis/GENOMICUS_SVN/data%d/"\
                         "GoodThreshold/tree.4F.cut.bz2" % ensembl_version
        withAncGenesNames = True
        field = 'family_name'
        output = os.path.join(datadir, '{genetree}', '{genetree}.nwk')
        fix_suffix = True
    else:
        # Take gene tree from Ensembl
        treeforestfile = "/users/ldog/glouvel/ws_alouis/GENOMICUS_SVN/data%d/"\
                         "tree.1.ensembl.bz2" % ensembl_version
        withAncGenesNames = False
        field = 'tree_name'
        output = os.path.join(datadir, '{genetree}', '{genetree}_ensembl.nwk')
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
    
    prune_outdir = os.path.join(datadir, gt_format, subtrees_dir)

    for gt in genetrees:
        p_outdir = prune_outdir.format(gt)
        print(p_outdir, end=' ')
        if not os.path.exists(p_outdir):
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


def run(outfile, genetrees, angle_style=0, ensembl_version=ENSEMBL_VERSION, 
        phyltreefile=None, colorize_clades=None, commonname=False,
        latinname=False, treebest=False, show_cov=False, ages=False,
        genenames=False):
    #global plt

    figsize = None
    gd = GenetreeDrawer(phyltreefile=phyltreefile,
                        ensembl_version=ensembl_version,
                        colorize_clades=colorize_clades,
                        commonname=commonname,
                        latinname=latinname,
                        treebest=treebest,
                        show_cov=show_cov)
    display = lambda: plt.show() # Display function for shell or notebook usage
    if __name__=='__main__' and outfile == '-':
        try:
            plt.switch_backend('Qt4Agg')
        except ImportError:
            plt.switch_backend('TkAgg')
        #mpl.use('Qt4Agg')
        #from importlib import reload; reload(plt)
    elif outfile.endswith('.pdf'):
        pdf = PdfPages(outfile)
        figsize = (8.2, 11.7)
        display = lambda: (pdf.savefig(bbox_inches='tight', papertype='a4'),
                           plt.close())
    else:
        assert len(genetrees) <= 1, "multipage output only supported for pdf"
        display = lambda: (plt.savefig(outfile, bbox_inches='tight'),
                           plt.close())

    for genetree in genetrees:
        genetree, *extratitles = genetree.split(',')
        extratitle = ', '.join(extratitles)
        print('INPUT FILE:', genetree, '(%s)' % extratitle)
        gd.draw(genetree, extratitle, angle_style=angle_style, ages=ages,
                figsize=figsize, genenames=genenames)

        display()

    if not genetrees:
        gd.draw_species_tree(figsize=figsize, angle_style=angle_style, ages=ages)
        display()

    if outfile.endswith('.pdf'):
        pdf.close()
    #plt.show()
    return gd


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description=__doc__,
                        #formatter_class=CustomHelpFormatter)
                        formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('outfile', help=("pdf file, or '-'. If '-', will use "
                                         "Qt to display the figure."))
    parser.add_argument('genetrees', nargs='*', default=[],
        help=("must be a genetree (nwk format with internal nodes labelling) "
            "reconciled with species tree, or a genetree formatted like "
            "`TreeBest` output. '-' means standard input."))
    parser.add_argument('--fromfile', action='store_true',
                        help='take genetree paths and description from a file')
    parser.add_argument('-e', '--ensembl-version', type=int,
                        default=ENSEMBL_VERSION)
    parser.add_argument('-p', '--phyltreefile', default=PHYLTREEFILE,
                        help='species tree in phyltree or newick format [%(default)s]')
    
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
    
    parser.add_argument('-c', '--colorize-clade', action='append',
                        dest='colorize_clades',
                        help='species in these clades will have a specific color')
    parser.add_argument('-t', '--treebest', action='store_true',
                        help='The input genetree is a treebest output')
    parser.add_argument('-s', '--show-cov', action='store_true',
                        help='Show genome coverage information (grey shading)')
    parser.add_argument('-A', '--ages', action='store_true',
                        help='Place species nodes at their real age.')
    parser.add_argument('-g', '--genenames',
            help='Display gene names: \n' \
                 '- comma-sep list of "leaf", "dup", "spe";\n'
                 '- "all" (identical to "leaf,dup,spe").')
    
    #parser.add_argument('-m', '--multiple-pdfs', action='store_true',
    #                    help='output one pdf file per genetree. [NOT implemented]')
    args = parser.parse_args()
    dictargs = vars(args)
    #if not dictargs.get('genetrees'):
    #    dictargs['genetrees'] = [TESTTREE]
    if dictargs.pop('fromfile'):
        assert len(dictargs['genetrees']) == 1
        genetreelistfile = dictargs.pop('genetrees')[0]
        with open(genetreelistfile) as stream:
            dictargs['genetrees'] = [line.rstrip() for line in stream
                                        if not line.startswith('#')]
            
    # TODO: add into run()
    #ANCGENE2SP = re.compile(r'([A-Z][A-Za-z_.-]+)%s' % ancgene_regex)

    gd = run(**dictargs)

