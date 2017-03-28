#!/usr/bin/env python3


"""
Draw a gene tree inside a species tree.

USAGE:

./genetree_drawer.py <genetreefile>

ARGUMENTS:
  - outputfile:   pdf file, or '-'. If '-', will use Qt to display the figure.
  - genetreefile: must be a genetree (nwk format with internal nodes labelling)
                  reconciled with species tree.
"""

import sys
import os
import os.path
import re
import argparse
import numpy as np
import matplotlib as mpl
#mpl.use('Agg')
if __name__=='__main__': mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.backends.backend_pdf import PdfPages, FigureCanvas
import ete3

import LibsDyogen.myPhylTree as PhylTree
#import LibsDyogen.myProteinTree as ProteinTree

from glou_duphist import dfw_descendants_generalized, ladderize
from codeml.select_leaves_from_specieslist import convert_gene2species


ENSEMBL_VERSION = 85
ANCGENE2SP = re.compile(r'([A-Z][A-Za-z_.-]+)ENS')

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

def get_taxon(node, ancgene2sp, ensembl_version=ENSEMBL_VERSION):
    """from a gene name in my newick gene trees, find the taxon:
        either:
            - node is a leaf (e.g ENSMUSG00...)
            - node is internal (e.g Mus.musculusENSGT...)"""
    if node.is_leaf():
        taxon = convert_gene2species(node.name, ensembl_version)
    else:
        try:
            taxon = ancgene2sp.match(node.name).group(1).replace('.', ' ')
        except AttributeError:
            raise RuntimeError("Can not match species name in %r" % node.name)
    return taxon


### Unused function
def get_taxa_set(etegenetree, ancgene2sp, ensembl_version=ENSEMBL_VERSION):
    taxa_set = set()
    for node in etegenetree.traverse():
        taxon = get_taxon(node, ancgene2sp, ensembl_version)
        taxa_set.add(taxon)
    return taxa_set


def walk_phylsubtree(phyltree, taxa):
    """Return an iterator that will progress on the tree from leaves to root.
    
    It also rotates forks such that the displayed tree will show branches
    that are most distant from the root (in number of nodes) on one side.

    Taxa returned by the iterator are space-separated.
    """
    root, subtree = phyltree.getSubTree(taxa)
    #if lower_root and phyltree.parent.get(root):
    #    newroot = phyltree.parent[root].name
    #    subtree[newroot] = [root]
    #    root = newroot

    # reorder branches in a visually nice manner:
    ladderize(subtree, root)
    get_children = lambda tree, node: tree.get(node, [])
    dfw = dfw_descendants_generalized(subtree, get_children,
                                      include_leaves=False, queue=[root])
    return reversed(list(dfw))


def iter_species_coords(phyltree, taxa, angle_style=0):
    """Assign a pair x,y of coordinates for each node in the tree.
    Yield (parent name, parent_xy, child name, child_xy).
    
    Just the topology. Branch lengths are ignored.

    lower_root: add the ancestor of the actual root to the figure.
    """
    coords = {}
    y0 = 0
    for parent, children in walk_phylsubtree(phyltree, taxa):
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

        if angle_style == 0:
            # Along the X axis: move one step to the left
            # Along the Y axis: take the average of the children Y coordinates.
            parent_x = min(children_xs) - 1
            parent_y = sum(children_ys) / len(children_ys)
        else:
            # TODO: dx < 0 ?
            dy = max(children_ys) - min(children_ys)
            if dy == 0: dy = 1 # when only one child
            dx = max(children_xs) - min(children_xs)

            if angle_style == 1:
                # branches are all at 45 degrees
                step = (dy - dx) / 2
                parent_x = min(children_xs) - step
                parent_y = max(children_ys) - step
            elif angle_style == 2:
                step = dy / (2+dx)
                parent_x = min(children_xs) - 1
                parent_y = max(children_ys) - step
            else:
                raise RuntimeError("Invalid angle_style value: %r" % angle_style)

        coords[parent] = (parent_x, parent_y)

        for child in children:
            yield parent, coords[parent], child, coords[child]


class GenetreeDrawer(object):
    """Draw a gene tree inside a species tree"""

    ensembl_version = ENSEMBL_VERSION
    phyltreefile = "/users/ldog/glouvel/ws_alouis/GENOMICUS_SVN/data{0}/" \
                   "PhylTree.Ensembl.{0}.conf"
    
    def __init__(self, phyltreefile=None, ensembl_version=None,
                 colorize_clades=None):
        if ensembl_version: self.ensembl_version = ensembl_version
        if phyltreefile: self.phyltreefile = phyltreefile

        self.phyltree = PhylTree.PhylogeneticTree(self.phyltreefile.format(
                                                        self.ensembl_version))
        
        self.colorize_species  = {}
        if colorize_clades is not None:
            cmap = plt.get_cmap('Set3', len(colorize_clades))
            for i, clade in enumerate(colorize_clades):
                self.colorize_species.update(**{sp: cmap(i) for sp in 
                                                self.phyltree.species[clade]})

        self.ancgene2sp = re.compile(r'('
                        + r'|'.join(list(self.phyltree.listSpecies) +
                                    sorted(self.phyltree.listAncestr,
                                           key=lambda a: len(a),
                                           reverse=True)).replace(' ','\.')
                        + r')(.*)$')

    def load_reconciled_genetree(self, filename, format=1, genetreename=None):
        """Load gene tree with all species nodes present.
        newick format with internal node labelling"""

        if genetreename:
            self.genetreename = genetreename
        else:
            self.genetreename = os.path.splitext(os.path.basename(filename))[0]
        self.genetree = ete3.Tree(filename, format=format)
        self.genetree.ladderize()
        # add only the meaningful taxa (not those with one child and age = 0)
        root = get_taxon(self.genetree, self.ancgene2sp, self.ensembl_version)
        alldescendants = self.phyltree.allDescendants[root]
        self.taxa = set()
        for taxon in alldescendants:
            if not (len(self.phyltree.items.get(taxon, [])) == 1 and
                    self.phyltree.ages[taxon] == 0):
                self.taxa.add(taxon)
        # Add the branch leading to the current root (if duplications in this branch)
        lower_root = self.phyltree.parent[root].name
        # This check is not especially necessary.
        while len(self.phyltree.items.get(lower_root, [])) == 1 \
                and self.phyltree.ages[lower_root] == 0:
            lower_root = self.phyltree.parent[lower_root].name

        self.taxa.add(lower_root)

        rerooted_genetree = ete3.TreeNode(name=self.genetree.name.replace(root,
                                                                lower_root))
        rerooted_genetree.add_child(child=self.genetree)
        self.genetree = rerooted_genetree



    def draw_species_tree(self, figsize=None, angle_style=0, branch_width=0.8,
                          colorize_clades = None):
        """Init figure + draw branches of the species tree.
        
        branch_width: proportion of vertical space between two branches taken
                      by the branch polygon.
        angle_style: - 0: parent node y = mean of children nodes y
                     - 1: branches always at 45 degrees
                     - 2: parent node positioned at x-1 but branch angles are equal"""
        self.species_coords   = {}
        self.species_branches = {}

        if figsize is None: figsize = (8, 8) # 20x20cm
        self.fig, ax0 = plt.subplots(figsize=figsize) #frameon=False) # set figsize later
        #ax0 = self.fig.add_axes([0.1,0.1,0.9,0.9]) #, adjustable='box-forced')
        #ax0 = self.fig.add_axes([0.1, 0.1, 0.9, 0.9])
        ax0.axis('off')

        ymin = 0
        for parent, (px, py), child, (cx, cy) in iter_species_coords(self.phyltree,
                                                                     self.taxa,
                                                                     angle_style):
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
            if child in self.phyltree.lstEsp6X: alpha = 0.6
            if child in self.phyltree.lstEsp2X: alpha = 0.3

            bgcolor = self.colorize_species.get(child, '#ffffff00')
            ax0.text(cx, cy, child, ha=ha, va=va, fontsize='x-small',
                     fontstyle='italic', family='serif', alpha=alpha,
                     backgroundcolor=bgcolor)

        # include root.
        self.species_coords[parent] = (px, py)
        ax0.text(px, py, parent, ha='right', fontsize='x-small',
                 fontstyle='italic', family='serif')
        
        ax0.set_xlim(px, 1)
        ax0.set_ylim(ymin - 1, 1)
        self.ax0 = ax0


    def set_gene_coords(self):
        ### TODO: enable repeted calls to set multiple gene trees
        ### TODO: propagate deleted lineage : draw invisible lines, just to
        ###       consume vertical space.
        ### TODO: genes whose children got deleted in one lineage should be
        ###       placed away from this lineage (avoid line crossing)
        
        self.gene_coords = {} # {species: [gene list]}
        #self.dup_branchings = [] # (genename, x, x_child1, x_child2, y_child1, y_child2)
        #self.spe_branchings = []
        self.branchings = {} # {genename: [species, x, y, dup/spe, (children)], ...]
        
        # temporary tree structure holding nodes that are inbetween 2 speciations
        interspecies_trees = {}

        # x and y are relative to the branch considered. (between 0 and 1)

        for nodeid, node in enumerate(self.genetree.traverse('postorder')):
            # set up a *unique* node ID, to avoid relying on node.name.
            node.id = nodeid

            taxon = get_taxon(node, self.ancgene2sp, self.ensembl_version)
            taxon_gene_coords = self.gene_coords.setdefault(taxon, [])

            if node.is_leaf():
                taxon_gene_coords.append(node.id)
                node_y = len(taxon_gene_coords) - 1
                interspecies_trees[node.id] = {'taxon': taxon,
                                               'ndup': 0,
                                               'x': 1,
                                               'y': node_y}
                self.branchings[node.id] = [taxon, 0, node_y, 'leaf', []]
            else:
                children_taxa = set(interspecies_trees[ch.id]['taxon'] for ch
                                    in node.children)
                if len(children_taxa) == 1 and len(node.children) > 1:
                    # It is a duplication
                    children_ys   = []
                    children_ndup = []
                    for ch in node.children:
                        children_ys.append(interspecies_trees[ch.id]['y'])
                        children_ndup.append(interspecies_trees[ch.id]['ndup'])

                    node_y    = (min(children_ys) + max(children_ys)) / 2
                    node_ndup = max(children_ndup) + 1

                    interspecies_trees[node.id] = {
                            'taxon': taxon,
                            'ndup': node_ndup,
                            'x': 1, # modified later
                            'y': node_y}
                else:
                    # It is a speciation
                    if taxon in children_taxa:
                        msg = "WARNING: the node %r -> %s" % \
                              (node.name, [ch.name for ch in node.children])
                        if len(node.children) > 1:
                            msg += " is a duplication + a speciation. Not " \
                                   "truly reconciled tree."
                        else:
                            msg += " is a duplication/speciation with one " \
                                   "descendant."
                        print(msg, file=sys.stderr)
                    taxon_gene_coords.append(node.id)

                    node_y = len(taxon_gene_coords) - 1
                    interspecies_trees[node.id] = {
                            'taxon': taxon,
                            'ndup': 0,
                            'x': 0,
                            'y': node_y}
                    self.branchings[node.id] = [taxon, 0, node_y, 'spe',
                                                [ch.id for ch in
                                                    node.children]]

                    for ch in node.children:
                        n_steps = interspecies_trees[ch.id]['ndup'] + 1

                        delta_x = 1 / n_steps

                        # climb up this subtree until next speciation
                        nextnodes = [ch]
                        while nextnodes:
                            nextnode = nextnodes.pop(0)
                            try:
                                nextnode_ndup = interspecies_trees[nextnode.id]['ndup']
                            except KeyError as err:
                                err.args += ('Non unique node name.',)
                                raise
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
                                                                nextnode_ch]
                                nextnodes.extend(nextnode.children)

                            interspecies_trees.pop(nextnode.id)
        print("root:", interspecies_trees)

                                                                
    def draw_gene_tree(self, extratitle='', branch_width=0.8):
        print(' --- Drawing genetree ---')
        # Duplicate the species tree axis to separate the plotting
        #if hasattr(self, 'ax1'):
        #    self.ax1.clear()
        #else:
        self.ax1 = self.ax0.twinx()
        self.ax1.set_ylim(self.ax0.get_ylim())
        self.ax1.axis('off')
        title = self.genetreename
        if extratitle: title += ' -- %s' % extratitle
        self.ax1.set_title(title)
        
        self.real_gene_coords = {}

        cmap = plt.get_cmap('Dark2', 10) # TODO: as many as duplications

        color_index = 0

        #seen_genenames = set()
        for node in self.genetree.traverse('postorder'):
            genename = node.id
            #if genename in seen_genenames:
            #    print('WARNING: %r already seen' % genename)
            #else:
            #    seen_genenames.add(genename)
            
            species, rel_x, rel_y, event, children = self.branchings[genename]
            pos_list = self.gene_coords[species]
            nranks = len(pos_list) + 1
            children_real_coords = [self.real_gene_coords[ch] for ch in children]
            
            #print(genename, event, species, children)
            if event == 'dup':

                children_rel_ys = [(self.branchings[ch][2]) for ch in children]

                parent_sp, Dx, Dy = self.species_branches[species]
                parent_x, parent_y = self.species_coords[parent_sp]
                real_x = parent_x + Dx * rel_x
                real_y = parent_y + Dy * rel_x - branch_width * (rel_y+1)/nranks
                nodecolor = 'red'

                branches_color = cmap(color_index)
                color_index = (color_index + 1) % cmap.N

                delta_ys = []
                for ch, ch_rel_y, (ch_real_x, ch_real_y) in \
                    zip(children, children_rel_ys, children_real_coords):
                    delta_y = (rel_y - ch_rel_y)/nranks * branch_width
                    delta_ys.append(delta_y)

                    self.ax1.plot((ch_real_x, real_x),
                                  (ch_real_y, real_y + delta_y),
                                  color=branches_color, alpha=0.5)
                
                # plot the vertical line of the fork.
                self.ax1.plot((real_x, real_x),
                              (real_y + min(delta_ys), real_y + max(delta_ys)),
                              color=branches_color, alpha=0.5)

            else:
                real_x, real_y = self.species_coords[species]
                pos = pos_list.index(genename) + 1
                real_y -= branch_width * pos/nranks
                nodecolor = 'blue'

                for ch_real_x, ch_real_y in children_real_coords:
                    self.ax1.plot((real_x, ch_real_x), (real_y, ch_real_y),
                                  color='black', alpha=0.5)
            
            if event != 'leaf':
                self.ax1.plot((real_x,), (real_y,), '.', color=nodecolor, alpha=0.5)
                #self.ax1.text(real_x, real_y, species, fontsize='xxx-small',
                #              color=nodecolor)

            self.real_gene_coords[genename] = (real_x, real_y)
            #self.fig.draw(self.fig.canvas.get_renderer())
            #plt.show()
            #input('Press Enter to continue.')
            ### TODO: add onpick action: display node name

    def draw(self, genetree, extratitle='', angle_style=0, figsize=None):
        """Once phyltree is loaded, perform all drawing steps."""
        self.load_reconciled_genetree(genetree)
        self.draw_species_tree(figsize=figsize, angle_style=angle_style) # Currently, must be called after load_reconciled
        self.set_gene_coords()
        self.draw_gene_tree(extratitle)


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
        colorize_clades=None):
    #global plt

    figsize = None
    gd = GenetreeDrawer(ensembl_version=ensembl_version,
                        colorize_clades=colorize_clades)
    if outfile != '-':
        pdf = PdfPages(outfile)
        figsize = (8.2, 11.7)
    elif __name__=='__main__':
        plt.switch_backend('Qt4Agg')
        #mpl.use('Qt4Agg')
        #from importlib import reload; reload(plt)

    for genetree in genetrees:
        genetree, *extratitles = genetree.split(',')
        extratitle = ', '.join(extratitles)
        print('INPUT FILE:', genetree, '(%s)' % extratitle)
        gd.draw(genetree, extratitle, angle_style=angle_style, figsize=figsize)

        if outfile == '-':
            plt.show()
            figsize = gd.fig.get_size_inches()
        else:
            pdf.savefig(bbox_inches='tight', papertype='a4')
            plt.close()
    if outfile != '-':
        pdf.close()
    #plt.show()
    return gd


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__,
                        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('outfile')
    parser.add_argument('genetrees', nargs='*')
    parser.add_argument('--fromfile', action='store_true',
                        help='take genetree paths and description from a file')
    parser.add_argument('-e', '--ensembl-version', type=int,
                        default=ENSEMBL_VERSION)
    parser.add_argument('-a', '--angle-style', type=int, choices=[0,1,2],
                        default=0,
                        help=("0: parent node y = mean of children nodes y\n"
                              "1: branches always at 45 degrees\n"
                              "2: parent node positioned at x-1 but branch "
                              "angles are equal"))
    parser.add_argument('-c', '--colorize-clade', action='append',
                        dest='colorize_clades',
                        help='species in these clades will have a specific color')
    args = parser.parse_args()
    dictargs = vars(args)
    if not dictargs.get('genetrees'):
        dictargs['genetrees'] = [TESTTREE]
    elif dictargs.pop('fromfile'):
        assert len(dictargs['genetrees']) == 1
        genetreelistfile = dictargs.pop('genetrees')[0]
        with open(genetreelistfile) as stream:
            dictargs['genetrees'] = [line.rstrip() for line in stream]
            
    # TODO: add into run()
    #ANCGENE2SP = re.compile(r'([A-Z][A-Za-z_.-]+)%s' % ancgene_regex)

    gd = run(**dictargs)

