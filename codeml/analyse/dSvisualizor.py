#!/usr/bin/env python3

# Evolution of `plot_dS.py`. Refactor code using classes.
# Use Pandas features more systematically.

"""Visualise values from given table (reconstructed ages) per taxonomic branch.

USAGE:
    ./dSvisualizor.py command <ages_file> <outfile>
"""

COMMANDS = ['lineage', 'tree', 'scatter', 'violin']

CMD_ARGS = {'lineage': [(('-l', '--lineage'),),
                        (('-p', '--phyltreefile'),),
                        (('-e', '--ensembl-version'),)],
            'tree':    [(('-v', '--vertical'),
                         dict(action='store_true')),
                        (('-p', '--phyltreefile'),),
                        (('-e', '--ensembl-version'),)],
            'scatter': [(('-x',),), (('-y',),), (('--xlim',),), (('--ylim',),)],
            'violin': [(('-x',),), (('-y',),), (('--xlim',),), (('--ylim',),)]}


import sys
import os.path
import re
import bz2
import pickle
import argparse
import matplotlib as mpl
mpl.use('Qt4Agg', warn=False) # for figures to show up when the script is called from the shell
import matplotlib.pyplot as plt
#plt.ion()
import numpy as np
import pandas as pd
import LibsDyogen.myPhylTree as PhylTree


from glou_duphist import dfw_descendants_generalized, ladderize

# Change all black to dark grey
grey10 = '#1a1a1a'
grey45 = '#737373'
grey80 = '#CCCCCC'
mpl.rcParams['text.color'] = grey10
mpl.rcParams['axes.edgecolor'] = grey10
mpl.rcParams['axes.labelcolor'] = grey10
mpl.rcParams['axes.spines.top'] = False
mpl.rcParams['axes.spines.right'] = False
mpl.rcParams['xtick.color'] = grey10
mpl.rcParams['ytick.color'] = grey10
mpl.rcParams['grid.color'] = grey80
mpl.rcParams['grid.linestyle'] = ':'
mpl.rcParams['patch.edgecolor'] = grey10
mpl.rcParams['boxplot.flierprops.markeredgecolor'] = grey10
mpl.rcParams['boxplot.capprops.color'] = grey10
mpl.rcParams['legend.facecolor'] = grey45
mpl.rcParams['legend.framealpha'] = 0.2
#mpl.rcParams['legend.edgecolor'] = grey10
mpl.rcParams['savefig.facecolor'] = 'none'
#mpl.rcParams['savefig.frameon'] = False  #background frame transparent
#mpl.rcParams['savefig.transparent'] = True # all background transparent
                                            # (including ggplot2 style)
#mpl.style.use('ggplot')
pd.set_option('display.max_colwidth', 85)

PHYLTREEFILE = "/users/ldog/glouvel/ws_alouis/GENOMICUS_SVN/data{0}/" \
                   "PhylTree.Ensembl.{0}.conf"
ENSEMBL_VERSION = 85

RE_TAXON = re.compile(r'[A-Z][A-Za-z_.-]+(?=ENSGT)')
#PAT_TAXON = r'^(ENS[A-Z]+G|[A-Z][A-Za-z_.-]+)(ENSGT[0-9]+|)(.*)$'
#PAT_TAXON = r'^([A-Z][A-Za-z_.-]+)(ENSGT[0-9]+)(.*)$'
PAT_TAXON = r'^(ENS[A-Z]+G|[A-Z][A-Za-z_.-]+(?=ENSGT))(|ENSGT[0-9]+)([0-9]+|[.A-Za-z`]*)$'

DEFAULT_NBINS = 50
DEFAULT_AGE_KEY = 'age_dS'

def bin_data(data, bins, binvar=DEFAULT_AGE_KEY, outvar=None):
    """given some data and some bins, return groups by bins."""
    #bins = data_bins[lab.replace('.', ' ')]
    nbins = len(bins)
    #print('NBINS:', nbins)
    if not outvar:
        outvar = data.keys() # select all columns.
    bin_positions = np.digitize(data[binvar], bins)
    # include values that are *equal* to the right-most edge.
    bin_positions[bin_positions == nbins] = nbins-1
    return data[outvar].groupby(bin_positions)

#def bin_data(dataseries, bins):
#    """Group data by bins"""
#    bin_positions = np.digitize(data, bins)
#    # include values that are *equal* to the right-most edge.
#    nbins = len(bins)
#    #bin_positions[bin_positions == nbins+1] = nbins
#    return data.groupby(bin_positions)


def intersect_size(serie, checked_values):
    is_edited = set(serie) & checked_values
    return float(len(is_edited)) / len(serie)


def splitname2taxongenetree(df, name_col="name", index_col=None):
    if index_col:
        names = df.index.get_level_values(index_col).to_series()
        names.index = df.index
    else:
        names = df[name_col]
    splitted_names = names.str.extract(PAT_TAXON, expand=True)
    splitted_names.columns = ['taxon', 'genetree', 'suffix']
    splitted_names['taxon'] = splitted_names['taxon'].str.replace('\.', ' ')
    #print(splitted_names.head(50))
    concat_cols = [col for col in splitted_names.columns if col not in df.columns]
    named = pd.concat([df, splitted_names[concat_cols]], ignore_index=True, axis=1)
    named.columns = df.columns.tolist() + concat_cols
    named.index = df.index
    return named


class DataVisualizor(object):

    ensembl_version = ENSEMBL_VERSION
    phyltreefile = PHYLTREEFILE
    default_nbins = DEFAULT_NBINS

    def load_edited_set(self, treeforest):
        print('Loading the set of edited nodes:', end=' ')
        pickled_file = os.path.basename(os.path.splitext(treeforest)[0]) + \
                        '.editedset.pickle'
        if os.path.exists(pickled_file):
            print('from pickle...')
            with open(pickled_file, 'rb') as pickle_in:
                self.edited_set = pickle.load(pickle_in)
        else:
            print('searching tree...')
            regex = re.compile("'family_name': '(ENSGT[0-9]+[A-Za-z`.]*)',.*'taxon_name': '([A-Z][A-Za-z _.-]+)'")
            edited_list = []
            with bz2.BZ2File(treeforest) as stream:
                for line in stream:
                    line = line.decode()
                    if "'Duplication': 3" in line:
                        fam_name, tax_name = regex.search(line).groups()
                        edited_list.append(tax_name.replace(' ','.') + fam_name)
            self.edited_set = set(edited_list)
            # save for later
            with open(pickled_file, 'wb') as pickle_out:
                pickle.dump(self.edited_set, pickle_out)

    def __init__(self, ages_file, no_edited=None, age_key=DEFAULT_AGE_KEY):
        """Load and format data"""

        self.edited_set = None
        self.not_edited = None
        self.no_edited  = no_edited
        self.age_key = age_key
        self.phyltree = None

        # graphical parameters:
        self.vertical = False

        self.all_ages = pd.read_table(ages_file, ) #, names=['name','age','type'])
        ### TODO: if at least one of these is missing.
        if not set(('taxon', 'genetree')) & set(self.all_ages.columns):
            self.all_ages = splitname2taxongenetree(self.all_ages, "name")

        self.dup_ages = self.all_ages[self.all_ages.type == 'dup'].copy()

        if no_edited:
            self.load_edited_set(no_edited)
            self.not_edited = self.dup_ages.name.apply(lambda x: x not in self.edited_set)
            self.dup_ages = self.dup_ages[self.not_edited]
        
        print('shape:', self.dup_ages.shape)
        self.dup_ages.drop_duplicates(inplace=True)
        print('shape after drop_dup:', self.dup_ages.shape)
        self.dup_ages.reset_index(drop=True, inplace=True)
        self.taxa_ages = self.dup_ages.groupby(['taxon'], sort=False)
        self.taxa = sorted(self.taxa_ages.groups.keys())
        #"""A dot is used to separate words (genre.species)"""
        
        #newdata = [taxa_ages.get_group(lab)[age_keys] for lab in labels]

    def load_phyltree(self, phyltreefile=None, ensembl_version=None):
        phyltreefile = phyltreefile if phyltreefile else self.phyltreefile
        ensembl_version = ensembl_version if ensembl_version else \
                          self.ensembl_version
        self.phyltree = PhylTree.PhylogeneticTree(phyltreefile.format(
                                                            ensembl_version))

    ### Plotting utility methods ###
    def save_or_show(self, outfile=None):
        """show data unless `outfile` is given"""
        if outfile:
            self.fig.savefig(outfile, bbox_inches='tight')
        elif not plt.isinteractive():
            #self.fig.show()
            plt.show()
            #plt.close(self.fig)


    def colorize_taxa(self, alpha=0.5, cmap_name='Dark2'):
        """Create colormap for taxa (add a column to the data + create a dict)"""
        cmap = plt.get_cmap(cmap_name, len(self.taxa))
        self.taxonalpha = alpha
        self.taxon2color = {taxon: cmap(i) for i, taxon in enumerate(self.taxa)}
        self.dup_ages['taxoncolor'] = self.dup_ages.taxon.apply(self.taxon2color.get)


    def scatter(self, x, y, xlim=None, ylim=None):
        """scatter plot of x~y, colorized by taxon."""
        #self.ax = self.dup_ages.plot.scatter(x, y, c=self.dup_ages.taxoncolor,
        #                                     alpha=self.taxonalpha)
        #self.fig = self.ax.figure
        self.fig, self.ax = plt.subplots()
        for taxon in self.taxa:
            data = self.taxa_ages.get_group(taxon)
            #print("Taxon: %s\n" % taxon, data.head())
            try:
                data.plot.scatter(x, y, ax=self.ax,
                                  color=self.taxon2color[taxon],
                                  alpha=self.taxonalpha, label=taxon)
            except KeyError as err:
                print(('Data column %r not available for scatter plot. Check '
                       'if your data can be converted to float'), file=sys.stderr)
                raise
        if xlim:
            self.ax.set_xlim(*(float(x) for x in xlim.split(',')))
        if ylim:
            self.ax.set_ylim(*(float(y) for y in ylim.split(',')))

        self.ax.legend()
        #self.save_or_show(outfile)
        return self.ax


    def violin(self, x, y, xlim=None, ylim=None):
        # Here, otherwise it changes my matplotlib rcParams
        from seaborn import violinplot
        self.fig, self.ax = plt.subplots()
        #for taxon in self.taxa:
        #    data = self.taxa_ages.get_group(taxon)
        #if ylim:
        #    ylim = [float(yl) for yl in ylim.split(',')]
        #    data = self.dup_ages[self.dup_ages[y] < ylim[1]]
        #else:
        data = self.dup_ages
        violinplot(x, y, data=data, scale='width', cut=0, ax=self.ax)
        plt.setp(self.ax.xaxis.get_majorticklabels(), rotation=45, ha='right',
                 va='top')
        if ylim:
            ylim = [float(yl) for yl in ylim.split(',')]
            self.ax.set_ylim(ylim)

        return self.ax


    ### Histogram plotting methods ###
    def make_hist_data(self, taxa=None):
        """Given the selected taxa, return the appropriate data, colors, legend"""
        label_len = max(len(lab) for lab in self.taxa)
        label_fmt = "%%-%ds" % label_len

        taxa = taxa if taxa else self.taxa

        data = [self.taxa_ages.get_group(lab)[self.age_key].dropna() for lab in taxa]
        colors      = [self.taxon2color[lab] for lab in taxa]
        labs_legend = [label_fmt % lab       for lab in taxa]
        return data, colors, labs_legend


    def lineage_hist(self, lineage=None, nbins=None):
        """lineage: species name"""
        nbins = nbins if nbins else self.default_nbins
        
        taxa = set(self.taxa)

        if lineage:
            full_lineage = set(self.phyltree.dicLinks[self.phyltree.root][lineage])
            taxa &= full_lineage 

        # set up hist_coords, in case `add_edited_prop` is used.
        self.hist_coords = {tax: (None, 0) for tax in taxa} # age not needed.

        data, colors, labels = self.make_hist_data(taxa)

        print("plotting histogram")
        
        self.fig, self.ax = plt.subplots()
        self.ax.hist(data, bins=nbins,
                     histtype='barstacked',
                     rwidth=1,
                     color=colors,
                     edgecolor='none',
                     label=labels)
        self.ax.legend(fontsize='x-small')
        self.ax.set_title("Distribution of the age of duplications")
        self.ax.set_xlabel(self.age_key)
        self.ax.set_ylabel("Number of duplications")
        #plt.gca().set_xlim(0,2)
        return self.ax


    def walk_phylsubtree(self):
        """Return an iterator that will progress on the tree from leaves to root.
        
        It also rotates forks such that the displayed tree will show branches
        that are most distant from the root (in number of nodes) on one side.

        Taxa returned by the iterator are space-separated
        """
        print("Loading species tree")
        root, subtree = self.phyltree.getSubTree(self.taxa)

        # reorder branches in a visually nice manner:
        ladderize(subtree, root)
        #print(" ---\n subtree: ")
        #for anc, children in subtree.items():
        #    print(label_fmt % anc, children)
        get_children = lambda tree, node: tree.get(node, [])
        dfw = dfw_descendants_generalized(subtree, get_children, queue=[root])
                                          #queue=get_children(subtree, root))
        return reversed(list(dfw))


    def assign_subplots(self):
        """Give the best subplot number to plot each data, such that it mirrors
        the phylogenetic tree.
        """
        label_fmt="%-30s"

        # TODO: move to the __init__ function ?
        self.hist_coords = {}
        """'species': (x, y) i.e (age, subplot)"""

        self.subs_taxa = []
        """at each index: tuple of taxa appearing on the corresponding subplot"""

        self.treeforks = []
        """segment to be drawn on figure to visualize the tree: (x, y0, y1)"""

        print(" ---\n Assigning ancestors to subplots:")
        for anc1, anc2list in self.walk_phylsubtree():
            print(label_fmt % anc1, anc2list)
            #dotanc1 = anc1.replace(' ', '.')
            subs = []
            for anc2 in anc2list:
                #dotanc2 = anc2.replace(' ', '.')
                _, sub = self.hist_coords.get(anc2, (None, None))
                if not sub:
                    sub = len(self.subs_taxa)
                    self.subs_taxa.append(set((anc2,)))
                    self.hist_coords[anc2] = (self.phyltree.ages[anc2], sub)
                subs.append(sub)
            age_anc1 = self.phyltree.ages[anc1]
            low_sub = max(subs)
            self.treeforks.append((age_anc1, low_sub, min(subs)))
            self.hist_coords[anc1] = (age_anc1, low_sub)
            self.subs_taxa[low_sub].add(anc1)

        # remove root from subs_taxa:
        self.subs_taxa[low_sub].remove(anc1)

        print(" ---\n Assigned coordinates (age, subplot)")
        for anc, coords in self.hist_coords.items():
            print(label_fmt % anc, ":", "%3d, %2d" % coords)
        for i, subdatalabels in enumerate(self.subs_taxa):
            print("%2d" %i, ', '.join(subdatalabels))

    # Add mouse events: print selected data
    def onpick_bar_printdata(self, pickevent):
        """Print names of data contained in a histogram bar when mouse clicked.
        """
        self.pickevent = pickevent
        self.picked_bar = pickevent.artist
        #bars = picked_bars.axes.containers
        #picked_taxa = [ltxt.get_text().rstrip().replace(' ', '.')
        #               for ltxt in picked_bar.axes.legend_.get_texts()]

        # get label (taxon)
        picked_taxa = self.subs_taxa[self.axes.tolist().index(self.picked_bar.axes)]

        print("bar x0:", self.picked_bar.get_x())
        print("bar x1:", self.picked_bar.get_x() + self.picked_bar.get_width())
        pick_age = pickevent.mouseevent.ydata if self.vertical else \
                    pickevent.mouseevent.xdata
        print("pick_age:", pick_age)
        axes_data = pd.concat((self.taxa_ages.get_group(tax) for tax in picked_taxa))
        axes_data = axes_data.dropna(subset=[self.age_key])
        #axes_data, _, _ = self.make_hist_data(taxa=picked_taxa)
        #axes_data = pd.concat(axes_data)
        
        print(self.data_bins.keys())
        picked_bins = self.data_bins[picked_taxa.pop()]
        picked_bin = (pick_age > picked_bins).sum()
        assert picked_bin < len(picked_bins) + 1
        print("picked_bins:", picked_bins)
        print("picked_bin:", picked_bin)
        binned_data = bin_data(axes_data, picked_bins, binvar=self.age_key)
        picked_data = binned_data.get_group(picked_bin)[['name', self.age_key]]
        #max_len = picked_data.name.apply(lambda x: len(x)).max()
        #fmt_name = lambda name: ' ' * (max_len - len(name)) + name
        print(picked_data.to_string(index=False)) #, formatters=(fmt_name, None)))
        print("executed onpick action")


    def tree_hist(self, nbins=None, vertical=False):
        # labels, newdata, self.taxa_ages, hist_coords, subs_labels, treeforks, 
        """Draw histograms on top of a given tree structure.
        
        It needs a list of data, and the list of tree forks positions.

        Arguments
        ---------
        labels:      list of data labels
        newdata:     list of data (each being a 1-dimensional array of numbers)
        hist_coords: positions of each histogram   (see `assign_subplots`)
        subs_labels: subplot number of each hist.  (see `assign_subplots`)
        treeforks  : coordinates of each tree fork (see `assign_subplots`)
        outfile
        vertical   : True/False, whether to draw tree from bottom to top.
        """
        
        nbins = nbins if nbins else self.default_nbins

        self.data_bins = {}
        """Store bins associated with each data"""

        self.vertical = vertical

        n_subs = len(self.subs_taxa)
        figsize = (15, 11)

        if vertical:
            nrows, ncols = (1, n_subs)
            sharex, sharey = False, True
            bar_orientation = 'horizontal'
            text_rotation = 'horizontal'
            get_hist_coords = lambda lab: (self.hist_coords[lab][1],
                                           self.hist_coords[lab][0])
            invisible_spines = ('top', 'left', 'right')
            label_ageaxis = lambda axes, axislabel: axes[0].set_ylabel(axislabel)
        else:
            nrows, ncols = (n_subs, 1)
            sharex, sharey = True, False
            bar_orientation = 'vertical'
            text_rotation = 'vertical'
            get_hist_coords = lambda lab: self.hist_coords[lab]
            invisible_spines = ('top', 'bottom', 'right')
            label_ageaxis = lambda axes, axislabel: axes[-1].set_xlabel(axislabel)

        # oldest_age needed to scale all subplots identically (range=(0, oldest_age))
        # warning: return the index **name**, not the row number.
        age_argmax = self.dup_ages[self.age_key].argmax()
        oldest_lab, oldest_age = self.dup_ages[['taxon', self.age_key]].loc[age_argmax]
        print("Oldest: %s (%s)" % (oldest_age, oldest_lab))

        fig, axes = plt.subplots(nrows, ncols, sharex=sharex, sharey=sharey,
                                 figsize=figsize)
        transfig = fig.transFigure.inverted()
        cmap = plt.cm.get_cmap("Dark2", len(self.taxa))

        for ax_pos, ax in enumerate(axes):
            labs        = self.subs_taxa[ax_pos]
            data, colors, labs_legend = self.make_hist_data(labs)
            print("N observations: %s" % ([d.shape for d in data],), file=sys.stderr)
            print("nbins: %r" % nbins, file=sys.stderr)
            #print("orientation: %r" % bar_orientation, file=sys.stderr)
            #print("colors: %s" % colors, file=sys.stderr)
            print("label: %s" % labs_legend, file=sys.stderr)
            _, bins, _ = ax.hist(data, bins=nbins, histtype='barstacked', rwidth=1,
                                 range=(0, oldest_age),
                                 orientation=bar_orientation,
                                 color=colors,
                                 edgecolor='none',
                                 label=labs_legend,
                                 picker=True) # allow mouse selecting bars.
            # Label tree forks.
            for lab in labs:
                self.data_bins[lab] = bins
                x, y = get_hist_coords(lab)
                ax.text(x, y, lab, rotation=text_rotation, va='bottom', ha='left', 
                        fontsize='x-small')
            #set_title(ax, "plot %s" % ax_pos)
            ax.legend(prop={'size': 'xx-small', 'family': 'monospace'})
            #ax.axis('off')
            for spine in invisible_spines:
                ax.spines[spine].set_visible(False)
            #invert_axis(ax)
            if vertical: ax.invert_yaxis()
        label_ageaxis(axes, self.age_key + ' (My)')

        # draw tree branches
        print("n_subs =", n_subs)
        for anc1_age, sub1, sub2 in self.treeforks:
            #print(label_fmt % anc1, "%3d %2d %2d" % (anc1_age, sub1, sub2))
            ax1 = axes[sub1]
            ax2 = axes[sub2]
            if vertical:
                coord1 = transfig.transform(ax1.transData.transform([0, anc1_age]))
                coord2 = transfig.transform(ax2.transData.transform([0, anc1_age]))
            else:
                coord1 = transfig.transform(ax1.transData.transform([anc1_age, 0]))
                coord2 = transfig.transform(ax2.transData.transform([anc1_age, 0]))
            line = mpl.lines.Line2D((coord1[0], coord2[0]), (coord1[1], coord2[1]),
                                   transform=fig.transFigure,
                                   color=grey10, alpha=0.4, linestyle='dashed')
            #ax1.text(anc1_age, 0, anc1)
            #ax1.plot(anc1_age, 0, 'ro', markersize=20)
            #ax2.plot(anc1_age, 0, 'ro', markersize=20)
            fig.lines.append(line)

        fig.canvas.mpl_connect('pick_event', self.onpick_bar_printdata)

        self.fig = fig
        self.axes = axes
        return self.axes


    def score_bins(self, score_name='score', scoring=None,
                    *scoring_args):
        """Give a score to each bin (from histogram) and plot it.
        
        A custom scoring function must be given, that will take
        a dataframe/dict in input."""
        
        self.scored_data = {}
        
        for lab, data in self.taxa_ages:
            data = data.dropna(subset=[self.age_key])
            #print("Nb of NA values:", data[self.age_key].isnull().sum())
            #print("min:", data[self.age_key].min())
            #print("max:", data[self.age_key].max())
            bins = self.data_bins[lab]
            binned_groups = bin_data(data, bins, binvar=self.age_key,
                                     outvar='name')
            #print("bin indices found:", binned_groups.groups.keys())
            scored_serie = binned_groups.aggregate(scoring, *scoring_args)
            scored_serie.name = score_name
            bin_serie = pd.Series([(bins[x-1] + bins[x])/2 for x in scored_serie.index],
                                  index=scored_serie.index, name='bin_x')
            scored_taxon = pd.concat([scored_serie, bin_serie], axis=1)
            self.scored_data[lab] = scored_taxon
            yield lab, scored_taxon


    def add_edited_prop(self, treeforest):
        """Add a point for the proportion of edited nodes in each histogram bar"""
        score_name = 'prop_edited'
        self.load_edited_set(treeforest)

        for taxon, data in self.score_bins(score_name, intersect_size,
                                           self.edited_set):
            ax = self.axes[self.hist_coords[taxon][1]]
            ax_2 = ax.twinx()
            ax_2.set_ylim(0,1) # TODO: what if score not in [0,1] ?
            ax_2.plot(data.bin_x, data, 'b_', alpha=0.7)
            ax_2.spines['bottom'].set_visible(False)
            ax_2.spines['right'].set_visible(True)
            #ax_2.spines['left'].set_visible(False)
            #ax_2.spines['top'].set_visible(False)


### Script commands ###
#COMMANDS = ['lineage', 'tree', 'scatter']
def run(command, ages_file, phyltreefile=None, ensembl_version=None,
        outfile=None, lineage=None,
        show_edited=None, no_edited=False, age_key=DEFAULT_AGE_KEY,
        nbins=DEFAULT_NBINS, vertical=False, x=None, y=None, xlim=None, ylim=None):

    if outfile:
        plt.switch_backend("Agg")

    dv = DataVisualizor(ages_file, no_edited=no_edited, age_key=age_key)
    
    if command in ('tree', 'lineage'):
        dv.colorize_taxa(alpha=1)
        dv.load_phyltree(phyltreefile, ensembl_version)
        if command == 'lineage':
            dv.lineage_hist(lineage, nbins)
        elif command == 'tree':
            dv.assign_subplots()
            dv.tree_hist(nbins, vertical)
            print("finished drawing tree & hist.")
        if show_edited:
            dv.add_edited_prop(show_edited)
    elif command == 'scatter':
        dv.colorize_taxa(alpha=0.5)
        dv.scatter(x, y, xlim, ylim)
    elif command == 'violin':
        dv.violin(x, y, xlim, ylim)

    dv.save_or_show(outfile)


# TODO: delete
def cmd_lineage(command, ages_file, outfile=None, lineage=None, show_edited=None,
        no_edited=False, age_key=DEFAULT_AGE_KEY, nbins=DEFAULT_NBINS):
    """histogram of duplications along a selected lineage"""

def cmd_tree(ages_file, outfile=None, nbins=DEFAULT_NBINS, vertical=False,
             show_edited=None, no_edited=False):
    """histogram on top of phylogenetic tree"""
#
def cmd_scatter(ages_file, x, y, outfile=None, show_edited=None, no_edited=False):
    """scatter plot"""

def cmd_violin(ages_file, x, y, outfile=None, show_edited=None, no_edited=False):
    """violin plot"""


def document_commands():
    """Add description of each script command to the module docstring."""
    global __doc__
    global_objects = globals()
    
    cmd_fmt = '  - {:%d}: {}' % max(len(cmd) for cmd_name in COMMANDS)
    
    __doc__ += '\n\nCOMMANDS:\n'
    for cmd_name in COMMANDS:
        cmd_func = global_objects['cmd_' + cmd_name]
        __doc__ += cmd_fmt.format(cmd_name, cmd_func.__doc__) + '\n'


CMD_FUNC = {cmd_name: globals()['cmd_' + cmd_name] for cmd_name in COMMANDS}


if __name__=='__main__':
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    ### Generic arguments
    parent_parser = argparse.ArgumentParser(add_help=False)
    parent_parser.add_argument('ages_file')
    parent_parser.add_argument('outfile', nargs='?')
    parent_parser.add_argument('-a', '--age-key', default=DEFAULT_AGE_KEY,
                               help='[%(default)s]')
    parent_parser.add_argument('-b', '--nbins', type=int, default=DEFAULT_NBINS,
                               help='[%(default)s]')
    
    process_edited_parser = parent_parser.add_mutually_exclusive_group()
    # these two options must be given the treeforest file.
    process_edited_parser.add_argument('--show-edited',
                                       metavar='TREEFOREST',
                                       help='plot the proportion of edited nodes')
    process_edited_parser.add_argument('--no-edited',
                                       metavar='TREEFOREST',
                                       help='delete edited gene trees from data')
    
    ### Command-specific arguments
    subparsers = parser.add_subparsers(dest='command')

    for cmd_name in COMMANDS:
        cmd_func = CMD_FUNC[cmd_name] # only used for the __doc__ ...
        cmd_parser = subparsers.add_parser(cmd_name,
                                           description=cmd_func.__doc__,
                                           parents=[parent_parser],
                        formatter_class=argparse.RawDescriptionHelpFormatter)
        for cmd_arg in CMD_ARGS[cmd_name]:
            #print(cmd_name)
            #print(cmd_arg)
            kwargs = cmd_arg[1] if len(cmd_arg) > 1 else {}
            cmd_parser.add_argument(*cmd_arg[0], **kwargs)

    args = parser.parse_args()
    
    dictargs = vars(args)

    #print(dictargs)

    run(**dictargs)
    #command = CMD_FUNC[dictargs.pop('command')]
    #command(**dictargs)
