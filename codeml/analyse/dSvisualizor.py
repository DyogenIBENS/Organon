#!/usr/bin/env python3

# Evolution of `plot_dS.py`. Refactor code using classes.
# Use Pandas features more systematically.

"""Visualise values from given table (reconstructed ages) per taxonomic branch."""
"""
USAGE:
    ./dSvisualizor.py command <ages_file> <outfile>
"""

PHYLTREEFILE = "~/GENOMICUS{0}/PhylTree.Ensembl.{0}.conf"
ENSEMBL_VERSION = 85
DEFAULT_NBINS   = 100
DEFAULT_AGE_KEY = 'age_dS'


COMMANDS = ['lineage', 'tree', 'scatter', 'violin']

CMD_ARGS = {
        'lineage': [
                    (('-a', '--age-key'), {'default': DEFAULT_AGE_KEY}),
                    (('-l', '--lineage'),),
                    (('-p', '--phyltreefile'),   {'default': PHYLTREEFILE}),
                    (('-e', '--ensembl-version'),{'default': ENSEMBL_VERSION}),
                    (('-x', '--xlim'),)],
        'tree':    [(('-a', '--age-key'), {'default': DEFAULT_AGE_KEY}),
                    (('-V', '--vertical'),       {'action':'store_true'}),
                    (('-p', '--phyltreefile'),   {'default': PHYLTREEFILE}),
                    (('-e', '--ensembl-version'),{'default': ENSEMBL_VERSION}),
                    (('-x', '--xlim'),),
                    (('-y', '--ylim'),),
                    (('-t', '--title'),),
                    (('--sharescale',), {'action':'store_true'})],
        'scatter': [(('-x',),), (('-y',),), (('--xlim',),), (('--ylim',),)],
        'violin':  [(('-x',),),
                    (('-y',), {'default': DEFAULT_AGE_KEY}),
                    (('--xlim',),), (('--ylim',),)],}

import os.path
import re
import bz2
import pickle
try:
    import argparse_custom as argparse
except ImportError:
    import argparse
from textwrap import dedent

import matplotlib as mpl
#print(mpl.get_backend())
mpl.use('Agg', warn=False) # for figures to show up when the script is called from the shell
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator #AutoLocator
#plt.ion()
import numpy as np
import pandas as pd
import LibsDyogen.myPhylTree as PhylTree

from dendron.climber import dfw_descendants_generalized
from dendron.sorter import ladderize

import logging
logger = logging.getLogger(__name__)
ch = logging.StreamHandler()
ch.setFormatter(logging.Formatter("%(levelname)s:l.%(lineno)d%:(funcName)s:%(message)s"))
logger.addHandler(ch)


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

RE_TAXON = re.compile(r'[A-Z][A-Za-z_.-]+(?=ENSGT)')
#PAT_TAXON = r'^(ENS[A-Z]+G|[A-Z][A-Za-z_.-]+)(ENSGT[0-9]+|)(.*)$'
#PAT_TAXON = r'^([A-Z][A-Za-z_.-]+)(ENSGT[0-9]+)(.*)$'
PAT_TAXON = r'^(ENS[A-Z]+G|[A-Z][A-Za-z_.-]+(?=ENSGT))(|ENSGT[0-9]+)([0-9]+|[.A-Za-z`]*)$'
RE_FILTER = re.compile(r'(\S+)\s*(>|>=|==|<=|<|!=)\s*(\S+)')
CMP_METHODS = {'>':  '__gt__',
               '>=': '__ge__',
               '==': '__eq__',
               '<=': '__le__',
               '<':  '__lt__',
               '!=': '__ne__',
               'in': 'isin'}
RE_QUOTED = re.compile(r'^[\'"](.*)[\'"]$')
RE_TUPLED = re.compile(r'^\((.*)\)$')  ### TODO


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
        logger.info('Loading the set of edited nodes:')
        pickled_file = os.path.basename(os.path.splitext(treeforest)[0]) + \
                        '.editedset.pickle'
        if os.path.exists(pickled_file):
            print('from pickle...')
            with open(pickled_file, 'rb') as pickle_in:
                self.edited_set = pickle.load(pickle_in)
        else:
            logger.info('Searching tree...')
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

    def __init__(self, ages_file, no_edited=None, age_key=DEFAULT_AGE_KEY,
                 filter=None):
        """Load and format data"""

        self.edited_set = None
        self.not_edited = None
        self.no_edited  = no_edited
        self.age_key = age_key
        self.phyltree = None

        # Parse the filter command:
        if filter is not None:
            try:
                col, comp, val = RE_FILTER.match(filter).groups()
            except AttributeError:
                raise

            quoted_match = RE_QUOTED.match(val)
            val = quoted_match.group(1) if quoted_match else float(val)

        # graphical parameters:
        self.vertical = False

        self.all_ages = pd.read_table(ages_file, sep='\t') #, names=['name','age','type'])

        # Check the integrity of the file (dtype object means some rows are wrong,
        # for example there is a header row in the middle of the document)
        assert (self.all_ages.dtypes != np.dtype(np.object)).any(), \
                "Data types of columns not understood, check their integrity. %s" % \
                    self.all_ages.dtypes
        #print(ages_file)
        #print(self.all_ages)
        ### TODO: if at least one of these is missing.
        if not set(('taxon', 'genetree')) & set(self.all_ages.columns):
            self.all_ages = splitname2taxongenetree(self.all_ages, "name")

        filtered = True if filter is None else getattr(self.all_ages[col], CMP_METHODS[comp])(val)
        self.ages = self.all_ages[filtered & (self.all_ages['calibrated'] == 0)].copy()
        logger.debug('shape: %s', self.ages.shape)

        if no_edited:
            self.load_edited_set(no_edited)
            self.not_edited = self.ages.name.apply(lambda x: x not in self.edited_set)
            self.ages = self.ages[self.not_edited]
        
        self.ages.drop_duplicates(inplace=True)
        assert self.ages.shape[0] > 0, "All data was filtered out."
        logger.debug('shape after drop_dup: %s', self.ages.shape)
        self.ages.reset_index(drop=True, inplace=True)
        self.taxa_ages = self.ages.groupby(['taxon'], sort=False)
        self.taxa_evt_ages = self.ages.groupby(['taxon', 'type'], sort=False)
        self.taxa_evt = sorted(k for k in self.taxa_evt_ages.groups.keys() if k[1] != 'leaf')
        self.taxa = sorted(set(taxon for taxon, evt in self.taxa_evt))
        self.all_taxa = self.all_ages.taxon.unique() # useful when considering extra groups
        #print('taxa:', ', '.join(self.taxa))

        #"""A dot is used to separate words (genre.species)"""
        
        #newdata = [taxa_ages.get_group(lab)[age_keys] for lab in labels]

    def load_phyltree(self, phyltreefile=None, ensembl_version=None):
        phyltreefile = phyltreefile if phyltreefile else self.phyltreefile
        ensembl_version = ensembl_version if ensembl_version else \
                          self.ensembl_version
        self.phyltree = PhylTree.PhylogeneticTree(
                            os.path.expanduser(
                                phyltreefile.format(ensembl_version)))

    ### Plotting utility methods ###
    def save_or_show(self, outfile=None):
        """show data unless `outfile` is given"""
        if outfile:
            self.fig.savefig(outfile, bbox_inches='tight')
        elif not plt.isinteractive():
            #self.fig.show()
            plt.show()
            #plt.close(self.fig)


    def colorize_taxa(self, alpha=0.5, cmap_name='tab20b'): #'Dark2'
        """Create colormap for taxa (add a column to the data + create a dict)"""
        cmaps = {'dup': plt.get_cmap(cmap_name, len(self.taxa)),
                 'spe': plt.get_cmap('tab20c',  len(self.taxa))}
        self.taxonalpha = alpha
        #self.taxon2color = {taxon: cmap(i) for i, taxon in enumerate(self.taxa)}
        #self.ages['taxoncolor'] = self.ages.taxon.apply(self.taxon2color.get)
        self.taxon_evt_2color = {(taxon, evt): cmaps[evt](i) \
                                    for i, taxon in enumerate(self.taxa) \
                                    for evt in ('dup', 'spe')}
        #use `hatch=` in plt.bar for patterning.

        self.ages['taxon_evt_color'] = self.ages[['taxon', 'type']].apply(
                lambda k: self.taxon_evt_2color[tuple(k)], axis=1, raw=True)


    def scatter(self, x, y, xlim=None, ylim=None):
        """Scatter plot of y~x, colorized by taxon."""
        #self.ax = self.ages.plot.scatter(x, y, c=self.ages.taxoncolor,
        #                                     alpha=self.taxonalpha)
        #self.fig = self.ax.figure
        self.fig, self.ax = plt.subplots()
        for taxon, evt in self.taxa_evt:
            #if evt == 'leaf': continue

            data = self.taxa_evt_ages.get_group((taxa, evt))
            datacolor = self.taxon_evt_2color[(taxon, evt)]
            #print("Taxon: %s\n" % taxon, data.head())
            try:
                data.plot.scatter(x, y, ax=self.ax,
                                  color=datacolor,
                                  alpha=self.taxonalpha, label=taxon)
            except KeyError as err:
                err.args += ('Data column %r not available for scatter plot. Check '
                             'if your data can be converted to float',)
                raise
        if xlim is not None:
            self.ax.set_xlim(xlim)
        if ylim is not None:
            self.ax.set_ylim(ylim)

        self.ax.legend()
        #self.save_or_show(outfile)
        return self.ax


    def violin(self, x, y, xlim=None, ylim=None):
        """Violin plot of y~x"""
        # Here, otherwise it changes my matplotlib rcParams
        from seaborn import violinplot
        self.fig, self.ax = plt.subplots()
        #for taxon in self.taxa:
        #    data = self.taxa_ages.get_group(taxon)
        #if ylim:
        #    ylim = [float(yl) for yl in ylim.split(',')]
        #    data = self.ages[self.ages[y] < ylim[1]]
        #else:
        data = self.ages
        violinplot(x, y, data=data, scale='width', cut=0, ax=self.ax)
        plt.setp(self.ax.xaxis.get_majorticklabels(), rotation=45, ha='right',
                 va='top')
        if ylim is not None:
            self.ax.set_ylim(ylim)

        return self.ax


    ### Histogram plotting methods ###
    def make_hist_data(self, taxa=None):
        """Given the selected taxa, return the appropriate data, colors, legend.
        This data is given to matplotlib to plot a single histogram."""
        
        label_len = max(len(lab) for lab in self.taxa)
        label_fmt = "%%-%ds (%%s)" % label_len

        taxa = taxa if taxa is not None else self.taxa

        #duptaxa = [k for k in self.taxa if k in self.taxa_ages.groups]
        data = [self.taxa_evt_ages.get_group((lab, evt))[self.age_key].dropna().values \
                    for lab, evt in self.taxa_evt \
                    if lab in taxa]
                    #if lab in self.taxa_ages.groups else [] \
        #data = [self.taxa_ages.groups.get(lab, {self.age_key: []})[self.age_key].dropna() \
        #           for lab in taxa
        colors      = [self.taxon_evt_2color[key] for key in self.taxa_evt \
                                                            if key[0] in taxa]
        labs_legend = [label_fmt % (lab, evt) for lab, evt in self.taxa_evt \
                                        if lab in taxa]
        #Quick fix when no rows in data (meaning, only NaN values)
        #Should dropna upstream!!!
        if not all(len(d) for d in data):
            data, colors, labs_legend = zip(*[(d,c,l) for d,c,l in \
                                                zip(data,colors,labs_legend) \
                                                if len(d)])
            logger.warning('Only NaN for some taxa')
        return data, colors, labs_legend


    def lineage_hist(self, lineage=None, nbins=None):
        """
        Histogram of duplications along a selected lineage.
        
        Lineage: a species name"""
        nbins = nbins if nbins else self.default_nbins
        
        taxa = set(self.taxa)

        if lineage:
            full_lineage = set(self.phyltree.dicLinks[self.phyltree.root][lineage])
            taxa &= full_lineage 

        # set up hist_coords, in case `add_edited_prop` is used.
        self.hist_coords = {tax: (None, 0) for tax in taxa} # age not needed.

        data, colors, labels = self.make_hist_data(taxa)

        logger.info("plotting histogram")
        
        self.fig, self.ax = plt.subplots()
        self.ax.hist(data, bins=nbins,
                     histtype='barstacked',
                     rwidth=1,
                     color=colors,
                     alpha=self.taxonalpha,
                     edgecolor='none',
                     label=labels)
        self.ax.legend(fontsize='x-small')
        self.ax.set_title("Distribution of %s" % self.age_key)
        self.ax.set_xlabel(self.age_key)
        self.ax.set_ylabel("Number of nodes")
        #plt.gca().set_xlim(0,2)
        return self.ax


    def walk_phylsubtree(self):
        """Return an iterator that will progress on the tree from leaves to root.
        
        It also rotates forks such that the displayed tree will show branches
        that are most distant from the root (in number of nodes) on one side.

        Taxa returned by the iterator are space-separated
        """
        logger.info("Loading species tree")
        #logger.debug(self.taxa)
        root, subtree = self.phyltree.getSubAncTree(self.taxa)

        # reorder branches in a visually nice manner:
        ladderize(subtree, root)
        #logger.debug(" ---\n subtree: ")
        #for anc, children in subtree.items():
        #    logger.debug(label_fmt % anc, children)
        get_children = lambda tree, node: tree.get(node, [])
        dfw = dfw_descendants_generalized(subtree, get_children, queue=[root])#,
                                          #include_leaves=True)
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

        logger.debug(" ---\n Assigning ancestors to subplots:")
        for anc1, anc2list in self.walk_phylsubtree():
            logger.debug((label_fmt % anc1) + str(anc2list))
            #dotanc1 = anc1.replace(' ', '.')
            subs = []
            for anc2 in anc2list:
                #dotanc2 = anc2.replace(' ', '.')
                _, sub = self.hist_coords.get(anc2, (None, None))
                if sub is None:
                    sub = len(self.subs_taxa)
                    self.subs_taxa.append(set((anc2,)))
                    self.hist_coords[anc2] = (self.phyltree.ages[anc2], sub)
                subs.append(sub)
            age_anc1 = self.phyltree.ages[anc1]
            low_sub = max(subs)
            self.treeforks.append((age_anc1, low_sub, min(subs)))
            self.hist_coords[anc1] = (age_anc1, low_sub)
            self.subs_taxa[low_sub].add(anc1)

        assert len(self.subs_taxa), \
                "No subplots for taxa [%s]. Use the 'lineage' command." % ", ".join(self.taxa)
        # remove root from subs_taxa:
        #self.subs_taxa[low_sub].remove(anc1)

        logger.debug(" ---\n Assigned coordinates (age, subplot):")
        for anc, coords in self.hist_coords.items():
            logger.debug((label_fmt % anc) + " : %3d, %2d" % coords)
        for i, subdatalabels in enumerate(self.subs_taxa):
            logger.debug("%2d %s", i, ', '.join(subdatalabels))

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

        logger.debug("bar x0:", self.picked_bar.get_x())
        logger.debug("bar x1:", self.picked_bar.get_x() + self.picked_bar.get_width())
        pick_age = pickevent.mouseevent.ydata if self.vertical else \
                    pickevent.mouseevent.xdata
        logger.info("pick_age:", pick_age)
        axes_data = pd.concat((self.taxa_ages.get_group(tax) for tax in picked_taxa))
        axes_data = axes_data.dropna(subset=[self.age_key])
        #axes_data, _, _ = self.make_hist_data(taxa=picked_taxa)
        #axes_data = pd.concat(axes_data)
        
        logger.debug(self.data_bins.keys())
        picked_bins = self.data_bins[picked_taxa.pop()]
        picked_bin = (pick_age > picked_bins).sum()
        assert picked_bin < len(picked_bins) + 1
        logger.debug("picked_bins:", picked_bins)
        logger.debug("picked_bin:", picked_bin)
        binned_data = bin_data(axes_data, picked_bins, binvar=self.age_key)
        picked_data = binned_data.get_group(picked_bin)[['name', self.age_key]]
        #max_len = picked_data.name.apply(lambda x: len(x)).max()
        #fmt_name = lambda name: ' ' * (max_len - len(name)) + name
        logger.debug(picked_data.to_string(index=False)) #, formatters=(fmt_name, None)))
        logger.debug("executed onpick action")


    def tree_hist(self, nbins=None, vertical=False, xlim=None, ylim=None, sharescale=False,
                  title=None):
        # labels, newdata, self.taxa_ages, hist_coords, subs_labels, treeforks, 
        """
        Draw histograms on top of a given phylogenetic tree (e.g species tree).
        
        It needs a list of data, and the list of tree forks positions.

        Arguments
        ---------
        nbins: total number of bins for a row.
        vertical   : True/False, whether to draw tree from bottom to top.
        sharescale : whether the duplication axis should have the same limits.

        Used attributes
        ---------------
        labels:      list of data labels
        newdata:     list of data (each being a 1-dimensional array of numbers)
        hist_coords: positions of each histogram   (see `assign_subplots`)
        subs_labels: subplot number of each hist.  (see `assign_subplots`)
        treeforks  : coordinates of each tree fork (see `assign_subplots`)
        outfile
        """
        
        nbins = nbins if nbins else self.default_nbins

        self.data_bins = {}
        """Store bins associated with each data"""

        self.vertical = vertical

        n_subs = len(self.subs_taxa)
        figsize = (15, 11)
        
        if ylim: sharescale = True

        ticklocator = lambda: MaxNLocator(integer=True, nbins='auto')
        if vertical:
            nrows, ncols = (1, n_subs)
            sharex, sharey = sharescale, True
            bar_orientation = 'horizontal'
            text_rotation = 'horizontal'
            get_hist_coords = lambda lab: (self.hist_coords[lab][1],
                                           self.hist_coords[lab][0])
            #invisible_spines = ('top', 'left', 'right')
            invisible_spines = ('top', 'right')
            label_ageaxis = lambda axes, axislabel: axes[0].set_ylabel(axislabel)
            fix_dupticks = lambda ax: ax.xaxis.set_major_locator(ticklocator())
        else:
            nrows, ncols = (n_subs, 1)
            sharex, sharey = True, sharescale
            bar_orientation = 'vertical'
            text_rotation = 'vertical'
            get_hist_coords = lambda lab: self.hist_coords[lab]
            #invisible_spines = ('top', 'bottom', 'right')
            invisible_spines = ('top', 'right')
            label_ageaxis = lambda axes, axislabel: axes[-1].set_xlabel(axislabel)
            fix_dupticks = lambda ax: ax.yaxis.set_major_locator(ticklocator())

        # oldest_age needed to scale all subplots identically (range=(0, oldest_age))
        # warning: return the index **name**, not the row number.
        #print(self.ages)
        age_argmax = self.ages[self.age_key].idxmax()
        oldest_lab, oldest_age = self.ages[['taxon', self.age_key]].loc[age_argmax]
        logger.info("Oldest: %s (%s)", oldest_age, oldest_lab)

        fig, axes = plt.subplots(nrows, ncols, sharex=sharex, sharey=sharey,
                                 squeeze=False, figsize=figsize)
        axes = axes.flatten()
        transfig = fig.transFigure.inverted()
        cmap = plt.cm.get_cmap("Dark2", len(self.taxa))

        for ax_pos, ax in enumerate(axes):
            labs        = self.subs_taxa[ax_pos]
            data, colors, labs_legend = self.make_hist_data(labs)
            logger.info("* Labels: %s; N observations: %s; nbins: %r",
                         labs_legend,
                         ([d.shape for d in data],),
                         nbins)
            #logger.debug("orientation: %r", bar_orientation)
            #logger.debug("colors: %s", colors)
            if not data:
                logger.warning('No data. Skip')
                continue
            try:
                if len(data)==1:
                    data, = data

                _, bins, _ = ax.hist(data, bins=nbins, histtype='barstacked', rwidth=1,
                                 range=(0, oldest_age),
                                 orientation=bar_orientation,
                                 color=colors,
                                 edgecolor='none',
                                 label=labs_legend,
                                 picker=True) # allow mouse selecting bars.
            except BaseException as err:
                err.args += tuple((type(d), repr(d)) for d in data)
                raise

            fix_dupticks(ax)
            # Label tree forks.
            for lab in labs:
                self.data_bins[lab] = bins
                x, y = get_hist_coords(lab)
                #logger.debug('Plotting label %r at (%s, %s)' % (lab, x, y))
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
        if xlim is not None:
            axes[-1].set_xlim(xlim)
        if ylim is not None:
            axes[-1].set_ylim(ylim)

        # draw tree branches
        logger.debug("n_subs =", n_subs)
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

        if title:
            fig.suptitle(title)
        fig.canvas.mpl_connect('pick_event', self.onpick_bar_printdata)

        self.fig = fig
        self.axes = axes
        return axes


    def score_bins(self, score_name='score', scoring=None,
                    *scoring_args):
        """Give a score to each bin (from histogram) and plot it.
        
        A custom scoring function must be given, that will take
        a dataframe/dict in input."""
        
        self.scored_data = {}
        
        for lab, data in self.taxa_ages:
            data = data.dropna(subset=[self.age_key])
            #logger.debug("Nb of NA values:", data[self.age_key].isnull().sum())
            #logger.debug("min:", data[self.age_key].min())
            #logger.debug("max:", data[self.age_key].max())
            bins = self.data_bins[lab]
            binned_groups = bin_data(data, bins, binvar=self.age_key,
                                     outvar='name')
            #logger.debug("bin indices found:", binned_groups.groups.keys())
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
        show_edited=None, no_edited=False, age_key=DEFAULT_AGE_KEY, filter=None,
        nbins=DEFAULT_NBINS, vertical=False, x=None, y=None, xlim=None, ylim=None,
        sharescale=False, title=None):

    curr_backend = mpl.get_backend()
    if not outfile and "inline" not in curr_backend and curr_backend != "nbagg":
        #plt.switch_backend("TkAgg")
        try:
            plt.switch_backend("Qt5Agg")
        except ImportError:
            try:
                plt.switch_backend("Qt4Agg")
            except ImportError:
                plt.switch_backend("TkAgg")

    dv = DataVisualizor(ages_file, no_edited=no_edited, age_key=age_key,
                        filter=filter)
    if xlim is not None:
        xlim = tuple(float(x) for x in xlim.split(','))
    if ylim is not None:
        ylim = tuple(float(y) for y in ylim.split(','))
    
    if command in ('tree', 'lineage'):
        dv.colorize_taxa(alpha=1)
        dv.load_phyltree(phyltreefile, ensembl_version)
        if command == 'lineage':
            dv.lineage_hist(lineage, nbins)
        elif command == 'tree':
            dv.assign_subplots()
            dv.tree_hist(nbins, vertical, xlim, ylim, sharescale, title)
            logger.info("Finished drawing tree & hist.")
        if show_edited:
            dv.add_edited_prop(show_edited)
    elif command == 'scatter':
        assert x is not None and y is not None
        dv.colorize_taxa(alpha=0.5)
        dv.scatter(x, y, xlim, ylim)
    elif command == 'violin':
        assert x is not None and y is not None
        dv.violin(x, y, xlim, ylim)

    dv.save_or_show(outfile)



#def document_commands():
#    """Add description of each script command to the module docstring."""
#    global __doc__
#    global_objects = globals()
#    
#    cmd_fmt = '  - {:%d}: {}' % max(len(cmd) for cmd_name in COMMANDS)
#    
#    __doc__ += '\n\nCOMMANDS:\n'
#    for cmd_name in COMMANDS:
#        cmd_func = global_objects['cmd_' + cmd_name]
#        __doc__ += cmd_fmt.format(cmd_name, cmd_func.__doc__) + '\n'


CMD_FUNC = {'lineage': DataVisualizor.lineage_hist,
            'tree': DataVisualizor.tree_hist,
            'scatter': DataVisualizor.scatter,
            'violin': DataVisualizor.violin}


if __name__=='__main__':

    #class RawDescArgDefFormatter(argparse.RawDescriptionHelpFormatter,
    #                             argparse.ArgumentDefaultsHelpFormatter):
    #    """Combine two argparse formatters: raw description + argument defaults"""
    #    pass

    parser = argparse.ArgumentParser(description=__doc__)

    ### Generic arguments
    parent_parser = argparse.ArgumentParser(add_help=False)
    parent_parser.add_argument('ages_file')
    parent_parser.add_argument('outfile', nargs='?')
    parent_parser.add_argument('-f', '--filter',
                               help="Filter rows on the value of a given "\
                                    "column.  e.g: 'type==\"spe\"' ")
    parent_parser.add_argument('-v', '--verbose', action='store_true')
    
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
                                           description=dedent(cmd_func.__doc__),
                                           parents=[parent_parser],
                        formatter_class=argparse.RawDescriptionHelpFormatter)#,
                                        #formatter_class=RawDescArgDefFormatter)
        for cmd_arg in CMD_ARGS[cmd_name]:
            kwargs = cmd_arg[1] if len(cmd_arg) > 1 else {}
            # Automatically append the default value in help
            default = kwargs.get('default')
            arghelp = kwargs.get('help', '')
            if default is not None and '%(default)s' not in arghelp:
                if arghelp: arghelp += ' '
                kwargs.update(help=(arghelp + '[%(default)s]'))
            # Add the argument
            cmd_parser.add_argument(*cmd_arg[0], **kwargs)

    args = parser.parse_args()
    
    dictargs = vars(args)
    
    loglevel = logging.INFO if dictargs.pop('verbose') else logging.WARNING
    logger.setLevel(loglevel)

    #print(dictargs)

    run(**dictargs)
    #command = CMD_FUNC[dictargs.pop('command')]
    #command(**dictargs)
