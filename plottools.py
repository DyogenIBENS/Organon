#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Code snippets for matplotlib/pandas plotting."""

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
import pandas as pd
from collections import namedtuple

import pandas as pd
import seaborn as sns

from dendron.climber import rev_dfw_descendants, iter_distleaves
import matplotlib.patches as patches
from matplotlib.path import Path
MOVETO, CURVE3, LINETO = Path.MOVETO, Path.CURVE3, Path.LINETO
### Plotting style in matplotlib ###

# Prefered_style: put in ~/.config/matplotlib/smoother
#mpl.style.use('smoother')

def softstyle():
    mpl.style.use('softer')
	## Change all black to dark grey
	#grey10 = '#1a1a1a'
	#grey45 = '#737373'
	#grey80 = '#CCCCCC'
	#mpl.rcParams['text.color'] = grey10
	#mpl.rcParams['axes.edgecolor'] = grey10
	#mpl.rcParams['axes.labelcolor'] = grey10
	#mpl.rcParams['axes.spines.top'] = False
	#mpl.rcParams['axes.spines.right'] = False
	#mpl.rcParams['xtick.color'] = grey10
	#mpl.rcParams['ytick.color'] = grey10
	#mpl.rcParams['grid.color'] = grey80
	#mpl.rcParams['boxplot.boxprops.color'] = grey10
	#mpl.rcParams['boxplot.capprops.color'] = grey10
	#mpl.rcParams['boxplot.flierprops.color'] = grey10
	#mpl.rcParams['boxplot.flierprops.markeredgecolor'] = grey10
	#mpl.rcParams['boxplot.whiskerprops.color'] = grey10
	#mpl.rcParams['hatch.color'] = grey10
    #mpl.rcParams['axes.prop_cycle'] = mpl.cycler('color', ['4C72B0', '55A868', 'C44E52', '8172B2', 'CCB974', '64B5CD'])
	#mpl.rcParams['grid.linestyle'] = ':'
	#mpl.rcParams['patch.edgecolor'] = grey10
	#mpl.rcParams['boxplot.flierprops.markeredgecolor'] = grey10
	#mpl.rcParams['boxplot.capprops.color'] = grey10
	#mpl.rcParams['legend.facecolor'] = grey45
	#mpl.rcParams['legend.framealpha'] = 0.2
	##mpl.rcParams['legend.edgecolor'] = grey10
	#mpl.rcParams['savefig.facecolor'] = 'none'
	##mpl.rcParams['savefig.frameon'] = False  #background frame transparent
	##mpl.rcParams['savefig.transparent'] = True # all background transparent
	#                                            # (including ggplot2 style)
	##mpl.style.use('ggplot')
    ##pd.set_option('display.max_colwidth', 85)


### Helping functions ###

def plotby(df, by, kind='line', ncols=3, sharex=True, sharey=True, **kwds):
    """Groupby and plot."""
    df_grouped = df.groupby(by)
    ngroups = len(df_grouped)
    nrows = ngroups // ncols + 1*(ngroups % ncols > 0)
    fig, axes = plt.subplots(nrows, ncols, sharex, sharey)
    #print(kde_axes, type(kde_axes))
    for i, (title, data) in enumerate(df_grouped):
        #print("title")
        ax = axes[i // ncols, i % ncols]
        subkwds = dict((key,v[i]) if isinstance(v, list)
                       else (key, v[title]) if isinstance(v, dict)
                       else(key,v)
                       for key,v in kwds.items())
        #print(subkwds)
        data.plot(kind=kind, ax=ax, title=title, **subkwds)

    # Hide unnecessary axes
    for j in range(i+1, nrows*ncols):
        ax = axes[j // ncols, j % ncols]
        ax.axis('off')

    return axes

# TODO:
def ordered_boxplot(df, x, y, order=None, **kwds):
    pass

# TODO:
def stackedbar(x, arr, ax=None, **kwds):
    assert isinstance(arr, np.ndarray)
    bar = plt.bar if ax is None else ax.bar
    stacked_bars = []
    bottom = None
    for row in arr:
        bars = bar(x, row, bottom=bottom, **kwds)
        stacked_bars.append(bars)
        if bottom is not None:
            bottom += row
        else:
            bottom = row
    return bars


def scatter_density(x, y, data=None, cmap='viridis', scale=None, ax=None, **kwargs):
    points = data[[x, y]].T if data is not None else np.array([x, y]) 
    density = gaussian_kde(points)(points)
    scatter = ax.scatter if ax is not None else plt.scatter
    collections = scatter(x, y, data=data, c=density, cmap=cmap, **kwargs)
    ax = ax if ax is not None else plt.gca()
    if scale is not None:
        ax.set_xscale(scale)
        ax.set_yscale(scale)
    if isinstance(x, str):
        xlabel = x
    elif isinstance(x, int):
        xlabel = 'column %d' % x
    elif isinstance(x, pd.Series):
        xlabel = x.name
    else:
        xlabel = 'X'
    if isinstance(y, str):
        ylabel = y
    elif isinstance(y, int):
        ylabel = 'column %d' % y
    elif isinstance(y, pd.Series):
        ylabel = y.name
    else:
        ylabel = 'Y'
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    return collections


def dendrogram():
    raise NotImplementedError
    #Z = hclust.linkage(df2, hc_method, hc_dist)
    #L = hclust.leaves_list(Z) # ordered list of leaf labels
    ddg = hclust.dendrogram(Z, orientation=orientation, labels=df.index)

    # Annotate nodes with the size of the union
    leaves = ddg['ivl']
    for x, leaf in zip(ax.get_xticks(), leaves):
        ax.text(x, 0, tree[leaf]['size'], va='bottom', ha='left',
                fontsize=7)
    ax.set_xticklabels(leaves, va='bottom', ha='right', rotation=90,
                       fontsize=8)
    ax.set_title(anc, style='italic', family='serif')
    ax.spines['top'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['right'].set_visible(False)


def plottree(tree, get_items, get_label, root=None, ax=None, *args, invert=True, topology_only=False, label_params=None, **kwargs):
    """Plot an ete3 tree, from left to right."""
    coord = namedtuple('coord', 'x y')

    if root is None:
        try:
            root = tree.clade  # .root
        except AttributeError:
            try:
                root = tree.root
            except AttributeError:
                root = tree  # .get_tree_root (ete3) but you usually want to take the current node.

    if topology_only:
        get_items_withdist = get_items
        def get_items(tree, nodedist):
            return [(child, 1) for child, _ in get_items_withdist(tree, nodedist)]

    #depth = tree.get_farthest_leaf()[1]
    leafdists = sorted(iter_distleaves(tree, get_items, root), key=lambda x: x[1])
    depth = leafdists[-1][1]  # furthest leaf.
    leafloc, leafstep = (0, 1) if invert is False else (len(leafdists)-1, -1)
    leafdists = dict(leafdists)
    try:
        rootdist = tree.dist  # ete3 instance
    except AttributeError:
        try:
            rootdist = tree.clade.branch_length  # Bio.Phylo
        except AttributeError:
            rootdist = getattr(tree, 'rootdist', None)  # myPhylTree
    if rootdist is None: rootdist = 0

    child_coords = {}  # x (node depth), y (leaf number)
    x = []
    y = []
    ticklabels = []

    extended_x = []  # Dashed lines to the right when tree is not ultrametric
    extended_y = []

    for (node,dist), items in rev_dfw_descendants(tree, get_items,
                                                  include_leaves=True,
                                                  queue=[(root, rootdist)]):
        if not items:
            # Is a leaf.
            child_coords[node] = coord(leafdists[node], leafloc)
            ticklabels.append(get_label(node))
            if child_coords[node].x < depth:
                extended_x.extend((depth, child_coords[node].x, None))
                extended_y.extend((leafloc, leafloc, None))
            leafloc += leafstep
        else:
            
            if len(items) == 1:
                (ch, chdist), = items
                child_coords[node] = nodecoord = coord(child_coords[ch].x - chdist,
                                                       child_coords[ch].y)
                x.extend((child_coords[ch].x, nodecoord.x, None))
                y.extend((child_coords[ch].y, nodecoord.y, None))
            else:
                sorted_items = sorted(items,
                                         key=lambda item: child_coords[item[0]].y)
                ch0, ch0dist = sorted_items[0]
                ch1, ch1dist = sorted_items[-1]
                child_coords[node] = nodecoord = coord(
                        child_coords[ch0].x - ch0dist,
                        (child_coords[ch0].y + child_coords[ch1].y)/2.)
                x.extend((child_coords[ch0].x,
                          nodecoord.x,
                          nodecoord.x,
                          child_coords[ch1].x, None))
                y.extend((child_coords[ch0].y,
                          child_coords[ch0].y,
                          child_coords[ch1].y,
                          child_coords[ch1].y, None))
                #forkwidth = child_coords[ch1] - child_coords[ch0].y
                #forkstep = forkwidth / (len(node.children) - 1.0)
                for extra_ch, _ in sorted_items[1:-1]:
                    x.extend((child_coords[extra_ch].x, nodecoord.x, None))
                    y.extend((child_coords[extra_ch].y,)*2 + (None,))
    if rootdist > 0:
        x.extend((0, -rootdist))
        y.extend((nodecoord.y, nodecoord.y))

    plot = plt.plot if ax is None else ax.plot

    if not kwargs.get('color') and (not args or not args[0][0] in 'bgrcmykw'):
        kwargs['color'] = mpl.rcParams['text.color']  # Default color to black.
    
    default_kwargs = {'clip_on': False}
    default_kwargs.update(kwargs)

    lines = plot(x, y, *args, **default_kwargs)
    if ax is None:
        ax = plt.gca()
    ax.plot(extended_x, extended_y, 'k--', alpha=0.4,
            linewidth=kwargs.get('linewidth', mpl.rcParams['lines.linewidth'])/2.)
    ax.set_xlim(min(0, 0-rootdist), depth)
    ax.spines['top'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.yaxis.tick_right()
    ax.tick_params('y', which='both', right=False)
    ax.set_yticks(range(len(ticklabels)))
    if label_params is None: label_params = {}
    if invert: ticklabels.reverse()
    ax.set_yticklabels(ticklabels, **label_params)
    ax.get_figure().tight_layout()  # extend interactive view to see labels.
    return lines



# Or joypy.joyplot
def kde_ridgeplot(x, by, data=None):
    """
    Overlapping densities ('ridge plot')
    ====================================
    
    From: https://github.com/mwaskom/seaborn/blob/master/examples/kde_ridgeplot.py
    """
    sns.set(style="white", rc={"axes.facecolor": (0, 0, 0, 0)})

    ## Create the data
    #rs = np.random.RandomState(1979)
    #x = rs.randn(500)
    #g = np.tile(list("ABCDEFGHIJ"), 50)
    #df = pd.DataFrame(dict(x=x, g=g))
    #m = df.g.map(ord)
    #df["x"] += m

    df = data

    # Initialize the FacetGrid object
    pal = sns.cubehelix_palette(10, rot=-.25, light=.7)
    g = sns.FacetGrid(df, row=by, hue=by, aspect=15, height=.5, palette=pal)

    # Draw the densities in a few steps
    g.map(sns.kdeplot, x, clip_on=False, shade=True, alpha=1, lw=1.5, bw=.2)
    g.map(sns.kdeplot, x, clip_on=False, color="w", lw=2, bw=.2)
    g.map(plt.axhline, y=0, lw=2, clip_on=False)


    # Define and use a simple function to label the plot in axes coordinates
    def label(x, color, label):
        ax = plt.gca()
        ax.text(0, .2, label, fontweight="bold", color=color,
                ha="left", va="center", transform=ax.transAxes)


    g.map(label, x)

    # Set the subplots to overlap
    g.fig.subplots_adjust(hspace=-.25)

    # Remove axes details that don't play well with overlap
    g.set_titles("")
    g.set(yticks=[])
    g.despine(bottom=True, left=True)
