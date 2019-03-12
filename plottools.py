#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Code snippets for matplotlib/pandas plotting."""

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
import scipy.cluster.hierarchy as hclust
import scipy.spatial.distance as spdist
import pandas as pd
from collections import namedtuple

import pandas as pd
import seaborn as sns

from .stattools import car2pol, cov2cor

from dendron.climber import rev_dfw_descendants, iter_distleaves
import matplotlib.patches as patches
from matplotlib.path import Path
MOVETO, CURVE3, LINETO = Path.MOVETO, Path.CURVE3, Path.LINETO
### Plotting style in matplotlib ###

import logging
logger = logging.getLogger(__name__)

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
    raise NotImplementedError('Use Seaborn')


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
    points = data[[x, y]].T.values if data is not None else np.array([x, y]) 
    if not np.isfinite(points).all():
        notfinite = (~np.isfinite(points)).any(axis=0)
        logger.warning('Dropping %d not finite points', notfinite.sum())
        points = points[~notfinite,:]
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


# Functions for the PCA

def plot_cov(ft_cov, features, cmap='seismic', figax=None, cax=None,
             ylabel="Features", cb_kw=None):
    """Plot a covariance matrix.
    
    Consider converting it to a correlation matrix if your PCA features were
    not normalized.
    """
    cmap = plt.get_cmap(cmap)
    norm = mpl.colors.Normalize(-1, 1)
    fig, ax = plt.subplots() if figax is None else figax
    img = ax.imshow(ft_cov, cmap=cmap, norm=norm, aspect='auto', origin='lower') #plt.pcolormesh
    ax.set_xticks(np.arange(len(features)))
    ax.set_yticks(np.arange(len(features)))
    ax.set_yticklabels(features)
    ax.set_xticklabels(features, rotation=45, ha='right')
    if ylabel:
        ax.set_ylabel(ylabel)
    ax.set_title("Feature covariance")
    if cb_kw is None: cb_kw = {}
    logger.debug('add heatmap colorbar')
    fig.colorbar(img, ax=None, #(ax if cax is None else None),
                 cax=cax, aspect=ft_cov.shape[0], **cb_kw)
    return img


def heatmap_cov(ft_cov, features=None, cmap='seismic',
                dendro_ratio=0.20, dendro_pad=0.1, cb_ratio=0.05, cb_pad=0.025):
    """plot_cov, but with hierarchical clustering on the side"""
    # Tested with figsize=(20, 12)
    if features is None:
        features = ft_cov.index.tolist()
        ft_cov = ft_cov.values
    fig, (ax_ddg, ax, ax_cb) = plt.subplots(1,3,
                                    subplot_kw={'facecolor': 'none'},
                                    gridspec_kw={'width_ratios': [
                                                  dendro_ratio,
                                                  1,
                                                  cb_ratio]})
    #(x0, y0), (w, h) = ax.get_position().get_points()
    ## absolute padding (in figure coordinate)
    ## correct ratio by taking pad into account
    
    ## position relatively to figure (percentages)
    #ax.set_position([x0 + (dendro_ratio+dendro_pad)*w, y0,
    #                 w*(1-dendro_ratio-dendro_pad-cb_ratio-cb_pad), h])
    ##width2 = width*ratio - float(pad)/w
    #ax_ddg = fig.add_axes([x0, y0, w*dendro_ratio, h], frameon=False,
    #                      sharey=ax)
    ax_ddg.get_shared_y_axes().join(ax)
    #ax_cb.get_shared_x_axes().join(ax)
    ax_ddg.set_title("hierarchical clustering (euclidean)")
    ax_ddg.axis('off')
    ax_ddg.xaxis.set_visible(False)
    #ax_cb = fig.add_axes([x0 + w*(1-cb_ratio), y0, w*cb_ratio, h])

    distmat = 1 - np.abs(cov2cor(ft_cov))
    tol=1e-15
    #assert (np.diag(distmat) < tol).all()
    #assert (np.abs(distmat - distmat.T) < tol).all()
    spdist.is_valid_dm(distmat, tol, throw=True)

    flatdist = spdist.squareform(distmat, force='tovector', checks=False)
    Z = hclust.linkage(flatdist, method='average', metric='euclidean')
    ddg = hclust.dendrogram(Z, orientation='left', no_labels=True, #labels=features,
                            ax=ax_ddg)

    clustered_ft_cov = ft_cov[ddg['leaves'],:][:,ddg['leaves']]
    #print(ddg['leaves'], ft_cov.shape)
    #print(clustered_ft_cov)
    logger.debug(np.array(features)[ddg['leaves']])
    logger.debug('%d artists: %s', len(ax_cb.get_children()), ax_cb.get_children())
    img = plot_cov(clustered_ft_cov,
             np.array(features)[ddg['leaves']], cmap, (fig, ax), ax_cb, '')#,
             #cb_kw={'fraction': cb_ratio/(1. - dendro_ratio), 'shrink': 0.5})
    #fig.colorbar(img, ax=None, #(ax if cax is None else None),
    #             cax=ax_cb, aspect=ft_cov.shape[0])
    logger.debug('%d artists: %s', len(ax_cb.get_children()), ax_cb.get_children())
    #ax_cb.xaxis.set_visible(True)
    ax_cb.set_ylabel('Correlation coefficient')
    box_cb = ax_cb.get_position()
    w_cb, h_cb = box_cb.size
    ax_cb.set_position(box_cb.translated(w_cb/2., h_cb/4.).shrunk(0.5,0.5).shrunk_to_aspect(ft_cov.shape[0]))

    #xmin_ylabel = min(yt.get_window_extent().x0 / fig.dpi for yt in ax.get_yticklabels())
    #logger.debug('xmin_ylabel: %s VS xmin (ax): %s', xmin_label,
    #             ax.get_position().x0)
    ax.set_position(ax.get_position().translated(1.5*w_cb, 0))
    ax_ddg.set_position(ax_ddg.get_position().translated(-w_cb, 0))
    #plt.tight_layout(pad=0)
    #plt.show()
    #return fig


def plot_loadings(components, cmap="PRGn"):
    """Not as great as the pandas dataframe styling, because the scale is on
    the entire image here."""
    cmap = plt.get_cmap(cmap)
    norm = mpl.colors.Normalize(components.min().min(), components.max().max())
    plt.imshow(components, cmap=cmap, norm=norm)
    ax = plt.gca()
    ax.set_xticks(np.arange(components.shape[1]))
    ax.set_yticks(np.arange(components.shape[0]))
    ax.set_yticklabels(components.index)
    ax.set_xticklabels(components.columns)
    plt.colorbar()
    ax.set_title("Feature contribution")


def annotate_features_radar(ax, components, features, PCs):

    rtxt = ax.get_rmax()
    seen_coords = np.zeros((components.shape[0], len(PCs)))
    tooclose = np.pi / 36  # 10 degrees.

    coords = pd.concat(car2pol(components[PCs[0]], components[PCs[1]]),
                       axis=1, keys=['a', 'r']).sort_values(['a', 'r'])

    # Get the density of the point angles:
    adensity = stats.gaussian_kde(coords.a, lambda gk: np.pi/18)(coords.a)
    # Spread angles

    for ft, coord in coords.iterrows():
        #ft_vect = components.loc[ft][PCs] * 0.1
        a, r = coord
        angle = a / (2*np.pi) * 360
        ha = 'left' if (-90 < angle <= 90) else 'right'
        va = 'bottom' if (angle>0) else 'top'
        # Text should not be upside down.
        rotation = angle if (-90 < angle <= 90) else angle-180
        ax.annotate(s=ft, xy=(a, r*1.05), xytext=(a, (r*1.05 + rtxt)/2), #xycoords='polar',
                    arrowprops={'arrowstyle':'->',
                                'linestyle':'dashed',
                                'alpha':0.5},
                    rotation=rotation, verticalalignment=va,
                    horizontalalignment=ha, alpha=0.8)
        #plt.text(ft_vect[0], ft_vect[1], ft)

    #ax.set_xlabel(PCs[0])
    #ax.set_ylabel(PCs[1])


def plot_features_PCspace(components, features, PCs=["PC1", "PC2"], ax=None):
    quiver = plt.quiver if ax is None else ax.quiver 
    quiver(0, 0, components[PCs[0]], components[PCs[1]],
           units='dots', width=1, scale_units='width')
           #units='xy', 
    if ax is None: ax = plt.gca()
    annotate_features_radar(ax, components, features, PCs)
    return ax


def plot_features_radar(components, features, PCs=['PC1', 'PC2'], ax=None):
    if ax is not None:
        assert ax.name == 'polar'
        polar = ax.plot
    else:
        polar = plt.polar
    polar(*car2pol(components[PCs[0]], components[PCs[1]]), '.')
    if ax is None: ax = plt.gca()
    annotate_features_radar(ax, components, features, PCs)
    return ax




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
