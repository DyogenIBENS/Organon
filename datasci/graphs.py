#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Code snippets for matplotlib/pandas plotting."""

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.collections as mc
from scipy.stats import gaussian_kde
import scipy.cluster.hierarchy as hclust
import scipy.spatial.distance as spdist
import pandas as pd
from collections import namedtuple

import pandas as pd
import seaborn as sns

from .stats import car2pol, cov2cor

from dendro.bates import rev_dfw_descendants, iter_distleaves
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
def stackedbar(x, arr, ax=None, zero=0, **kwds):
    assert isinstance(arr, np.ndarray)
    if kwds.get('orientation', None) == 'horizontal':
        bar = plt.barh if ax is None else ax.barh
        base = 'left'
    else:
        bar = plt.bar if ax is None else ax.bar
        base = 'bottom'

    stacked_bars = []
    bottom = np.zeros(len(x))
    bottom -= arr[:zero].sum(axis=0)
    for row in arr:
        kwds.update({base: bottom})
        bars = bar(x, row, **kwds)
        stacked_bars.append(bars)
        if bottom is not None:
            bottom += row
    return stacked_bars


def scatter_density(x, y, data=None, cmap='viridis', scale=None, ax=None, **kwargs):
    points = data[[x, y]].T.values if data is not None else np.array([x, y]) 
    logger.debug("Shape of points: %s", points.shape)
    if not np.isfinite(points).all():
        notfinite = (~np.isfinite(points)).any(axis=0)
        logger.warning('Dropping %d not finite points', notfinite.sum())
        points = points[:,~notfinite]
    density = gaussian_kde(points)(points)
    scatter = ax.scatter if ax is not None else plt.scatter
    collections = scatter(points[0], points[1], c=density, cmap=cmap, **kwargs)
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
    norm = mpl.colors.Normalize(-1 if (ft_cov<0).any() else 0, 1)
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
    adensity = gaussian_kde(coords.a, lambda gk: np.pi/18)(coords.a)
    # Spread angles

    for ft, coord in coords.iterrows():
        #ft_vect = components.loc[ft][PCs] * 0.1
        a, r = coord
        angle = a / (2*np.pi) * 360
        ha = 'left' if (-90 < angle <= 90) else 'right'
        va = 'bottom' if (angle>0) else 'top'
        # Text should not be upside down.
        rotation = angle if (-90 < angle <= 90) else angle-180
        ax.annotate(ft, xy=(a, r*1.05), xytext=(a, (r*1.05 + rtxt)/2), #xycoords='polar',
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


## Create colormaps

def colortuple_to_hex(coltup):
    return '#' + ('%02x' * len(coltup)) % tuple(int(round(c*255)) for c in coltup)

def value2color(values, cmap='afmhot', extend=1):
    # Setup the colormap for evolutionary rates
    cmap = plt.get_cmap('afmhot')
    value_range = values.max() - values.min()
    norm = mpl.colors.Normalize(values.max() - extend*value_range,
                                values.min() + extend*value_range)

    return values.apply(lambda v: colortuple_to_hex(cmap(norm(v))))


def plottree(tree, get_items, get_label, root=None, ax=None, invert=True,
             age_from_root=False,
             topology_only=False,
             label_params=None, edge_colors=None,
             edge_cmap='afmhot', add_edge_axes=None, **kwargs):
    """Plot an ete3 tree, from left to right.
    
    param: edge_colors dict-like object with keys being the nodes, and values ~~a color string~~
            a scalar value mapped to a color using a cmap.
    param: add_edge_axes can be None, "top", or "middle".
    """
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

    if edge_colors is not None:
        extend = 1.1
        edge_range = edge_colors.max() - edge_colors.min()
        #edge_norm = mpl.colors.Normalize(edge_colors.max() - extend*edge_range,
        #                            edge_colors.min() + extend*edge_range)


    #depth = tree.get_farthest_leaf()[1]
    leafdists = sorted(iter_distleaves(tree, get_items, root), key=lambda x: x[1])
    depth = leafdists[-1][1]  # furthest leaf.
    leafloc, leafstep = (0, 1) if invert is False else (len(leafdists)-1, -1)
    if age_from_root:
        root_age = 0
        leafdists = dict(leafdists)
    else:
        root_age = -depth
        leafdists = {l: ld-depth for l, ld in leafdists}
    try:
        rootdist = tree.dist  # ete3 instance
    except AttributeError:
        try:
            rootdist = tree.clade.branch_length  # Bio.Phylo
        except AttributeError:
            rootdist = getattr(tree, 'rootdist', None)  # myPhylTree
    if rootdist is None: rootdist = 0

    child_coords = {}  # x (node depth), y (leaf number)
    #xy = []  # coords to be unpacked and given to plot: plt.plot(*xy)
    segments = []
    line_edge_values = []

    ticklabels = []

    extended_x = []  # Dashed lines to the right when tree is not ultrametric
    extended_y = []

    if ax is None:
        fig, ax = plt.subplots()
    else:
        fig = ax.get_figure()

    axes_to_add = []  # List of args to be given to fig.add_axes:
                      # rectangle positions in data coordinates:
                      # [left, bottom, width, height]

    for (node,dist), items in rev_dfw_descendants(tree, get_items,
                                                  include_leaves=True,
                                                  queue=[(root, rootdist)]):
        if not items:
            # Is a leaf.
            child_coords[node] = coord(leafdists[node], leafloc)
            ticklabels.append(get_label(node))
            if child_coords[node].x < depth+root_age:
                extended_x.extend((depth+root_age, child_coords[node].x, None))
                extended_y.extend((leafloc, leafloc, None))
            leafloc += leafstep
        else:
            
            if len(items) == 1:
                (ch, chdist), = items
                child_coords[node] = nodecoord = coord(child_coords[ch].x - chdist,
                                                       child_coords[ch].y)
                #xy += [(child_coords[ch].x, nodecoord.x)
                #       (child_coords[ch].y, nodecoord.y)]
                # segments for LinesCollection
                segments.append(
                            [(child_coords[ch].x, child_coords[ch].y),
                             (nodecoord.x,        nodecoord.y)])
                if edge_colors is not None:
                #    xy.append(edge_colors[ch])
                    line_edge_values.append(edge_colors[ch])
            else:
                sorted_items = sorted(items,
                                      key=lambda item: child_coords[item[0]].y)
                ch0, ch0dist = sorted_items[0]
                ch1, ch1dist = sorted_items[-1]
                child_coords[node] = nodecoord = coord(
                        child_coords[ch0].x - ch0dist,
                        (child_coords[ch0].y + child_coords[ch1].y)/2.)
                for ch,chdist in sorted_items:
                    #xy += [(child_coords[ch].x, nodecoord.x, nodecoord.x),
                    #       (child_coords[ch].y, child_coords[ch].y, nodecoord.y)]
                    segments.append(
                             [(child_coords[ch].x, child_coords[ch].y),
                              (nodecoord.x, child_coords[ch].y),
                              (nodecoord.x, nodecoord.y)])
                    if edge_colors is not None:
                        line_edge_values.append(edge_colors[ch])
                    if add_edge_axes:
                        # Assuming drawing tree left->right and bottom->top
                        shift = -0.5 if add_edge_axes == 'middle' else 0
                        axes_to_add.append((ch,
                                            [nodecoord.x,
                                             child_coords[ch].y + shift,
                                             chdist,
                                             1]))
    if rootdist > 0:
        #xy += [(0, -rootdist),
        #       (nodecoord.y, nodecoord.y)]
        segments.append([(root_age,          nodecoord.y),
                         (root_age-rootdist, nodecoord.y)])
        if edge_colors is not None:
        #    xy.append(edge_colors.get(root))
            line_edge_values.append(edge_colors.get(root, np.NaN))
        if add_edge_axes:
            # Assuming drawing tree left->right and bottom->top
            shift = -0.5 if add_edge_axes == 'middle' else 0
            axes_to_add.append((root,
                                [root_age-rootdist,
                                 nodecoord.y + shift,
                                 root_age,
                                 nodecoord.y + 1 + shift]))

    #if edge_colors is None:
    #    # or not args[0][0] in 'bgrcmykw'):
    #    edge_cmap = None
    line_color = kwargs.pop('color', mpl.rcParams['text.color'])
    #else:
    #    line_color = None

    #lines = plot(x, y, *args, **default_kwargs)
    #lines = ax.plot(*xy, **default_kwargs)
    lines = mc.LineCollection(segments, colors=line_color, cmap=edge_cmap)#, norm=edge_norm)

    if edge_colors is not None:
        lines.set_array(np.array(line_edge_values))  # Value to be converted to colors.

    lines.set_clip_on(False)
    ax.add_collection(lines)
    ax.set(**kwargs)

    ax.plot(extended_x, extended_y, 'k--', alpha=0.4,
            linewidth=kwargs.get('linewidth', mpl.rcParams['lines.linewidth'])/2.)

    ax.set_xlim(min(root_age, root_age-rootdist), depth+root_age)
    ax.set_ylim(-0.5, len(leafdists))
    #ax.autoscale_view()
    
    ax.spines['top'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.yaxis.tick_right()
    ax.tick_params('y', which='both', right=False)
    ax.set_yticks(range(len(ticklabels)))
    if label_params is None: label_params = {}
    if invert: ticklabels.reverse()
    ax.set_yticklabels(ticklabels, **label_params)
    fig = ax.get_figure()

    #if lines.get_array() is not None:
    #    cax,kw = mpl.colorbar.make_axes_gridspec(ax, panchor=(0,1))  # no effect of panchor.
    #    fig.colorbar(lines, cax=cax, **kw)
    #fig.tight_layout()  # extend interactive view to see labels, but before subplots.
    
    #if edges_colors is not None:
        # Would be better directly using the cmap and norm.
        #scalarmappable = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        #sm.set_array([])
        # Inside the ax would be nicer.
        #fig.colorbar(sm)

    # Transform coords from data to display, then from display to figure.
    data2fig = ax.transData + fig.transFigure.inverted()

    subaxes = {}
    for ch, rect in axes_to_add:
        xlim = rect[0], rect[0] + rect[2]
        subax = fig.add_axes(rect, xlim=xlim, position=data2fig.transform(rect),
                             frame_on=False, autoscale_on=False
                             ) #share_y #clip_on=False, visible=False
        subax.set_clip_on(False)
        subax.axis('off')
        subaxes[ch] = subax

    return lines, child_coords, subaxes


def toothpaste(x, y, hue, data=None, order=None, hue_order=None,
               cmap='viridis', **kwargs):
    # Suggested colormaps: 'viridis', 'gist_rainbow', 'rainbow'
    default_kwargs = {'widths': 1.}
    default_kwargs.update(kwargs)
    
    g = data.groupby(x)
    if order is None:
        order = sorted(g.groups)
    if hue_order is None:
        hue_order = sorted(data.groupby(hue).groups)

    cmap = plt.get_cmap(cmap, len(hue_order))

    dataset = [g.get_group(name)[y].values for name in order]
    violins = plt.violinplot(dataset=dataset, **default_kwargs)
    ax = plt.gca()

    facecolor='none'
    edgecolor='k'

    #plt.setp(ax.lines, zorder=100)
    plt.setp(ax.collections, edgecolor=edgecolor, alpha=0.8)
    
    for i, polycollection in enumerate(violins['bodies']):
        label = order[i]
        polycollection.set_facecolor(facecolor)
        polycollection.set_edgecolor(edgecolor)
        polycollection.set_alpha(0.8)
        
        violinpath, = polycollection.get_paths()
        coords = violinpath.vertices
        #x = i+1
        #left_coords = coords[coords[:,0]  < x]
        #right_coords = coords[coords[:,0] > x]
        mid = coords.shape[0] // 2
        left_coords = coords[:mid]
        if coords.shape[0] % 2:
            # Skip middle coords if it draws a horizontal line.
            mid += 1
        right_coords = coords[mid:]

        hue_left_x = left_coords[:,0].copy()
        violin_y = left_coords[:,1]
        
        # Compute the density for each hue value
        violindata = g.get_group(label)

        # Estimate the width-scaling (Did not find a way to avoid recomputing the density)
        scaling = float(default_kwargs['widths']) / max(gaussian_kde(violindata[y])(violin_y))
        scaling /= violindata.shape[0]

        violindata_hues = violindata.groupby(hue)

        for j, hue_value in enumerate(hue_order):
            logger.debug('violin #%d: hue #%d %r', i, j, hue_value)
            try:
                hue_data = violindata_hues.get_group(hue_value)
            except KeyError as err:
                logger.debug(str(err))
                continue

            npoints = hue_data.shape[0]
            hue_density = npoints*scaling * gaussian_kde(hue_data[y])(violin_y)
            ax.fill_betweenx(violin_y, hue_left_x, hue_left_x + hue_density, color=cmap(j))
            hue_left_x += hue_density

        #ax.plot(left_coords[:,0], violin_y, color=edgecolor)
        #ax.plot(hue_left_x, violin_y, color=edgecolor)


    return violins


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
