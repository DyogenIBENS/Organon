#!/usr/bin/env python3

import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.patches as patches
import matplotlib.pyplot as plt
from IPython.display import display_html


def centered_background_gradient(s, cmap='PRGn', center=0, extend=0):
    """Color background in a range according to the data, centered around the given value.
    Adapted from pandas.io.formats.style.Styler.background_gradient()"""
    smin, smax = s.min(), s.max()
    assert smin <= center and center <= smax
    most_distant_absval = max(center - smin, smax - center)
    rng = 2 * most_distant_absval
    # extend lower / upper bounds, compresses color range
    norm = mpl.colors.Normalize(center - most_distant_absval - (rng * extend),
                                center + most_distant_absval + (rng * extend))
    # matplotlib modifies inplace?
    # https://github.com/matplotlib/matplotlib/issues/5427
    normed = norm(s.values)
    c = [mpl.colors.rgb2hex(x) for x in mpl.cm.get_cmap(cmap)(normed)]
    return ['background-color: {color}'.format(color=color) for color in c]


def magnify():
    """Increase the size of the table cells."""
    return [dict(selector="th",
                 props=[("font-size", "4pt"),
                        ("writing-mode", "bt-rl")]),
            dict(selector="td",
                 props=[('padding', "0em 0em")]),
            dict(selector="th:hover",
                 props=[("font-size", "10pt")]),
            dict(selector="tr:hover td:hover",
                 props=[('max-width', '150px'),
                        ('font-size', '10pt')])
            ]


def tolong(df, x=None, subset=None, varname=None, idname=None):
    """Stack operation while keeping one column as "x".
    
    `x`:      original column(s) to keep in the long dataframe.
    `subset`: columns to be stacked.
    `varname`: name of the new stacked column.
    `idname`:  name of the new column identifying the original column name.

    NOTE: to keep a variable telling the row of origin of the observation, 
    do a `df.reset_index()` first (and add its name to `x`).
    """
    Nx = 1
    if x is not None:
        df = df.set_index(x)
        Nx = 1 if isinstance(x, str) else len(x)
        #else: the original index will become the "level_0" column.
    if subset is not None:
        df = df[subset]
    
    renames = {}
    if varname is not None:
        renames[0] = varname
    if idname is not None:
        renames['level_%d' % Nx] = idname

    return df.stack().reset_index().rename(columns=renames)


def matplotlib_stylebar(data, y=None, color='#d65f5f', horizontal=True,
                        err=None, ticks=False):
    """Like pandas.DataFrame.style.bar, but with matplotlib.
    
    + Also: error bars
    
    Return axes (flat).
    """
    if y is None:
        try:
            y = data.columns.tolist()
        except AttributeError:  # it is a Series.
            y = data.name
            if y is None:
                y = 'Y'
            data = data.to_frame(y)
    if not isinstance(y, list):
        y = [y]
    if horizontal:
        kind = 'barh'
        sharex, sharey = False, True
        layout = (1, len(y))
        text = lambda ax, pos, val: ax.text(ax.get_xlim()[1], pos, '%.5f' % val,
                                            va='center', ha='right')
        xerr, yerr = err, None
    else:
        kind = 'bar'
        sharex, sharey = True, False
        layout = (len(y), 1)
        text = lambda ax, pos, val: ax.text(pos, ax.get_ylim()[1], '%.5f' % val, ha='center',
                                            va='top', rotation=90)
        xerr, yerr = None, err

    axes = data.plot(kind=kind, y=y, color='#d65f5f', width=0.95,
                     subplots=True, legend=False, sharex=sharex, sharey=sharey,
                     layout=layout, yerr=yerr, xerr=xerr,
                     error_kw={'capsize': 3, 'elinewidth': 1, 'alpha': 0.7})#'ecolor': mpl.rcParams['text.color']})
    if horizontal:
        axes.flat[0].invert_yaxis()
    for i, (ax, yvar) in enumerate(zip(axes.flat, y)):
        ymin, ymax = data[yvar].min(), data[yvar].max()
        if err is not None and np.ndim(err)==3:
            max_low_err = err[i,0].max()
            max_up_err = err[i,1].max()
            if not np.isnan(max_low_err):
                ymin -= max_low_err
            if not np.isnan(max_up_err):
                ymax += max_up_err
        #if align=='left':
        #    datalim = 0, ymax-ymin
        #    #bottom = ymin, heights = data[yvar] - ymin
        #elif align=='zero':
        #    maxabs = max(abs(ymin), abs(ymax))
        #    datalim = (-maxabs, maxabs)
        #elif align=='mid':
        #    datalim = ymin, ymax

        if horizontal:
            ax.set_xlim(ymin, ymax)
            if not ticks: ax.set_xticklabels([])
            #ax.grid(False)
            ax.spines['top'].set_visible(True)
            ax.spines['bottom'].set_visible(False)
            ax.xaxis.set_label_position('top')
        else:
            ax.set_ylim(ymin, ymax)
            if not ticks: ax.set_yticklabels([])
            ax.set_title(None)
            ax.set_ylabel(yvar)
        ax.grid(False)
        ax.tick_params(left=False, bottom=False)
        # Add the table content (Once the xlim/ylim is set!)
        for j, val in enumerate(data[yvar].values):
            text(ax, j, val)
    return axes.flat


def matplotlib_background_gradient(data, cmap='YlGn', axis=None,
                                   float_fmt='%.4f',
                                   surround=None, cbar=True, sep=True, ax=None,
                                   cbar_kw=None, **style_kwargs):
    #data = styled_data.data.xs('adjR2', axis=1, level=1)  # select from multiIndex
    orig_data = data
    if axis in (0, 1):
        broadcast_axis = 1 - axis
        data = data.subtract(data.min(axis=axis), broadcast_axis)\
                   .div(np.ptp(data.values, axis=axis), axis=broadcast_axis)
    #im = plt.imshow(data, cmap=cmap, aspect='auto', interpolation='none'); ax = plt.gca()
    #pcolor is the most precise on rectangle boundaries
    pcolor = plt.pcolor if ax is None else ax.pcolor
    im = pcolor(data, cmap=cmap, edgecolors='', snap=True)
    if ax is None:
        ax = plt.gca()
    ax.invert_yaxis()
    yticks = np.arange(data.shape[0]) + 0.5
    xticks = np.arange(data.shape[1]) + 0.5
    ax.xaxis.tick_top()
    ax.xaxis.set_label_position('top')

    ax.set_yticks(yticks)
    ax.set_yticklabels(data.index, {'fontweight': 'bold'})

    ax.set_xticks(xticks)
    ax.set_xticklabels(data.columns, {'fontweight': 'bold'})
    ax.tick_params(axis='y', labelsize='large')

    ax.spines['top'].set_visible(True)
    ax.spines['bottom'].set_visible(False)
    ax.grid(False)
    ax.tick_params(bottom=False, top=False, right=False, left=False)
    #ax.set_xlabel('Measure of rate heterogeneity')
    #ax.set_xticklabels(['%s\n%s' % t for t in all_adjR2.columns.values])

    #fig.colorbar(im)

    # Now, fix the text color according to the background color (copy from Pandas style)
    styled_data = orig_data.style.background_gradient(cmap=cmap, axis=axis, **style_kwargs)
    styled_data.render()  # To get the text colors.
    #print(styled_data.ctx)

    cell_txts = []
    for x,xt in enumerate(xticks):
        for y,yt in enumerate(yticks):
            attributes = styled_data.ctx[(y, x)]  # if it's not a MultiIndex
            dict_attr = dict(t.split(':') for t in attributes)
            try:
                textcolor = mpl.colors.to_rgba(dict_attr['color'].strip())
            except KeyError as err:
                err.args += ((y, x), attributes)
                raise
            cell_txts.append(
                    ax.text(xt, yt, float_fmt % orig_data.iloc[y,x],
                            color=textcolor, va='center', ha='center'))

    if cbar:
        orient = 'horizontal' if axis==1 else 'vertical'
        cbar_kw = {'fraction':0.05, 'pad':0.03, 'aspect':40,
                    **({} if cbar_kw is None else cbar_kw)}
        cbar = plt.colorbar(im, orientation=orient, **cbar_kw)
        if axis is not None:
            cbar.set_ticks([0, 1])
            cbar.set_ticklabels(['Row minimum', 'Row maximum'] if axis==1 else
                                ['Column minimum', 'Column maximum'])
        #if cbar_label is not None:
            #cbar.set_label('Relative dispersion (%s-wise)')
    #else:
    #    cbar = None

    if sep and axis is not None:
        if axis==1:
            sep_coords = range(1, orig_data.shape[0])
            plot_sep = ax.axhline
        else:
            sep_coords = range(1, orig_data.shape[1])
            plot_sep = ax.axvline
        for coord in sep_coords:
            plot_sep(coord, color=mpl.rcParams['axes.edgecolor'], linewidth=0.5,
                     alpha=1)

    if surround in ('min', 'max'):
        if axis is not None:
            get_coords = np.argmax if surround=='max' else np.argmin
            surround_coords = get_coords(orig_data.values, axis=axis)
            #print('DEBUG: surround_coords =', surround_coords)
            iter_coords = (zip(surround_coords, np.arange(orig_data.shape[1]))
                           if axis==0 else
                           zip(np.arange(orig_data.shape[0]), surround_coords))
            for i, j in iter_coords:
                #print('DEBUG: i, j = %s, %s' % (i, j))
                ax.add_patch(patches.Rectangle((j+0.002, i+0.005), 0.99, 0.95, fill=False,
                             edgecolor='gold', linestyle='--', linewidth=1, hatch='/'))
                             #zorder=3)) # edgecolor='khaki'

    # Get cell width/height
    #w, h = cell_txts[-1].get_window_extent()  # In display coordinates
    #print(w, h)
    return ax, cbar

