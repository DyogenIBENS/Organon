#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Code snippets for matplotlib/pandas plotting."""

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt


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
    df_grouped = df.groupby(by)
    ngroups = len(df_grouped)
    nrows = ngroups // ncols + 1*(ngroups % ncols > 0)
    fig, axes = plt.subplots(nrows, ncols, sharex, sharey)
    #print(kde_axes, type(kde_axes))
    for i, (title, data) in enumerate(df_grouped):
        ax = axes[i // ncols, i % ncols]
        subkwds = dict((key,v[i]) if isinstance(v, list)
                       else (key, v[title]) if isinstance(v, dict)
                       else(key,v)
                       for key,v in kwds.items())
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
        
