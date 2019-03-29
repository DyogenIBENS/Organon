#!/usr/bin/env python3
# -*- coding: utf-8 -*-


"""Statistics utilities."""
#  Should be called `basics`. `stats` for more complex operations.
import numpy as np
import scipy.stats as stats


def normal_fit(var):
    """Return normal density curve with same mean and variance as the given distribution"""
    sorted_var = np.sort(var)
    return sorted_var, stats.norm.pdf(sorted_var, sorted_var.mean(), sorted_var.std())
    # * np.isfinite(sorted_var).sum()


def cov2cor(cov):
    """Converts covariance matrix into correlation matrix"""
    var = np.diagonal(cov)[:,np.newaxis]
    return cov / np.sqrt(var.dot(var.T))


def car2pol(x,y):
    """Convert cartesian coordinates to polar (angle in radians)."""
    a = np.arctan2(x, y)
    r = np.sqrt(x*x + y*y)
    return a, r


def car2pol_deg(x,y):
    """Convert cartesian coordinates to polar (angle in degrees)."""
    a = np.arctan2(x, y) / np.pi * 180
    r = np.sqrt(x*x + y*y)
    return a, r


def f_test(x, y, data=None):
    """Do not use. Use Levene's test."""
    if data is not None:
        x, y = data[x], data[y]
    F = np.var(x, ddof=1) / np.var(y, ddof=1)
    dfx = len(x) - 1
    dfy = len(y) - 1
    f_mean = stats.f.mean(dfx, dfy)
    if F >= f_mean:
        p_further = stats.f.sf(F, dfx, dfy)
    else:
        p_further = stats.f.cdf(F, dfx, dfy)

    return p_further


def rescale_groups(y, hue, by, data=None):
    ### Works on `long` format (same as seaborn.violinplot).
    hue_levels = data[hue].unique()
    #assert hue_levels.size <= 2, "Not implemented for more than 2 hue values."
    print(hue_levels)
    
    hue0 = hue_levels[0]
    #grouped_y0 = data.loc[data[hue] == hue_levels[0], y].groupby(by)
    #Y_means = grouped_y0.agg('mean')
    #Y_stds  = grouped_y0.agg('std')
    
    #transformed_Y = [(data.loc[data[hue] == hue_level, y].groupby(by) - Y_means) / Y_stds
    #                 for hue_level in hue_levels]
    
    return data.groupby(by)\
                    .apply(
                        lambda v: v.assign(**{y: (v[y] - v.loc[v[hue] == hue0, y].mean())\
                                                / v.loc[v[hue] == hue0, y].std()})
                    )

# Must do a test of variance between 2 groups (blue/green):
# -> group all blue groups into one single one, normalizing each, and normalizing the
# corresponding green group by the same factor.
def multi_vartest(y, hue, by, data=None):
    # There must be only 2 `hue` values.
    # `by` is the "x" of the violinplot.
                    
    transformed_Y = rescale_groups(y, hue, by, data)
    return stats.levene(*(Y_group for _, Y_group in transformed_Y.groupby(hue)[y]))
    
