#!/usr/bin/env python3
# -*- coding: utf-8 -*-


"""Statistics utilities."""
#  Should be called `basics`. `stats` for more complex operations.

from functools import partial
import numpy as np
import scipy.stats as stats
import matplotlib.pyplot as plt


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


def iqr(v):
    """Inter-quartile range"""
    q1, q3 = v.quantile([0.25, 0.75])
    return q3 - q1


def mad(v):
    """Median absolute deviation from the median"""
    return (v - v.median()).abs().median()


def partial_zscore(v, condition, ddof=1):
    return (v - v[condition].mean()) / v[condition].std(ddof=ddof)

def partial_quartiledispersion(v, condition):
    """Data centered around the median, scaled by the
    Inter-quartile range. Alternatively, one could normalise by the
    Median Absolute Deviation from the median."""
    q1, med, q3 = v[condition].quantile([0.25, 0.5, 0.75])
    return (v - med) / (q3 - q1)
    
def partial_maddispersion(v, condition):
    med = v[condition].median()
    mad = (v[condition] - med).abs().median()
    return (v - med) / mad

def rescale_groups(y, hue, by, data=None, ref_hue=None, center='mean',
                   scale='std', **scale_kw):
    """
    Rescale groups of measures Y to a comparable size:
    
    Data is split into groups (`by`), with each group subdivided by `hue`:
    the mean and stddev of `y` in each group (in `by`) is rescaled so that it
    is equal for the first hue value.

    The idea is to detect `hue` effects on `y` accross multiple groups (`by`).

    Works on `long` format (same as seaborn.violinplot).
    
    :param: `center`: 'mean' or 'median', or anything that .transform can accept.
    :param: `scale`: 'std', iqr, mad, or anything that .transform can accept.

    N.B: std default ddof is 0.
    """
    
    hue_levels = data[hue].unique()
    #assert hue_levels.size <= 2, "Not implemented for more than 2 hue values."
    #print(hue_levels)
    
    hue0 = hue_levels[0] if ref_hue is None else ref_hue

    # for each group (`by`), subtract the mean and divide by the SD of the subgroup 0.
    #return data.groupby(by, group_keys=False)\
    #                .apply(
    #                    lambda v: v.assign(**{y: scale(v[y], v[hue]==hue0, **scale_kw)})
    #                )
    
    # Not good for getting the overall variance of hues!=hue0 (because of intergroup variance)
    g = data.groupby([by, hue])
    centered = (data[y] - g[y].transform(center))

    if scale == 'std':
        scale_kw.setdefault('ddof', 0)

    scales = (data[y].where(data[hue].isin(hue0))
              if isinstance(hue0, (tuple, list, set, frozenset))
              else data[y].where(data[hue]==hue0))\
            .groupby(data[by]).transform(scale, **scale_kw)
    return centered / scales


# Must do a test of variance between 2 groups (blue/green):
# -> group all blue groups into one single one, normalizing each, and normalizing the
# corresponding green group by the same factor.
def multi_vartest(y, hue, by, data=None, ref_hue=None, center='mean',
                  scale='std', **scale_kw):
    # There must be only 2 `hue` values.
    # `by` is the "x" of the violinplot.

    if data is None:
        raise NotImplementedError
    rescaled_y = rescale_groups(y, hue, by, data, ref_hue, center, scale, **scale_kw)
    transformed_Y = data.assign(**{y:rescaled_y}).dropna(subset=[y])

    # Print details
    #if scale is partial_zscore:
    #    func = 'std'
    #    name = 'std'
    #elif scale is partial_quartiledispersion:
    #    func = lambda v: v.quantile(0.75) - v.quantile(0.25)
    #    name = 'iqr'
    #elif scale is partial_maddispersion:
    #    func = lambda v: (v - v.median()).abs().median()
    #    name = 'mad'
    scalename = scale.__name__ if callable(scale) else scale
    print(transformed_Y.groupby(hue)[y].agg(scale, **scale_kw).rename('%s(%s)' % (scalename, y)))
    #subgroups = transformed_Y.groupby([by, y]).std(ddof=0))
    return stats.levene(*(Y_group for _, Y_group in transformed_Y.groupby(hue)[y]))
    

def MStot(Y):
    """Mean of squared deviation from the average"""
    return Y.var(ddof=0)
    #return ( (Y - Y.mean())**2 ).mean()

def MSres(Y, pred):
    """Residual mean of squared errors (deviation from the predicted value)"""
    return ((Y - pred)**2).mean()


def r_squared(Y, pred):
    return 1 - (MSres(Y, pred) / MStot(Y))


def adj_r_squared(Y, pred, p):
    """p: Number of explanatory variables in the model"""
    n = len(Y)
    return 1 - (MSres(Y, pred) / MStot(Y)) * (n - 1)/(n - p -1)


def multicol_test(X):
    """Values >20 mean high colinearity."""
    norm_xtx = np.dot(X.T, X)

    eigs = np.linalg.eigvals(norm_xtx)
    #condition_number
    return np.sqrt(eigs.max() / eigs.min())


# Gamma distribution utilities
def expectation_gamma(shape, scale, loc=0):
    return loc + shape*scale


def conf_int_gamma(shape, scale, loc=0, conf=0.95):
    out = (1. - conf)/2
    return (stats.gamma.cdf(out, shape, loc, scale=scale),
            stats.gamma.sf(out, shape, loc, scale=scale))

def inv_conf_int_gamma(shape, scale, loc=0, conf=0.95):
    #out = (1. - conf)/2
    #return stats.gamma.ppf([out, conf+out], shape, loc, scale=scale)
    return stats.gamma.interval(conf, shape, loc, scale=scale)


def dist_from_conf_int_gamma(ci, shape, scale, loc, conf=0.95):
    """With loc and/or scale being an array"""
    inf, sup = stats.gamma.interval(conf, shape, loc, scale)
    return ((ci - np.stack((inf, sup)).T)**2).sum(axis=1)


def gamma_params_from_mean_and_ci(mean, ci, conf=0.95, shape_bounds=None,
                                  loc_bounds=None, nshapes=200, nlocs=None):
    # function to minimize: sum of squared differences from the requested interval boundaries
    ci = np.array(ci)
    func = partial(dist_from_conf_int_gamma, ci, conf=conf)

    if shape_bounds is None:
        shape_bounds = [0.025, 5.] #(ci[1] - ci[0])]
    shape_values = np.linspace(*shape_bounds, nshapes)

    if loc_bounds is None:
        loc_bounds = [ci[0]-(ci[1] - ci[0]), ci[0]]

    print('Loc is in %s' % loc_bounds)
    nlocs = nshapes if nlocs is None else nlocs
    loc_values = np.linspace(*loc_bounds, nlocs, endpoint=False)

    epsilon = 0.005
    objective = np.zeros((nshapes, nlocs))
    for i,x in enumerate(shape_values):
        scale_values = (mean-loc_values)/x
        assert np.allclose(stats.gamma.mean(x, loc_values, scale_values), mean), stats.gamma.mean(x, loc_values, scale_values)
        print('Shape = %s' % (x,))
        objective[i,:] = func(x, scale_values, loc_values)

    return objective, shape_values, loc_values


def show_obj(objective, shape_values, loc_values):
    #corners = (shape_values.min(), shape_values.max(),
    #           loc_values.min(), loc_values.max())
    im = plt.pcolormesh(loc_values, shape_values, objective,
                        cmap='viridis_r')
    ax = plt.gca()
    ax.set_xlabel('Loc')
    ax.set_ylabel('shape')
    #print(ax.get_xticks())
    #ax.set_xticklabels(loc_values[ax.get_xticks().astype(int)[:-1]])
    #ax.set_yticklabels(shape_values[ax.get_yticks().astype(int)[:-1]])
    ax.tick_params('x', labelrotation=45)
    plt.colorbar(im)
    return ax


def superimpose_gamma_params(mean, ci, shape, loc, conf=0.95, scale=None):
    if scale is None:
        scale = (mean-loc)/shape
    w = ci[1] - ci[0]
    x = np.linspace(ci[0] - w/10, ci[1] + w/10, 200)
    y = stats.gamma.pdf(x, shape, loc, scale)

    if shape < 1:
        x = np.concatenate(([x[0]], x))
        y = np.concatenate(([0], y))
    
    print('scale = %s' % scale,
          stats.gamma.mean(shape, loc, scale),
          stats.gamma.cdf(ci, shape, loc, scale),
          stats.gamma.interval(conf, shape, loc, scale))
    y2 = [y.min(), y.max()]
    plt.plot(x, y, '-',
             [ci[0]]*2, y2, '-',
             [mean]*2, y2, '-',
             [ci[1]]*2, y2, '-')

