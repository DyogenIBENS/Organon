#!/usr/bin/env python3
# -*- coding: utf-8 -*-


"""Statistics utilities"""
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

