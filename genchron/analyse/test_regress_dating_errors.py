#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import numpy as np
from genchron.analyse.regress_dating_errors import *
import logging
logger = logging.getLogger(__name__)


def test_old_group_average():
    #(g, var, weight_var="median_brlen")
    pass

def test_raw_group_average():
    #raw_group_average(g, var, weight_var)
    pass

def test_npraw_group_average():
    #npraw_group_average(gvalues, var_idx, weight_idx)
    pass

def test_npraw_group_average_w0():
    #npraw_group_average_w0(gvalues)
    gvalues = np.array([[1,1,1],
                        [1,2,3],
                        [4,1,9],
                        [0.1, 0.7, 0.7]]).T
    expected = [np.mean([1,2,3]), np.mean([4,1,9]), np.mean([0.1, 0.7, 0.7])]
    r = npraw_group_average_w0(gvalues)
    assert r.shape == (3,)
    assert np.allclose(r, expected)
    assert all(r[i] == expected[i] for i in range(3))
    return True


def test_group_average():
    #group_average(g, var, weight_var)
    pass

def test_group_weighted_std():
    #group_weighted_std(g, var, weight_var="median_brlen")
    pass

# for calculating rates: with dS, dN, t
def test_tree_dist_2_rate():
    #tree_dist_2_rate(g, dist_var, norm_var="median_brlen")
    pass


def test_compute_dating_errors():
    #compute_dating_errors(ages_controled, control='median', measures=['dS'],
    #                      rescale=None)
    pass


def test_compute_branchrate_std():
    #compute_branchrate_std(ages_controled, dist_measures,
    #                       branchtime='median_brlen_dS', taxon_age=None,
    #                       mean_condition=None,
    #                       std_condition=None)
    pass

def test_triplet_aggfunc():
    #triplet_aggfunc(func, func_args, func_kwargs, parentgrouped, row)
    pass

def test_raw_triplet_std():
    #raw_triplet_std(row, parentgrouped)
    pass


def test_compute_correlated_rate():
    #compute_correlated_rate(ages_controled, dist_measures,
    #                        branchtime='median_brlen_dS', taxon_age=None,
    #                        mean_condition=None,
    #                        std_condition=None):
    pass


