#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd
from datasci.routines import leave1out_eval, loop_leave1out

# The function to evaluate:
def ascii_alphabet_index(letter):
    """Return 0 for 'a', 1 for 'b', ... 25 for 'z'."""
    if letter not in 'abcdefghijklmnopqrstuvwxyz':
        raise ValueError('Letter should be lower case alphabetic letter.')
    return ord(letter) - ord('a')


def test_eval_concat():
    features = list('abc')
    left1out = leave1out_eval(features, lambda features: ''.join(features), dtype=str)
    assert isinstance(left1out, pd.Series)
    assert left1out.shape == (3,)
    assert left1out['a'] == 'bc'
    assert left1out['b'] == 'ac'
    assert left1out['c'] == 'ab'

def test_eval_sumbytes():
    features = list('abc')
    data = pd.Series([4, 2, 1], index=features)
    left1out = leave1out_eval(features, lambda features: data[features].sum())
    assert isinstance(left1out, pd.Series)
    assert left1out.shape == (3,)
    assert left1out.index.tolist() == features
    assert left1out.values.tolist() == [3, 5, 6]

def test_loop_eval_nloops1_drops1():
    features = list('abc')
    data = pd.Series([4, 2, 1], index=features)
    dropped, _ = loop_leave1out(features, lambda features: data[features].sum(),
                                nloops=1, format_df=None)
    assert len(dropped)==1

def test_loop_eval_min_is_a_is_dropped():
    features = list('abc')
    data = pd.Series([4, 2, 1], index=features)
    dropped, _ = loop_leave1out(features, lambda features: data[features].sum(),
                                nloops=1, criterion='min', format_df=None)
    assert dropped[0] == 'a'


def test_loop_eval_min_all_equal_drops_first():
    features = list('abc')
    data = pd.Series([1, 1, 1], index=features)
    dropped, _ = loop_leave1out(features, lambda features: data[features].sum(),
                                nloops=1, criterion='min', format_df=None)
    assert dropped[0] == 'a'

def test_loop_eval_min_drops_nan():
    features = list('abc')
    data = pd.Series([1, np.NaN, 1], index=features)
    dropped, _ = loop_leave1out(features, lambda features: np.sum(data[features].values),
                                nloops=1, criterion='min', format_df=None)
    assert dropped[0] == 'a'
