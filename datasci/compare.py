#!/usr/bin/env python3

"""Tools for comparing datasets"""

# Anc. Greek for `comparison`: 'Parabole'

from itertools import combinations
import numpy as np


def pairwise_intersections(*data, **nameddata):
    nameddata.update({i: datum for i, datum in enumerate(data)})
    
    difference_left = []
    intersections = []
    difference_right = []
    labels = []

    for (l1,d1), (l2,d2) in combinations(nameddata.items(), 2):
        s1, s2 = set(d1), set(d2)
        diff1 = s1.difference(s2)
        diff2 = s2.difference(s1)
        inter = s1 & s2

        if len(diff2) < len(diff1):
            # Reorder: put the smallest difference on the left.
            diff1, diff2 = diff2, diff1
            s1, s2 = s2, s1
            l1, l2 = l2, l1

        labels.append((l1, l2))
        difference_left.append(len(diff1))
        difference_right.append(len(diff2))
        intersections.append(len(inter))

    return labels, difference_left, intersections, difference_right


def align_sorted(x, y):
    aligned_x, aligned_y = [], []
    x_matched = np.searchsorted(y, x)
    # Spread y at consecutive identical indices
    #x_matchcounts = np.bincount(x_matched)

    #FIXME: few misalignments
    prev_m = 0
    for i,m in enumerate(x_matched):
        if m == prev_m:
            aligned_y.append(np.NaN)
            aligned_x.append(x[i])
        else:
            if m > prev_m+1:
                aligned_y.extend(y[(prev_m+1):m])
                aligned_x.extend([np.NaN]*(m - (prev_m+1)))
            aligned_y.append(y[m])
            aligned_x.append(x[i])
        prev_m = m
    return aligned_x, aligned_y

