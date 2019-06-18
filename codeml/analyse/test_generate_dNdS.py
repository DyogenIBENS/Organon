#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#%run -n ~/scripts/codeml/analyse/generate_dNdStable.py
from codeml.analyse.generate_dNdStable import *
import pandas as pd

logging.basicConfig(format=logging.BASIC_FORMAT)
logger.setLevel(logging.DEBUG)

# Output formatter for interactive examinations
def agesdf(*args, **kwargs):
    measures = kwargs.get('measures', ['dist'])
    node_info = kwargs.get('node_info', [])
    columns = ['name'] + \
              ['%s_%s' %(s,m) for s in ('branch', 'age', 'p_clock')
                              for m in measures] + \
              ['cal'] + \
              [attr for attr,_ in node_info] + ['parent', 'root', 'treename']
    ages, _ = bound_average(*args, **kwargs)
    return pd.DataFrame(ages, columns=columns).set_index('name').drop('treename', 1)


# Test the MPL algo
def tocalibrate(node, *a, **kw):
    return node.name.startswith('x')


### BEGIN TESTS ###


tree = ete3.Tree('((A:1.2,B:0.8)x:1,C:2)R;', format=1)
calib = {'R': 2}

ages = agesdf(tree, calib,
              measures=['dist'],
              unweighted=False,
              method2=False,
              tocalibrate=tocalibrate,
              keeproot=False,
              allow_unequal_children_age=0,
              calib_selecter='name')

# Because of `keeproot=False`, nothing gets estimated
assert set('ABC') == set(ages.index)

ages = agesdf(tree, calib,
              measures=['dist'],
              unweighted=False, method2=False,
              tocalibrate=tocalibrate, keeproot=True,
              allow_unequal_children_age=0,
              calib_selecter='name')

assert ages.age_dist['x'] == 1

calib = {'R': 3}
tree = ete3.Tree('(((A:0.8,B:1)x0:1,C:2.2)x:1,D:3)R;', format=1)

ages = agesdf(tree, calib,
              measures=['dist'],
              unweighted=True, method2=False,
              tocalibrate=tocalibrate, keeproot=True,
              allow_unequal_children_age=0,
              calib_selecter='name')
assert ages.age_dist['x'] == 2
ages = agesdf(tree, calib,
              measures=['dist'],
              unweighted=False, method2=False,
              tocalibrate=tocalibrate, keeproot=True,
              allow_unequal_children_age=0,
              calib_selecter='name')
assert ages.age_dist['x'] > 2

calib['x'] = 1
ages = agesdf(tree, calib,
              measures=['dist'],
              unweighted=True, method2=False,
              tocalibrate=tocalibrate, keeproot=True,
              allow_unequal_children_age=0,
              calib_selecter='name')
assert ages.age_dist['x'] == 2
ages = agesdf(tree, calib,
              measures=['dist'],
              unweighted=False, method2=False,
              tocalibrate=tocalibrate, keeproot=True,
              allow_unequal_children_age=0,
              calib_selecter='name')
assert ages.age_dist['x'] > 2

tree = ete3.Tree('(((A:0.8[&&NHX:dS=0.01],B:1[&&NHX:dS=0.01])x0:1[&&NHX:dS=0.01],C:2.2[&&NHX:dS=0.02])x:1[&&NHX:dS=0.01],D:3[&&NHX:dS=0.03])R[&&NHX:S=100];', format=1)

for node in tree.iter_descendants():
    node.dS = float(node.dS)
tree.S = float(tree.S)

# Succeed all broadcasting operations?
ages = agesdf(tree, calib,
              measures=['dist', 'dS'],
              unweighted=False, method2=False,
              tocalibrate=tocalibrate, keeproot=True,
              allow_unequal_children_age=0,
              calib_selecter='name')
ages = agesdf(tree, calib,
              measures=['dist', 'dS'],
              unweighted=True, method2=False,
              tocalibrate=tocalibrate, keeproot=True,
              allow_unequal_children_age=0,
              calib_selecter='name')
ages = agesdf(tree, calib,
              measures=['dist', 'dS'],
              unweighted=False, method2=True,
              tocalibrate=tocalibrate, keeproot=True,
              allow_unequal_children_age=0,
              calib_selecter='name')
ages = agesdf(tree, calib,
              measures=['dist', 'dS'],
              unweighted=True, method2=True,
              tocalibrate=tocalibrate, keeproot=True,
              allow_unequal_children_age=0,
              calib_selecter='name')

# Compare method2 (scaling by keeping the original supporting branch length) VS method1
calib = {'R':2}
tree = ete3.Tree('((A:1,B:1)x:0.25,C:2.75)R;', format=1)
ages = agesdf(tree, calib,
              measures=['dist'],
              unweighted=False, method2=False,
              tocalibrate=tocalibrate, keeproot=True,
              allow_unequal_children_age=0,
              calib_selecter='name')
assert ages.age_dist['x'] == 1.0

ages = agesdf(tree, calib,
              measures=['dist'],
              unweighted=False, method2=True,
              tocalibrate=tocalibrate, keeproot=True,
              allow_unequal_children_age=0,
              calib_selecter='name')
assert ages.age_dist['x'] == 1.6  # 1/1.25 * 2


# Avoid negative branch lengths

tree = ete3.Tree('((A:1,B:1)x:0.25,C:0.3)R;', format=1)
ages = agesdf(tree, calib,
              measures=['dist'],
              unweighted=False, method2=False,
              tocalibrate=tocalibrate, keeproot=True,
              allow_unequal_children_age=0,
              calib_selecter='name')
assert ages.age_dist['x'] == 2.0

# Numeric example of Britton 2002
tree = ete3.Tree('(Uvularia perfoliata:15,(Uvularia Pudica:8,Disporum:20)x0:3)R:6;', format=1)
## Check p-values of unclocklikeliness
ages = agesdf(tree, calib,
              measures=['dist'],
              unweighted=True, method2=False,
              tocalibrate=tocalibrate, keeproot=True,
              allow_unequal_children_age=0,
              calib_selecter='name')
# s²(x0) = (8+20)/4 = 7  =>  s(x0) = 2.64  OK.
# s²(R) = ((8+20 + 15) + 3×2²) / 3² = 55 / 9 = 6.11
#
# Also the sd of the MPL of 'R' should be 2.47  OK.

# Also s²(R) = (s²(x0) + 3)×2² + 15) / 3²

# s²(R->x0) = (8+20+3×2²) / 2² = 10   ---> ERROR in the paper!!
# s²(R->C) = 15
# s(x1 - x2) = sqrt(15+10) = 5
# P(|15-17|/sqrt(15+10) > 1.96) = 0.689
assert round(ages.p_clock_dist['R'], 3) == 0.689

ages = agesdf(tree, calib,
              measures=['dist'],
              unweighted=False, method2=False,
              tocalibrate=tocalibrate, keeproot=True,
              allow_unequal_children_age=0,
              calib_selecter='name')
# No difference for the unweighted algo if ≤3 leaves
