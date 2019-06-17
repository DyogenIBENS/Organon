#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#%run -n ~/scripts/codeml/analyse/generate_dNdStable.py
from codeml.analyse.generate_dNdStable import *
import pandas as pd

logging.basicConfig(format=logging.BASIC_FORMAT)


# Output formatter for interactive examinations
def agesdf(*args, **kwargs):
    measures = kwargs.get('measures', ['dist'])
    node_info = kwargs.get('node_info', [])
    columns = ['name'] + ['branch_%s' %m for m in measures] + \
              ['age_%s' % m for m in measures] + ['cal'] + \
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

tree = ete3.Tree('(((A:0.8[&&NHX:dS=1],B:1[&&NHX:dS=1])x0:1[&&NHX:dS=1],C:2.2[&&NHX:dS=2])x:1[&&NHX:dS=1],D:3[&&NHX:dS=3])R;', format=1)

for node in tree.iter_descendants():
    node.dS = float(node.dS)

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


