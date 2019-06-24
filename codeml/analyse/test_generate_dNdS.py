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
    node_ft_set = kwargs.get('node_feature_setter', [])
    columns = ['name'] + \
              ['%s_%s' %(s,m) for s in ('branch', 'age', 'p_clock')
                              for m in measures] + \
              ['cal', 'parent'] + \
              [attr for attr,_ in node_info] + [ft for ft,_ in node_ft_set] + \
              ['root', 'treename']
    ages, _ = bound_average(*args, **kwargs)
    return pd.DataFrame(ages, columns=columns).set_index('name').drop('treename', 1)


# Test the MPL algo
def todate(node, *a, **kw):
    return node.name.startswith('x')

### BEGIN TESTS ###

#### Utility functions ####

tree = ete3.Tree('((ENSG0000,ENSG0001)Homo.sapiensENSGT,Pan.paniscusENSGT)HomoPanENSGT;', format=1)

assert 'Homo sapiens' == get_taxon(tree&'ENSG0000', 93)
assert 'Homo sapiens' == get_taxon(tree&'ENSG0001', 93)
assert 'Homo sapiens' == get_taxon(tree&'Homo.sapiensENSGT', 93)
assert 'Pan paniscus' == get_taxon(tree&'Pan.paniscusENSGT', 93)
assert 'HomoPan' == get_taxon(tree&'HomoPanENSGT', 93)

subtree = {
        'ENSG0000': {'taxon': 'Homo sapiens'},
        'ENSG0001': {'taxon': 'Homo sapiens'},
        'Homo.sapiensENSGT': {'taxon': 'Homo sapiens'},
        'Pan.paniscusENSGT': {'taxon': 'Pan paniscus'},
        'HomoPanENSGT': {'taxon': 'HomoPan'}
        }

assert not isdup(tree&'ENSG0000', subtree)
assert not isdup(tree&'ENSG0001', subtree)
assert isdup(tree&'Homo.sapiensENSGT', subtree)
assert not isdup(tree&'Pan.paniscusENSGT', subtree)
assert not isdup(tree&'HomoPanENSGT', subtree)

def get_eventtype(node, subtree):
    if node.is_leaf():
        return 'leaf'
    return 'dup' if isdup(node, subtree) else 'spe'

assert 'leaf' == get_eventtype(tree&'ENSG0000', subtree)
assert 'leaf' == get_eventtype(tree&'ENSG0001', subtree)
assert 'dup' == get_eventtype(tree&'Homo.sapiensENSGT', subtree)
assert 'leaf' == get_eventtype(tree&'Pan.paniscusENSGT', subtree)
assert 'spe' == get_eventtype(tree&'HomoPanENSGT', subtree)

#### Recursive algo ####
#def get_taxon(node):
#    return node.name

tree = ete3.Tree('((A:1.2,B:0.8)x:1,C:2)R;', format=1)
calib = {'R': 2}

ages = agesdf(tree, calib,
              measures=['dist'],
              todate=todate,
              unweighted=False,
              original_leading_paths=False,
              correct_unequal_calibs='default',
              keeproot=False,
              calib_selecter='name')

# Because of `keeproot=False`, nothing gets estimated
assert set('ABC') == set(ages.index)

ages = agesdf(tree, calib,
              measures=['dist'],
              unweighted=False, original_leading_paths=False,
              todate=todate, keeproot=True,
              correct_unequal_calibs='default',
              calib_selecter='name')

assert ages.age_dist['x'] == 1
assert round(ages.p_clock_dist['x'],3) == 0.777  # 2*st.norm.sf(0.4, 0, np.sqrt(2))

# Compare with other corrections:
for correct_uneq_calibs in ('ignore', 'mine', 'pathd8'):
    assert agesdf(tree, calib,
                  measures=['dist'],
                  unweighted=False, original_leading_paths=False,
                  todate=todate, keeproot=True,
                  correct_unequal_calibs=correct_uneq_calibs,
                  calib_selecter='name').age_dist['x'] == 1
# => No effect when there is no extra calibration point.

calib = {'R': 3, 'y': 2, 'y0': 1}
tree = ete3.Tree('(((A:0.8,B:1)x0:1,C:2.2)x:1,D:3)R;', format=1)

def get_name(node, subtree):
    return node.name

# Test the `node_feature_setter`
ages = agesdf(tree, calib,
              measures=['dist'],
              unweighted=True, original_leading_paths=False,
              todate=todate, keeproot=True,
              correct_unequal_calibs='default',
              calib_selecter='name', node_feature_setter=[('name2', get_name)])
assert (ages.name2 == ages.index).all(), "Node features improperly set."
assert ages.age_dist['x'] == 2
assert round(ages.p_clock_dist['x'],4) == 0.8752  # 2*st.norm.sf(2.2-1.9, 0, np.sqrt(2.2+(0.8+1+1*4)/4))

for correct_uneq_calibs in ('ignore', 'mine', 'pathd8'):
    assert agesdf(tree, calib,
                  measures=['dist'],
                  unweighted=True, original_leading_paths=False,
                  todate=todate, keeproot=True,
                  correct_unequal_calibs=correct_uneq_calibs,
                  calib_selecter='name').age_dist['x'] == 2

ages = agesdf(tree, calib,
              measures=['dist'],
              unweighted=False, original_leading_paths=False,
              todate=todate, keeproot=True,
              correct_unequal_calibs='default',
              calib_selecter='name')
assert ages.age_dist['x'] > 2

tree = ete3.Tree('(((A:0.8,B:1)x0:1,C:2.2)y:0.8,D:3)R;', format=1)
ages = agesdf(tree, calib,
              measures=['dist'],
              unweighted=True, original_leading_paths=False,
              todate=todate, keeproot=True,
              correct_unequal_calibs='default',
              calib_selecter='name')
assert ages.cal['y'] == 1
assert ages.age_dist['y'] == 2
assert ages.cal['x0'] == 0
assert ages.age_dist['x0'] == 0.9

ages = agesdf(tree, calib,
              measures=['dist'],
              unweighted=False, original_leading_paths=False,
              todate=todate, keeproot=True,
              correct_unequal_calibs='default',
              calib_selecter='name')
assert ages.age_dist['y'] == 2
ages = agesdf(tree, calib,
              measures=['dist'],
              unweighted=True, original_leading_paths=True,
              todate=todate, keeproot=True,
              correct_unequal_calibs='default',
              calib_selecter='name')
assert ages.age_dist['y'] == 2

tree = ete3.Tree('(((A:0.8,B:1)y0:1,C:2.2)x:0.8,D:3)R;', format=1)
ages = agesdf(tree, calib,
              measures=['dist'],
              unweighted=True, original_leading_paths=False,
              todate=todate, keeproot=True,
              correct_unequal_calibs='default',
              calib_selecter='name')
assert ages.cal['y0'] == 1
assert ages.age_dist['y0'] == calib['y0']
assert ages.cal['x'] == 0

my0 = 0.9
mx_y0 = my0 + 1  # 1.9
mx_C = 2.2
mx = ((my0+1)*2 + 2.2)/3  # = 2
mR = ((mx+0.8)*3 + 3)/4   # = 2.85

# method 'mine':
ax = ( 2*(3 * (mx - my0)/mR*(3-0)/(3-1)) + 3*mx_C/mR) / 3  # = 1.9298
# method 'default':
ax = ( 2*(1 + (3 - 1)*(mx_y0-my0)/(mR - my0)) + 3*mx_C/mR ) / 3  # = 2.122357
assert round(ages.age_dist['x'], 4) == 2.1224
assert agesdf(tree, calib,
              measures=['dist'],
              unweighted=True, original_leading_paths=False,
              todate=todate, keeproot=True,
              correct_unequal_calibs='default',
              calib_selecter='name').age_dist.round(4)['x'] == 2.5  ## NOPE!

ages = agesdf(tree, calib,
              measures=['dist'],
              unweighted=False, original_leading_paths=False,
              todate=todate, keeproot=True,
              correct_unequal_calibs='default',
              calib_selecter='name')
mx = ((my0+1) + 2.2)/2  # = 2.05
mR = ((mx+0.8) + 3)/2   # = 2.85

# method 'mine': ax = ( 2*(1 + (mx - my0)*(3-0)/(3-1)) + 2.2) / 3  # = 2.5
ax = ( 1 + (3 - 1)*(mx_y0-my0)/(mR - my0) + 3*mx_C/mR ) / 2  # = 2.122357
assert round(ages.age_dist['x'], 4) == 2.1220

# Multiple measures at once.
tree = ete3.Tree('(((A:0.8[&&NHX:dS=0.01],B:1[&&NHX:dS=0.01])x0:1[&&NHX:dS=0.01],C:2.2[&&NHX:dS=0.02])x:1[&&NHX:dS=0.01],D:3[&&NHX:dS=0.03])R[&&NHX:S=100];', format=1)

for node in tree.iter_descendants():
    node.dS = float(node.dS)
tree.S = float(tree.S)

# Succeed all broadcasting operations?
ages = agesdf(tree, calib,
              measures=['dist', 'dS'],
              unweighted=False, original_leading_paths=False,
              todate=todate, keeproot=True,
              correct_unequal_calibs='default',
              calib_selecter='name')
ages = agesdf(tree, calib,
              measures=['dist', 'dS'],
              unweighted=True, original_leading_paths=False,
              todate=todate, keeproot=True,
              correct_unequal_calibs='default',
              calib_selecter='name')
ages = agesdf(tree, calib,
              measures=['dist', 'dS'],
              unweighted=False, original_leading_paths=True,
              todate=todate, keeproot=True,
              correct_unequal_calibs='default',
              calib_selecter='name')
ages = agesdf(tree, calib,
              measures=['dist', 'dS'],
              unweighted=True, original_leading_paths=True,
              todate=todate, keeproot=True,
              correct_unequal_calibs='default',
              calib_selecter='name')

# Compare method2 (scaling by keeping the original supporting branch length) VS method1
calib = {'R':2}
tree = ete3.Tree('((A:1,B:1)x:0.25,C:2.75)R;', format=1)
ages = agesdf(tree, calib,
              measures=['dist'],
              unweighted=False, original_leading_paths=False,
              todate=todate, keeproot=True,
              correct_unequal_calibs='default',
              calib_selecter='name')
assert ages.age_dist['x'] == 1.0

ages = agesdf(tree, calib,
              measures=['dist'],
              unweighted=False, original_leading_paths=True,
              todate=todate, keeproot=True,
              correct_unequal_calibs='default',
              calib_selecter='name')
assert ages.age_dist['x'] == 1.6  # 1/1.25 * 2


# Avoid negative branch lengths

tree = ete3.Tree('((A:1,B:1)x:0.25,C:0.3)R;', format=1)
ages = agesdf(tree, calib,
              measures=['dist'],
              unweighted=False, original_leading_paths=False,
              correct_unequal_calibs='default',
              todate=todate, keeproot=True,
              calib_selecter='name')
assert ages.age_dist['x'] == 2.0

tree = ete3.Tree('(((A:0.8,B:1)y0:1,C:0.2)x:0.2,D:0.3)R;', format=1)
ages = agesdf(tree, calib,
              measures=['dist'],
              unweighted=False, original_leading_paths=False,
              correct_unequal_calibs='default',
              todate=todate, keeproot=True,
              calib_selecter='name')
assert ages.age_dist['x'] == 1

# Numeric example of Britton 2002
tree = ete3.Tree('(Uvularia perfoliata:15,(Uvularia Pudica:8,Disporum:20)x0:3)R:6;', format=1)
## Check p-values of unclocklikeliness
ages = agesdf(tree, calib,
              measures=['dist'],
              unweighted=True, original_leading_paths=False,
              correct_unequal_calibs='default',
              todate=todate, keeproot=True,
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
              unweighted=False, original_leading_paths=False,
              correct_unequal_calibs='default',
              todate=todate, keeproot=True,
              calib_selecter='name')
# No difference for the unweighted algo if ≤3 leaves
