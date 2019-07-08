#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#%run -n ~/scripts/codeml/analyse/generate_dNdStable.py
from codeml.analyse.generate_dNdStable import *
import os.path as op
import pandas as pd
from types import GeneratorType
import matplotlib.pyplot as plt


logging.basicConfig(format=logging.BASIC_FORMAT)
logger.setLevel(logging.DEBUG)

# Numpy sugar
def a(*args, **kw):
    if len(args) == 1 and isinstance(args[0], (GeneratorType, range)):
        args = list(args[0])
    return np.array(args, **kw)

def c(*args, **kw):
    if len(args) == 1 and isinstance(args[0], (GeneratorType, range)):
        args = list(args[0])
    return np.concatenate(args, **kw)

m = np.nanmean


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
mR_y0 = 1.8

# method 'default':
ax = ( 2*(1 + (3 - 1)*(mx-my0)/(mR - my0)) + 3*mx/mR ) / 3  # = 2.12056  # = 2.122357
assert round(ages.age_dist['x'], 4) == 2.1206

# method 'mine':
ax = ( 2*((mx - my0)/(mR-my0)*(3-0)/(3-1)) + mx/mR*3) / 3  # = 1.9298

assert agesdf(tree, calib,
              measures=['dist'],
              unweighted=True, original_leading_paths=False,
              todate=todate, keeproot=True,
              correct_unequal_calibs='mine',
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


# Check that there is no *overweighting* with correct='default' formula.
tree = ete3.Tree('(((A:0.9,B:1.1)y0:0.7,C:2)x:1.2,D:3)R;', format=1)
ages = agesdf(tree, calib,
              measures=['dist'],
              unweighted=True, original_leading_paths=False,
              correct_unequal_calibs='default',
              todate=todate, keeproot=True,
              calib_selecter='name')
# We must set constraints:
# mR/aR == my0/ay0  ->  method 'ignore' *should* give the same result.
# my0 + dy0 != dC   ->  We will detect any overweighting of the y0 side.
# 1. 'ignore' method:
px = a(1.6, 1.8, 2)
mx = 1.8  # m(px)
pR = a(2.8, 3, 2.2, 3)  # c(px + 1.2, [3])
mR = 3.0  # m(pR)
ax = mx/mR * 3  # == mx

# 2. 'default' method  with constant rates, should always equal 'ignore'.
my0 = 1
mx = 1.8   # == ((my0 + 0.7)*2 + 2) / 3
mR = 3.0   # == ((my0 + 0.7 + 1.2)*2 + 3.2 + 3) / 4
# => MPL values are the same. What about ages?

# Paths going through y0 are shorter. If overweighting, this should attract the mean more.

# Method 'default':
# ax = (2 * (1 + (mx - my0)/(mR - my0)*(3-1)) + 3*mx/mR ) / 3
#    = (2 * (1 + (1.8 - 1)/(3-1)*(3-1)) + 3*1.8/3 ) / 3
assert ages.age_dist['x'] == 1.8
assert agesdf(tree, calib,
              measures=['dist'],
              unweighted=True, original_leading_paths=False,
              correct_unequal_calibs='ignore',
              todate=todate, keeproot=True,
              calib_selecter='name').age_dist['x'] == 1.8
# Check that there is equality in the weighted computation
assert agesdf(tree, calib,
              measures=['dist'],
              unweighted=False, original_leading_paths=False,
              correct_unequal_calibs='ignore',
              todate=todate, keeproot=True,
              calib_selecter='name'
              ).age_dist['x'] == agesdf(tree, calib,
                                      measures=['dist'],
                                      unweighted=False, original_leading_paths=False,
                                      correct_unequal_calibs='ignore',
                                      todate=todate, keeproot=True,
                                      calib_selecter='name').age_dist['x']
(tree&'y0').name = 'x0'
assert agesdf(tree, calib,
              measures=['dist'],
              unweighted=True, original_leading_paths=False,
              correct_unequal_calibs='default',
              todate=todate, keeproot=True,
              calib_selecter='name').age_dist['x'] == 1.8

assert agesdf(tree, calib,
              measures=['dist'],
              unweighted=False, original_leading_paths=False,
              correct_unequal_calibs='ignore',
              todate=todate, keeproot=True,
              calib_selecter='name'
              ).age_dist.round(4)['x'] == 1.8347


### Numeric example of Britton 2002

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

### Simulations

tree = ete3.Tree('(((A,B)y0,C)x,D)R;', format=1)
calib = {'R': 3, 'y': 2, 'y0': 1}

def simulate_brlen(tree, ages=None, rates=None, mean=None, stddev=None, law='normal'):
    # Assuming normal law for branch lengths with sd=length*rate
    nbranches = sum(1 for n in tree.iter_descendants())
    
    expected_subst = None
    if mean is None or stddev is None or law != 'normal':
        if ages is None:
            ages = {n.name: (0 if n.is_leaf()
                             else 1+n.get_farthest_leaf(topology_only=True)[1])
                    for n in tree.traverse()}
        if rates is None:
            rates = {}
        real_brlen = a(ages[n.up.name] - ages[n.name] for n in tree.iter_descendants())
        br_rates = a(rates.get(n.name, 1) for n in tree.iter_descendants())

        expected_subst = real_brlen * br_rates

    if mean is None:
        mean = expected_subst if law != 'poisson' else 1
    if stddev is None:
        stddev = np.sqrt(expected_subst)

    size = nbranches if (np.isscalar(mean) and np.isscalar(stddev)) else None

    if law == 'normal':
        sim_brlen = np.random.normal(mean, stddev, size=size)
        sim_brlen[sim_brlen<0] = 0
    elif law == 'lognormal':
        sim_brlen = np.random.lognormal(mean, stddev, size=size)
    elif law == 'poisson':
        sim_brlen = np.random.poisson(expected_subst * mean)
    else:
        raise ValueError('Invalid `law` parameter %r' % law)

    for n, brlen in zip(tree.iter_descendants(), sim_brlen):
        n.dist = brlen

    return tree


def date_tree(tree, ages, age_col='age_dist'):
    for n in tree.iter_descendants():
        n.dist = ages.loc[n.up.name, age_col] - ages.loc[n.name, age_col]
    return tree

sim_ages_x = []
for i in range(500):
    sim_ages_x.append(agesdf(simulate_brlen(tree), calib, measures=['dist'],
                             unweighted=False, original_leading_paths=False,
                             correct_unequal_calibs='default',
                             todate=todate, keeproot=True,
                             calib_selecter='name').age_dist['x'])
sim_ages_x = []
for i in range(500):
    sim_ages_x.append(agesdf(simulate_brlen(tree), calib, measures=['dist'],
                             unweighted=False, original_leading_paths=False,
                             correct_unequal_calibs='default',
                             fix_conflict_ages=False,
                             todate=todate, keeproot=True,
                             calib_selecter='name').age_dist['x'])

sim_ages_x = []
for i in range(500):
    sim_ages_x.append(agesdf(simulate_brlen(tree, stddev=1), calib, measures=['dist'],
                             unweighted=False, original_leading_paths=False,
                             correct_unequal_calibs='default',
                             fix_conflict_ages=True,
                             todate=todate, keeproot=True,
                             calib_selecter='name').age_dist['x'])

from datasci.graphs import plottree
from dendro.any import ete3 as ete3_methods
get_items = ete3_methods.get_items
get_label = ete3_methods.get_label

plottree(tree, get_items, get_label, style='V', alpha=0.05); ax = plt.gca()
for i in range(50):
    simtree = simulate_brlen(tree, stddev=0.5)
    maxlength = simtree.get_farthest_leaf()[1]
    for n in simtree.iter_descendants():
        n.dist *= 3./maxlength
    plottree(simtree, get_items, get_label, style='V', ax=ax, alpha=0.05)

plottree(tree, get_items, get_label, style='V', alpha=0.05); ax = plt.gca()
for i in range(50):
    simtree = simulate_brlen(tree, stddev=0.5)
    maxlength = simtree.get_farthest_leaf()[1]
    for n in simtree.iter_descendants():
        n.dist *= 3./maxlength
    plottree(simtree, get_items, get_label, style='V', ax=ax, alpha=0.05)


rates = {'A': 2, 'B': 2}

niter = 2000
alpha = 5./niter
sim_ages_x = np.full(niter, np.NaN)  # np.empty(niter)
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14,10))
ax2.plot([-2, -2], [-1, len(tree)], 'r--', label='Real age', alpha=0.6)
for i in range(niter):
    simtree = simulate_brlen(tree, rates=rates, mean=1, law='poisson')
    maxlength = simtree.get_farthest_leaf()[1]
    for n in simtree.iter_descendants():
        n.dist *= 3./maxlength
    plottree(simtree, get_items, get_label, style='V', ax=ax1, alpha=alpha)

    ages = agesdf(simtree, calib, measures=['dist'],
                                  unweighted=False, original_leading_paths=False,
                                  correct_unequal_calibs='default',
                                  fix_conflict_ages=True,
                                  todate=todate, keeproot=True,
                                  calib_selecter='name')
    sim_ages_x[i] = -ages.age_dist['x']
    dated_tree = date_tree(simtree, ages)
    plottree(dated_tree, get_items, get_label, style='V', ax=ax2, alpha=alpha)
ax1.set_xlabel('Nucleotide substitutions (simul)')
ax2.set_xlabel('Estimated age')
ax3 = ax2.twinx()
ax2.yaxis.set_label_position('right')
ax2.yaxis.set_ticks_position('right')
heights, _, _ = ax3.hist(sim_ages_x, bins=50, density=True, alpha=0.6)
ax3.set_ylim(0, heights.max()*(1 + len(tree)))
ax3.axis('off');

#fig.savefig(op.expanduser('~/ws2/DUPLI_data93/alignments_analysis/fig/simul_date-um1-default_normsd05.svg'))
fig.savefig(op.expanduser('~/ws2/DUPLI_data93/alignments_analysis/fig/simul_date-um1-default_lognorm.svg'))
fig.savefig(op.expanduser('~/ws2/DUPLI_data93/alignments_analysis/fig/simul_date-um1-default_lognormsd05.svg'))
fig.savefig(op.expanduser('~/ws2/DUPLI_data93/alignments_analysis/fig/simul_date-um1-default_poisson10.svg'))
fig.savefig(op.expanduser('~/ws2/DUPLI_data93/alignments_analysis/fig/simul_date-um1-default_poisson1.svg'))

fig.savefig(op.expanduser('~/ws2/DUPLI_data93/alignments_analysis/fig/simul_date-um1-default_rate2_normsd05.svg'))
fig.savefig(op.expanduser('~/ws2/DUPLI_data93/alignments_analysis/fig/simul_date-um1-default_rate2_poisson10.svg'))
fig.savefig(op.expanduser('~/ws2/DUPLI_data93/alignments_analysis/fig/simul_date-um1-default_rate2_poisson1.svg'))
