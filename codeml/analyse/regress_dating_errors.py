#!/usr/bin/python3
# coding: utf-8

# Source:
# jupyter nbconvert --to python \
#   ~/ws2/DUPLI_data85/alignments_analysis/subtrees_stats/subtrees_stats.ipynb

import warnings
from io import StringIO
from collections import OrderedDict
import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.use('TkAgg', warn=False)
get_ipython().magic('matplotlib inline')
import matplotlib.pyplot as plt
import seaborn as sb

import os.path as op

from codeml.analyse.dSvisualizor import splitname2taxongenetree
from seqtools.compo_freq import weighted_std
from plottools import scatter_density
from dendron.climber import dfw_pairs_generalized

mpl.style.use("softer")
pd.set_option("display.max_columns", 50)
pd.set_option("display.width", 115)
pd.set_option("display.max_colwidth", 50)
pd.set_option("display.show_dimensions", True)  # even on non truncated dataframes

mpl.rcParams['figure.figsize'] = (22, 14) # width, height

from scipy import stats
import scipy.cluster.hierarchy as hclust
import scipy.spatial.distance as spdist
#stats.skew, stats.kurtosis

from LibsDyogen import myPhylTree

from sklearn.decomposition import PCA
from sklearn.linear_model import LinearRegression
from sklearn.preprocessing import StandardScaler
import statsmodels.api as sm
#import statsmodels.formula.api as smf

from IPython.display import display_html
import logging
logger=logging.getLogger(__name__)
logging.basicConfig()
logger.setLevel(logging.INFO)


# Convert "time used" into seconds.
def time2seconds(time_str):
    factors = [1, 60, 3600, 3600*24]
    s = 0
    for factor, n_units in zip(factors, reversed(time_str.split(':'))):
        s += factor * int(n_units)
    return s


def load_subtree_stats(template, stattypes=('al', 'tree', 'codeml')):
    """
    param: template to the csv/tsv files containing al/tree/codeml stats.

    Example template: 'subtreesRawTreeBestO2_{stattype}stats_Simiiformes.tsv'
    """
    
    # ## Load tree/alignment statistics
    alfile = template.format(stattype='al')
    treefile = template.format(stattype='tree')
    codemlfile = template.format(stattype='codeml')

    print('Load', alfile)
    aS = pd.read_table(alfile, index_col=0)
    print('Load', treefile)
    ts = pd.read_table(treefile, index_col=0,
                       dtype={'leaves_robust':      bool,
                              'single_child_nodes': bool,
                              'nodes_robust':       bool,
                              'only_treebest_spe':  int,
                              'aberrant_dists':     int})
                              #'root2tip_var': float
    ts['really_robust'] = ts.leaves_robust & ts.nodes_robust & ~ts.only_treebest_spe

    # ## Load codeml output statistics
    print('Load', codemlfile)
    cs = pd.read_table(codemlfile, index_col=0)
    cs['seconds'] = cs['time used'].apply(time2seconds)

    return aS, ts, cs

def check_load_subtree_stats(aS, ts, cs):
    print("shapes: aS %s, ts %s, cs %s" % (aS.shape, ts.shape, cs.shape))
    print("aS has dup:", aS.index.has_duplicates)
    print("ts has dup:", ts.index.has_duplicates)
    print("cs has dup:", cs.index.has_duplicates)
    common_subtrees = set(aS.index) & set(ts.index) & set(cs.index)
    print("%d common subtrees" % len(common_subtrees))
    only_al = aS.index.difference(ts.index.union(cs.index))
    only_tr = ts.index.difference(aS.index.union(cs.index))
    only_co = cs.index.difference(aS.index.union(ts.index))
    l_al = len(only_al)
    l_tr = len(only_tr)
    l_co = len(only_co)
    print("%d only in al stats: %s" % (l_al, list(only_al)[:min(5, l_al)]))
    print("%d only in tree stats: %s" % (l_tr, list(only_tr)[:min(5, l_tr)]))
    print("%d only in codeml stats: %s" % (l_co, list(only_co)[:min(5, l_co)]))
    # Todo: pyupset plot



# ## Function to merge additional subgenetree information into `ages`
def merge_criterion_in_ages(criterion_serie, ages=None, ages_file=None,
                            criterion_name=None):
    """Merge a column into the *node ages* table: the common field is the *subgenetree*."""
    
    assert (ages is not None or ages_file) and not (ages_file and ages is not None), "At least `ages` (dataframe) or `ages_file` (filename) must be given."
    if ages is None:
        ages = pd.read_table(ages_file, sep='\t')
    
    print("Input shape:", ages.shape)
    criterion = criterion_serie.name if not criterion_name else criterion_name
    ages_subgenetrees = ages.subgenetree.unique()
    
    assert len(criterion_serie) >= ages_subgenetrees.size, \
            "Not all genetrees have a criterion value: %d < %d" % (len(criterion_serie), ages_subgenetrees.size)
    assert set(criterion_serie.index) >= set(ages_subgenetrees)
    
    #criterion_df = pd.DataFrame({criterion: criterion_serie})
    criterion_df = criterion_serie.to_frame(criterion)
    assert (criterion_df.index == criterion_serie.index).all()
    
    ages_c = pd.merge(ages, criterion_df, how='left', left_on='subgenetree', right_index=True)
    # If all subgenetrees have a criterion, equivalent to:
    #ages[criterion] = criterion[ages.subgenetree]
    return ages_c


def add_robust_info(ages_p, ts):
    # Compute Number of duplications/speciation per tree.
    Ndup = ages_p.groupby('subgenetree').type.agg(lambda v: sum(v == "dup"))
    Ndup.name = 'Ndup'
    print(Ndup.describe())
    Nspe = ages_p.groupby('subgenetree').type.agg(lambda v: sum(v == "spe"))
    Nspe.name = 'Nspe'
    print(Nspe.describe())
    # merge tree stats to select robust trees
    ages_treestats = pd.merge(ages_p.drop('_merge', axis=1),
                              pd.concat((Ndup, Nspe,
                                         ts[['really_robust','aberrant_dists']]),
                                        axis=1,
                                        join='outer', sort=False),
                              how='left', left_on='subgenetree', right_index=True,
                              indicator=True, validate='many_to_one')
    print("Ndup", Ndup.shape, "Nspe", Nspe.shape, "ages_treestats", ages_treestats.shape)
    print(ages_treestats.groupby('_merge')['_merge'].count())
    return ages_treestats, Ndup, Nspe
    

def load_prepare_ages(ages_file, ts):
    ages = pd.read_table(ages_file, sep='\t', index_col=0)

    print("Shape ages: %s; has dup: %s" % (ages.shape, ages.index.has_duplicates))
    n_nodes = ages.shape[0]
    print("Shape ages internal nodes:", ages[ages.type != 'leaf'].shape)
    n_nodes_int = (ages.type != 'leaf').sum()

    # Fetch parent node info
    ages_p = pd.merge(ages,
                      ages[ages.type != 'leaf'][['taxon', 'type', 'age_t',
                              'age_dS', 'age_dN', 'age_dist', 'calibrated']],
                      how="left", left_on="parent", right_index=True,
                      suffixes=('', '_parent'), indicator=True,
                      validate='many_to_one')
    print("Shape ages with parent info: %s" % (ages_p.shape,))
    n_nodes_p = ages_p.shape[0]

    if n_nodes_p < n_nodes:
        logger.warning('%d nodes were lost when fetching parent information.',
                       n_nodes - n_nodes_p)

    orphans = ages_p._merge=='left_only'
    ages_orphans = ages_p[orphans]

    print("\nOrphans: %d\n" % ((orphans).sum()))
    #display_html(ages_orphans.head())

    assert ((ages_orphans.parent == ages_orphans.root) |
            (ages_orphans.type == 'leaf')).all()
    print("All orphans are expected (leaf, or child of the root).")
    
    e_nochild = set(ages.index) - set(ages.parent)  # expected
    parent_nodata = set(ages[ages.type!='leaf'].parent) - set(ages_p.index)
    n_nochild = len(e_nochild)
    print("\nExpected nodes without children (leaves): %d" % (n_nochild,))
    print("Observed nodes not found as parents: %d" % \
            len(set(ages.index) - set(ages_p[ages_p._merge=='both'].index)))
    print("Parent nodes without data: %d" % (len(parent_nodata),))
    #assert len(nochild) == n_nochild, \
    #    "Found %d unexpected nodes without child:\n%s" % \
    #        (len(nochild) - n_nochild,
    #         ages.loc[nochild - e_nochild])

    ages_spe2spe = ages_p[(ages_p.type.isin(('spe', 'leaf'))) & \
                          (ages_p.type_parent == 'spe')]
    print("\nShape ages speciation to speciation branches (no dup):", ages_spe2spe.shape)
    ages_treestats, Ndup, Nspe = add_robust_info(ages_p, ts)
    ages_robust = ages_treestats[ages_treestats.really_robust & \
                                 (ages_treestats.aberrant_dists == 0)]\
                            .drop(['really_robust', 'aberrant_dists'], axis=1)
    print("\n%d nodes from robust trees" % (ages_robust.shape[0],))
    return ages_treestats, ages_robust, Ndup, Nspe


# for averaging by taking into account branch length: with Omega.
# NOTE: columns will be reordered following `var`.
def group_average(g, var, weight_var="median_brlen"):
    return pd.Series(np.average(g[var], axis=0, weights=g[weight_var]))

def group_weighted_std(g, var, weight_var="median_brlen"):
    return pd.Series(weighted_std(g[var], axis=0, weights=g[weight_var]))

# for calculating rates: with dS, dN, t
def tree_dist_2_rate(g, dist_var, norm_var="median_brlen"):
    # in pandas, sum is done per columns.
    return g[dist_var].sum() / g[norm_var].sum()



#def add_control_dates_lengths(ages_robust, phyltree, timetree_ages_CI=None):
#    pass
#    return ages_controled, median_taxon_ages, median_brlen

def add_control_dates_lengths(ages, ages_robust, phyltree, timetree_ages_CI=None):
    # Merge control dates
    median_taxon_ages = ages[ages.type.isin(("spe", "leaf"))]\
                               .groupby("taxon").age_dS.median()
                               #& (ages_robust.taxon != 'Simiiformes')]\
    median_taxon_ages.name = 'median_taxon_age'
    
    timetree_ages = median_taxon_ages.index.to_series().apply(phyltree.ages.get)
    timetree_ages.name = 'timetree_age'

    control_ages = pd.concat((median_taxon_ages, timetree_ages,
                              timetree_ages_CI), axis=1, sort=False)

    print(control_ages.sort_values('timetree_age', ascending=False))

    ages_controled = pd.merge(ages_robust, control_ages,
                              left_on="taxon", right_index=True, validate="many_to_one")
    ages_controled_withnonrobust = pd.merge(ages, control_ages,
                              left_on="taxon", right_index=True, validate="many_to_one")

    # Merge control branch lengths
    invalid_taxon_parent = ages_controled.taxon_parent.isna()
    # Should be nodes whose parent node is the root.
    if invalid_taxon_parent.any():
        assert (ages_controled[invalid_taxon_parent].parent == \
                ages_controled[invalid_taxon_parent].root).all()
        debug_columns = ['parent', 'subgenetree', 'taxon', 'taxon_parent', 'median_taxon_age']
        logger.error("%d invalid 'taxon_parent':\n%s\n"
                     "The following taxa have no parent taxa information, "
                     "please check:\n%s\n**DROPPING** this data!",
                       invalid_taxon_parent.sum(),
                       ages_controled[invalid_taxon_parent][
                           debug_columns
                           ].head(),
                       ', '.join(ages_controled[invalid_taxon_parent].taxon.unique()))
        ages_controled.dropna(subset=['taxon_parent'], inplace=True)
    
    ages_controled['median_brlen'] = \
        ages_controled.taxon_parent.apply(control_ages.median_taxon_age.get) \
        - ages_controled.median_taxon_age

    ages_controled['timetree_brlen'] = \
        ages_controled.taxon_parent.apply(control_ages.timetree_age.get) \
        - ages_controled.timetree_age

    # Resulting branch lengths
    control_brlen = ages_controled[
                        ~ages_controled.duplicated(["taxon_parent", "taxon"])
                        ][
                           ["taxon_parent", "taxon", "median_brlen",
                            "median_taxon_age", "timetree_brlen", "timetree_age"]
                        ].sort_values("taxon_parent", ascending=False)

    branch_info = ["taxon_parent", "taxon"]
    control_brlen.index = pd.MultiIndex.from_arrays(
                                            control_brlen[branch_info].values.T,
                                            names=branch_info)
    control_brlen.drop(branch_info, axis=1, inplace=True)
    display_html(control_brlen)
    return ages_controled, ages_controled_withnonrobust, control_ages, control_brlen


def check_control_dates_lengths(control_brlen, phyltree, root):
    get_phylchildren = lambda phyltree, ancdist: phyltree.items.get(ancdist[0], [])

    expected_branches, expected_dists = zip(*(((p[0],ch[0]),ch[1]) for p,ch in
                             dfw_pairs_generalized(phyltree,
                                                   get_phylchildren,
                                                   queue=[(None, (root,0))])))

    #logger.debug(expected_branches)
    #logger.debug(expected_dists)

    median_brlen_sum = control_brlen.median_brlen.sum()
    print("Sum of median branch lengths =", median_brlen_sum, "My")
    timetree_brlen_sum = control_brlen.timetree_brlen.sum()
    print("Sum of timetree branch lengths =", timetree_brlen_sum, "My")
    real_timetree_brlen_sum = sum(expected_dists)
    print("Real sum of TimeTree branch lengths (in phyltree) =",
          real_timetree_brlen_sum, "My")

    unexpected_branches = set(control_brlen.index) - set(expected_branches)
    if unexpected_branches:
        logger.error("Extraneous branches not seen in phyltree:\n%s",
                     unexpected_branches)
    lost_branches = set(expected_branches) - set(control_brlen.index)
    if lost_branches:
        logger.error("Forgotten branches in phyltree:\n%s",
                     lost_branches)
    
    median_treelen_phyltree = control_brlen.reindex(list(expected_branches)).median_brlen.sum()
    timetree_treelen_phyltree = control_brlen.reindex(list(expected_branches)).timetree_brlen.sum()
    print("Sum of median branch lengths for branches found in phyltree =",
          median_treelen_phyltree)
    print("Sum of timetree branch lengths for branches found in phyltree =",
          timetree_treelen_phyltree)


def compute_dating_errors(ages_controled, control='median'):
    if control == 'timetree':
        raise NotImplementedError("control ages/brlen with timetree.")

    ages_controled["abs_age_error"] = \
                (ages_controled.age_dS - ages_controled.median_taxon_age).abs()
    ages_controled["signed_age_error"] = \
                (ages_controled.age_dS - ages_controled.median_taxon_age)
    ages_controled["abs_brlen_error"] = \
                (ages_controled.age_dS_parent - ages_controled.age_dS - \
                 ages_controled.median_brlen).abs()
    ages_controled["signed_brlen_error"] = \
                (ages_controled.age_dS_parent - ages_controled.age_dS - \
                 ages_controled.median_brlen)

    mean_errors = ages_controled[
                     ['subgenetree', 'abs_age_error', 'signed_age_error',
                      'abs_brlen_error', 'signed_brlen_error']
                    ].groupby("subgenetree").mean()
    return mean_errors


def compute_branchrate_std(ages_controled, dist_measures):
    
    groupby_cols = ["subgenetree", "median_brlen", "median_taxon_age",
                    "taxon_parent", "taxon"] + dist_measures
    #ages_controled["omega"] = ages_controled.branch_dN / ages_controled.branch_dS

    sgg = subgenetree_groups = ages_controled[groupby_cols].groupby('subgenetree')

    # ### Average (substitution) rates over the tree:
    #     sum of all branch values / sum of branch lengths

    # Sum aggregation + division broadcasted on columns
    cs_rates = sgg[dist_measures].sum().div(sgg.median_brlen.sum(), axis=0)
    #cs_rates["omega"] = (sgg.branch_dN / sgg.branch_dS).apply() 
    rate_measures = [(m.replace('branch_', '') + '_rate') for m in dist_measures]
    cs_rates.columns = rate_measures

    # ### Weighted standard deviation of substitution rates among branches

    tmp = pd.merge(ages_controled[["subgenetree", "median_brlen"] + dist_measures],
                   cs_rates, left_on="subgenetree", right_index=True)

    #rate_dev = pd.DataFrame({
    #            "branch_dist": (tmp.branch_dist/ tmp.median_brlen - tmp.dist_rate)**2,
    #            "branch_t":    (tmp.branch_t   / tmp.median_brlen - tmp.t_rate)**2,
    #            "branch_dS":   (tmp.branch_dS  / tmp.median_brlen - tmp.dS_rate)**2,
    #            "branch_dN":   (tmp.branch_dN  / tmp.median_brlen - tmp.dN_rate)**2,
    #            #"omega":       (ages_controled.dN / ages_controled.dS - tmp.omega)**2,
    #            "subgenetree": tmp.subgenetree,
    #            "median_brlen": tmp.median_brlen})

    # subtract branch rate with mean rate, then square.
    rate_dev_dict = {d: (tmp[d] / tmp.median_brlen - tmp[r])**2
                     for d,r in zip(dist_measures, rate_measures)}
    rate_dev_dict.update(subgenetree=tmp.subgenetree,
                         median_brlen=tmp.median_brlen)
    rate_dev = pd.DataFrame(rate_dev_dict)

    cs_wstds = rate_dev.groupby("subgenetree").apply(
                (lambda x, var, weight_var:
                                    sqrt(group_average(x, var, weight_var))),
                dist_measures, "median_brlen")
                
    cs_wstds.columns = [(r + '_std') for r in cs_rates.columns]

    return cs_rates, cs_wstds


def subset_on_criterion_tails(criterion_serie, ages=None, ages_file=None,
                              outbase=None, criterion_name=None, nquantiles=4,
                              thresholds=None, save=False):
    """From the input data, output two files:
    - one with the lower quantile of criterion values,
    - one with the upper quantile.
    
    Otherwise thresholds can be used as a tuple (lower_tail_max, upper_tail_min)
    
    output files will be named as `outbase + ("-lowQ%d" % nquantiles) + ".tsv"`
    """
    ages_c = merge_criterion_in_ages(criterion_serie, ages, ages_file, criterion_name)
        
    print(ages.columns)
    
    if not outbase and ages_file:
        outbase, _ = os.path.splitext(ages_file)
    elif not outbase and not ages_file:
        outbase = "dataset"
    
    if not thresholds:
        q = 1 / nquantiles
        low_lim = criterion_serie.quantile(q)
        high_lim = criterion_serie.quantile(1. - q)
        outlow = "%s-%slowQ%d.tsv" % (outbase, criterion_name, nquantiles)
        outhigh = "%s-%shighQ%d.tsv" % (outbase, criterion_name, nquantiles)
    else:
        low_lim, high_lim = thresholds
        outlow = outbase + "-%slow%1.3f.tsv" % (criterion_name, low_lim) 
        outhigh = outbase + "-%shigh%1.3f.tsv" % (criterion_name, high_lim)
        
    print("Output %s values outside %s (from %d quantiles)" % (criterion_name,
                                    [low_lim, high_lim],
                                    (nquantiles if not thresholds else None)))
    if save:
        ages_c[ages_c[criterion_name] <= low_lim].to_csv(outlow, sep='\t')
        ages_c[ages_c[criterion_name] >= high_lim].to_csv(outhigh, sep='\t')
    print("Output files:", outlow, outhigh, sep='\n')


def get_tails_on_criterion(df, criterion_name, nquantiles=4):
    q = 1. / nquantiles
    low_lim, high_lim = df[criterion_name].quantile([q, 1. - q])
    df_low = df[df[criterion_name] <= low_lim].copy()
    df_high = df[df[criterion_name] >= high_lim].copy()
    print(low_lim, high_lim)
    return df_low, df_high


def annot_quantiles_on_criterion(ages, criterion_serie, criterion_name=None,
                                 nquantiles=4, transform=None):
    assert criterion_name not in ages.columns
    
    # Exclude the points 0. and 1.
    q = np.linspace(1./nquantiles, 1, nquantiles-1, endpoint=False)
    
    if transform is not None:
        criterion_serie = transform(criterion_serie)
        criterion_name = transform.__name__ + '_' + criterion_name

    quantiles = criterion_serie.quantile(q)
    print("Quantiles (n=%d) values for %s:" % (nquantiles, criterion_name),
          quantiles, sep='\n')

    ages_c = merge_criterion_in_ages(criterion_serie, ages, criterion_name=criterion_name)
    ages_c["Q_" + criterion_name] = ages_c[criterion_name].apply(
                                lambda v: np.searchsorted(quantiles.values, v))
    return ages_c


def _violin_spe_ages_vs_criterion_quantiles(annot_df, criterion_name, isin=None,
                                           split=True):
    Q_col = "Q_" + criterion_name
    if isin is None:
        # Look at extreme quantiles only
        isin = (annot_df[Q_col].min(), annot_df[Q_col].max())
        print(isin)
        
    ax = sb.violinplot(x="taxon", y="age_dS", hue=Q_col, data=annot_df[(annot_df.type == "spe") & annot_df[Q_col].isin(isin)], split=split)
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45)
    return ax


def violin_spe_ages_vs_criterion(ages, criterion_serie, criterion_name=None,
                                 nquantiles=10, split=False):
    criterion_name = criterion_name or criterion_serie.name
    annot_ages = annot_quantiles_on_criterion(ages, criterion_serie,
                                              criterion_name, nquantiles)
    isin = None if split else list(range(nquantiles))
    ax = _violin_spe_ages_vs_criterion_quantiles(
                                annot_ages[annot_ages.taxon != 'Simiiformes'],
                                criterion_name, isin, split)


# Functions to check if variables need transformation

def normal_fit(var):
    """Return normal density curve with same mean and variance as the given distribution"""
    sorted_var = var.sort_values()
    return sorted_var, stats.norm.pdf(sorted_var, sorted_var.mean(), sorted_var.std()) # * np.isfinite(sorted_var).sum()


def all_test_transforms(alls, variables):
    #fig, axes = plt.subplots(len(variables),3, figsize=(22, 5*len(variables)))
    nbins = 50

    for i, ft in enumerate(variables):
        
        var = alls[ft]
        
        # Plot original variable distribution 
        if var.dtype != float:
            print("Variable %r not continuous: %s" % (ft, var.dtype))
            if var.dtype != int:
                print("Variable %r does not seem to be numeric. Skipping" % ft)
                continue
        
        fig, axes = plt.subplots(1, 3, figsize=(22, 5))
        #axes[i, 0].set_ylabel(ft)
        axes[0].set_ylabel(ft)

        axes[0].hist(var, bins=nbins, density=True)
        axes[0].plot(*normal_fit(var), '-')
        axes[0].set_title("original")
        _, xmax0 = axes[0].get_xlim()
        _, ymax0 = axes[0].get_ylim()
        
        text = "Skew: %g\nKurtosis: %g\n" % (var.skew(), var.kurt())
        if (var < 0).any():
            if (var > 0).any():
                print("Variable %r has negative and positive values. Shifting to positive." % ft)
                text += "Negative and positive values. Shifting to positive.\n"
                var -= var.min()
            else:
                print("Variable %r converted to positive values" % ft)
                text += "Converted to positive values.\n"
                var = -var
        
        # Plot log-transformed distribution
        with warnings.catch_warnings(record=True) as w:
            logtransformed_var = np.log10(var)
            if w:
                assert issubclass(w[-1].category, RuntimeWarning)
                assert "divide by zero encountered in log10" in str(w[-1].message)
                zero_vals = True
                n_zero_vals = (var == 0).sum()

        n_infinite_vals = (~np.isfinite(logtransformed_var)).sum()
        logtransformed_var = logtransformed_var[np.isfinite(logtransformed_var)]
        
        logskew, logkurt = logtransformed_var.skew(), logtransformed_var.kurt()
        logtext = "Skew: %g\nKurtosis: %g\n" % (logskew, logkurt)

        axes[1].hist(logtransformed_var, bins=nbins, density=True, alpha=0.5)
        axes[1].plot(*normal_fit(logtransformed_var), '-', alpha=0.5)
        
        if n_infinite_vals:
            suggested_increment = logtransformed_var.quantile(0.05)
            print("%s: Nb of not finite values: %d. Suggested increment: 10^(%g)"\
                        % (ft, n_infinite_vals, suggested_increment))
            text += "%d not finite values. Suggested increment: 10^(%g)"\
                        % (n_infinite_vals, suggested_increment)
            
            logtransformed_inc_var = np.log10(var + 10**suggested_increment) #1.) 
            twin_ax1 = axes[1].twinx()
            twin_ax1.hist(logtransformed_inc_var, bins=nbins, density=True,
                          alpha=0.5, label="Incremented", color="#78a86a") # Green
            twin_ax1.plot(*normal_fit(logtransformed_inc_var), '-', alpha=0.5)
            
            logtransformed_inc1_var = np.log10(var + 1)
            twin2_ax1 = axes[1].twinx()
            twin2_ax1.hist(logtransformed_inc1_var, bins=nbins, density=True,
                           alpha=0.5, label="Increment+1", color="#a86a78") # Red
            twin2_ax1.plot(*normal_fit(logtransformed_inc1_var), '-', alpha=0.5)
            
            logtext_inc = ("Skew     (+inc): %g\nKurtosis (+inc): %g\n"
                           % (logtransformed_inc_var.skew(),
                              logtransformed_inc_var.kurt()))
            logtext_inc1 = ("Skew     (+1): %g\nKurtosis (+1): %g\n"
                            % (logtransformed_inc1_var.skew(),
                               logtransformed_inc1_var.kurt()))

        xmin1, xmax1 = axes[1].get_xlim()
        _, ymax1 = axes[1].get_ylim()

        axes[1].set_title("log10 transformed")
        
        sqrttransformed_var = np.sqrt(var)
        sqrttext = "Skew: %g\nKurtosis: %g\n" % (sqrttransformed_var.skew(),
                                                 sqrttransformed_var.kurt())
        axes[2].hist(sqrttransformed_var, bins=nbins, density=True)
        axes[2].plot(*normal_fit(sqrttransformed_var), '-')
        axes[2].set_title("Square root transformed")
        _, xmax2 = axes[2].get_xlim()
        _, ymax2 = axes[2].get_ylim()
        
        axes[0].text(xmax0, ymax0, text, va="top", ha="right")

        xpos1 = xmax1 if logskew>0 else xmin1
        ha1 = 'right' if logskew>0 else 'left'
        if n_infinite_vals:
            axes[1].text(xpos1, ymax1, logtext_inc, va="top", ha=ha1, color='#78a86a')
            axes[1].text(xpos1, 0.85*ymax1, logtext_inc1, va="top", ha=ha1, color='#a86a78')
        else:
            axes[1].text(xpos1, ymax1, logtext, va="top", ha=ha1)

        axes[2].text(xmax2, ymax2, sqrttext, va="top", ha="right")

        #fig.show(warn=False)
        plt.show()


# Variable transformation

notransform = lambda x: x
#notransform.__name__ = "%s"
sqrt = np.sqrt
log = np.log10

logneg = lambda x: np.log10(-x)
logneg.__name__ = "log10(-%s)"

def make_logtransform_inc(inc=0):
    loginc = lambda x: np.log10(x + inc)
    loginc.__name__ = "log10(%g+%%s)" % inc
    return loginc

def make_logpostransform_inc(inc=0):
    loginc = lambda x: np.log10(x + x.min() + inc)
    loginc.__name__ = "log10(%g+min+%%s)" % (inc)
    return loginc


# Functions for the PCA

def plot_cov(ft_cov, features, cmap='seismic', figax=None, cax=None):
    cmap = plt.get_cmap(cmap)
    norm = mpl.colors.Normalize(-1, 1)
    fig, ax = plt.subplots() if figax is None else figax
    img = ax.imshow(ft_cov, cmap=cmap, norm=norm, aspect='auto', origin='lower') #plt.pcolormesh
    ax.set_xticks(np.arange(len(features)))
    ax.set_yticks(np.arange(len(features)))
    ax.set_yticklabels(features)
    ax.set_xticklabels(features, rotation=45, ha='right')
    ax.set_ylabel("Features")
    ax.set_title("Feature covariance")
    fig.colorbar(img, ax=ax, cax=cax, aspect=ft_cov.shape[0])

def cov2cor(cov):
    """Converts covariance matrix into correlation matrix"""
    var = np.diagonal(cov)[:,np.newaxis]
    return cov / np.sqrt(var.dot(var.T))


def heatmap_cov(ft_cov, features, cmap='seismic',
                dendro_ratio=0.20, dendro_pad=0.1, cb_ratio=0.025, cb_pad=0.025):
    """plot_cov, but with hierarchical clustering on the side"""
    # Tested with figsize=(20, 12)
    fig, ax = plt.subplots()
    (x0, y0), (w, h) = ax.get_position().get_points()
    # absolute padding (in figure coordinate)
    # correct ratio by taking pad into account
    
    # position relatively to figure (percentages)
    ax.set_position([x0 + (dendro_ratio+dendro_pad)*w, y0,
                     w*(1-dendro_ratio-dendro_pad-cb_ratio-cb_pad), h])
    #width2 = width*ratio - float(pad)/w
    ax_ddg = fig.add_axes([x0, y0, w*dendro_ratio, h], frameon=False)
    #ax_ddg.set_title("hierarchical clustering (euclidean)")
    #ax_ddg.axis('off')
    ax_ddg.xaxis.set_visible(False)
    ax_cb = fig.add_axes([x0 + w*(1-cb_ratio), y0, w*cb_ratio, h])

    distmat = 1 - np.abs(cov2cor(ft_cov))
    tol=1e-15
    #assert (np.diag(distmat) < tol).all()
    #assert (np.abs(distmat - distmat.T) < tol).all()
    spdist.is_valid_dm(distmat, tol, throw=True)

    flatdist = spdist.squareform(distmat, force='tovector', checks=False)
    Z = hclust.linkage(flatdist, method='average', metric='euclidean')
    ddg = hclust.dendrogram(Z, orientation='left', no_labels=True, #labels=features,
                            ax=ax_ddg)

    clustered_ft_cov = ft_cov[ddg['leaves'],:][:,ddg['leaves']]
    #print(ddg['leaves'], ft_cov.shape)
    #print(clustered_ft_cov)
    logger.debug(np.array(features)[ddg['leaves']])
    plot_cov(clustered_ft_cov,
             np.array(features)[ddg['leaves']], cmap, (fig, ax), ax_cb)
    plt.show()


def centered_background_gradient(s, cmap='PRGn', center=0, extend=0):
    """Color background in a range according to the data, centered around the given value.
    Adapted from pandas.io.formats.style.Styler.background_gradient()"""
    smin, smax = s.min(), s.max()
    assert smin <= center and center <= smax
    most_distant_absval = max(center - smin, smax - center)
    rng = 2 * most_distant_absval
    # extend lower / upper bounds, compresses color range
    norm = mpl.colors.Normalize(center - most_distant_absval - (rng * extend),
                                center + most_distant_absval + (rng * extend))
    # matplotlib modifies inplace?
    # https://github.com/matplotlib/matplotlib/issues/5427
    normed = norm(s.values)
    c = [mpl.colors.rgb2hex(x) for x in plt.cm.get_cmap(cmap)(normed)]
    return ['background-color: {color}'.format(color=color) for color in c]


def magnify():
    """Increase the size of the table cells."""
    return [dict(selector="th",
                 props=[("font-size", "12pt")]),
            dict(selector="td",
                 props=[('padding', "0em 0em")]),
            dict(selector="th:hover",
                 props=[("font-size", "12pt")]),
            dict(selector="tr:hover td:hover",
                 props=[('max-width', '200px'),
                        ('font-size', '12pt')])
            ]


def plot_loadings(components, cmap="PRGn"):
    """Not as great as the pandas dataframe styling, because the scale is on
    the entire image here."""
    cmap = plt.get_cmap(cmap)
    norm = mpl.colors.Normalize(components.min().min(), components.max().max())
    plt.imshow(components, cmap=cmap, norm=norm)
    ax = plt.gca()
    ax.set_xticks(np.arange(components.shape[1]))
    ax.set_yticks(np.arange(components.shape[0]))
    ax.set_yticklabels(components.index)
    ax.set_xticklabels(components.columns)
    plt.colorbar()
    ax.set_title("Feature contribution")


def plot_features_PCspace(components, features, PCs=["PC1", "PC2"], ax=None):
    quiver = plt.quiver if ax is None else ax.quiver 
    quiver(0, 0, components[PCs[0]], components[PCs[1]], units='dots', width=1,
          scale_units='width')
    ax = ax or plt.gca()
    
    for ft in features:
        ft_vect = components.loc[ft][PCs] * 0.1
        angle = ((np.arctan2(ft_vect[1], ft_vect[0]) / (2*np.pi)) * 360)
        ax.annotate(s=ft, xy=ft_vect, xytext=1.25*ft_vect, rotation=angle)
        #plt.text(ft_vect[0], ft_vect[1], ft)
    
    ax.set_xlabel(PCs[0])
    ax.set_ylabel(PCs[1])
    return ax


def detailed_pca(alls_normed, features):

    ft_pca = PCA(n_components=15)
    ft_pca_components = ft_pca.fit_transform(alls_normed[features])

    # ### PC contribution to variance

    print("Components dimensions: %s" % (ft_pca_components.shape,))

    # Explained variance of each principal component.
    PC_contrib = ft_pca.explained_variance_ratio_
    print("Feature contributions:\n", PC_contrib)

    # Plot cumulative contribution of PC
    fig, ax = plt.subplots()
    ax.bar(np.arange(PC_contrib.size), PC_contrib.cumsum())
    ax.set_title("Cumulative ratio of variance explained by the Principal Components")
    ax.set_xlabel("Principal Components")
    ax.set_ylabel("Cumulative ratio of variance explained");
    plt.show()

    # Coefficients of the linear combination of each parameter in the resulting components
    print("Components dimensions:", ft_pca.components_.shape)

    # ### PC loadings

    # weights of each feature in the PCs

    print("### PC loadings")
    print("%-17s:\t%10s\t%10s" % ('Feature', 'coef PC1', 'coef PC2'))
    for ft, coef1, coef2 in sorted(zip(features, ft_pca.components_[0,:],
                                       ft_pca.components_[1,:]),
                                   key=lambda x: (abs(x[1]), abs(x[2])),
                                   reverse=True):
        print("%-17s:\t%10.6f\t%10.6f" % (ft, coef1, coef2))

    components = pd.DataFrame(ft_pca.components_.T, index=features,
                              columns=["PC%d" % (i+1) for i in
                                           range(ft_pca.components_.shape[0])])

    styled_components = components.sort_values(["PC1", "PC2"]).style.\
            apply(centered_background_gradient, cmap="PRGn", extend=0.15).\
            set_caption("Principal Components loadings").\
            set_properties(**{'max-width': '80px', 'font-size': '1pt'}).\
            set_table_styles(magnify())
    print("Rendered_components:", type(styled_components), styled_components)
    display_html(styled_components)

    # ### Feature covariance

    ft_cov = ft_pca.get_covariance()
    print("Covariance dimensions:", ft_cov.shape)
    plot_cov(ft_cov, features, cmap='seismic')
    plt.show()

    # Plot feature vectors in the PC space
    fig, (ax0, ax1) = plt.subplots(1, 2, figsize=(18,9))

    plot_features_PCspace(components, features, ax=ax0)
    plot_features_PCspace(components, features, PCs=["PC1", "PC3"], ax=ax1)
    fig.suptitle("Features in Principal Component space"); 
    plt.show()
    return ft_pca


# Functions for checking colinearity between variables

def check_decorrelate(var, correlated_var, data, logdata=None):
    _, axes = plt.subplots(2, 1 if logdata is None else 2)
    axes = axes.flat
    
    if logdata is None:
        ax0, ax2 = axes
    else:
        ax0, ax1, ax2, ax3 = axes
        
    scatter_density(correlated_var, var, alpha=0.5, s=9, data=data, ax=ax0)
    
    # Plot the uncorrelated variables
    scatter_density(data[correlated_var], data[var] / data[correlated_var], s=9, alpha=0.5, ax=ax2)
    
    ax0.set_title("Original scale")
    ax0.set_ylabel("Original %s" % var)
    ax2.set_ylabel("Uncorrelated %s" % var)
    ax2.set_xlabel(correlated_var)
    
    if logdata is not None:
        # Plot the log of uncorrelated variables
        scatter_density(correlated_var, var, s=9, alpha=0.5, data=logdata, ax=ax1)
        scatter_density(data[correlated_var], logdata[var] - logdata[correlated_var], s=9, alpha=0.5, ax=ax3)
        ax3.set_xlabel("log(%s)" % correlated_var)
        ax1.set_title("Log-scale")


# Functions for linear models

def lm_summary(lm, features, response, data):
    """Display the summary statistics of the sklearn multiple linear regression."""
    print("R^2       =", lm.score(data[features], data[response]))
    print("Intercept =", lm.intercept_)
    print("\nSlopes\n------")

    features_by_coef = sorted(zip(features, lm.coef_), key=lambda x: np.abs(x[1]), reverse=True) 

    for ft, coef in features_by_coef:
        print("%-17s: %10.6f" % (ft, coef))


def sm_ols_summary(olsfit):
    r_coefs = pd.read_csv(StringIO(olsfit.summary().tables[1].as_csv()),
                          sep='\s*,\s*', index_col=0, engine='python')

    # Sort
    r_coefs['abs_coef'] = r_coefs.coef.abs()
    r_coefs.sort_values("abs_coef", ascending=False, inplace=True)
    r_coefs.drop("abs_coef", axis=1, inplace=True)
    #r_coefs_styled = r_coefs.style.apply(centered_background_gradient, axis=0, subset="coef")
    r_coefs_styled = r_coefs.style.bar(
                        subset=pd.IndexSlice[
                            olsfit.params.index.drop(["Intercept", "const"],
                                                     errors="ignore"),
                            "coef"],
                        axis=0,
                        align="zero")
    return r_coefs_styled


if __name__ == '__main__':

    # # Load data

    aS, ts, cs = load_subtree_stats("subtrees_{stattype}stats-Simiiformes.tsv")
    check_load_subtree_stats(aS, ts, cs)

    al_params = ["ingroup_glob_len", "ingroup_mean_GC", "ingroup_mean_N",
                 "ingroup_mean_gaps", "ingroup_mean_CpG", "ingroup_std_len",
                 "ingroup_std_GC",  "ingroup_std_N",  "ingroup_std_gaps",
                 "ingroup_std_CpG"]
    tree_params = ["robust"]
    cl_params = cs.columns.tolist()[1:]
    cl_params.remove('time used')

    # ### Group columns by type

    s = pd.merge(aS, ts.drop(['subgenetree'], axis=1), how='inner')

    ingroup_cols = s.columns.str.startswith('ingroup')
    outgroup_cols = ~ingroup_cols
    # Keep the genetree column
    ingroup_cols[0] = True

    s_out = s[s.columns[outgroup_cols]]
    s_in  = s[s.columns[ingroup_cols]]

    glob_cols = s_out.columns.str.startswith('glob')
    mean_cols = s_out.columns.str.startswith('mean')
    med_cols  = s_out.columns.str.startswith('med')
    std_cols  = s_out.columns.str.startswith('std')
    w_mean_cols = s_out.columns.str.startswith('w_mean')
    w_std_cols  = s_out.columns.str.startswith('w_std')

    s_out_glob = s_out[['genetree'] + s_out.columns[glob_cols].tolist()]
    s_out_mean = s_out[['genetree'] + s_out.columns[mean_cols].tolist()]
    s_out_med =  s_out[['genetree'] + s_out.columns[med_cols].tolist()]
    s_out_std =  s_out[['genetree'] + s_out.columns[std_cols].tolist()]
    s_out_w_mean = s_out[['genetree'] + s_out.columns[w_mean_cols].tolist()]
    s_out_w_std =  s_out[['genetree'] + s_out.columns[w_std_cols].tolist()]

    s_in_glob = s_in[['genetree'] + s_in.columns[glob_cols].tolist()]
    s_in_mean = s_in[['genetree'] + s_in.columns[mean_cols].tolist()]
    s_in_med =  s_in[['genetree'] + s_in.columns[med_cols].tolist()]
    s_in_std =  s_in[['genetree'] + s_in.columns[std_cols].tolist()]
    s_in_w_mean = s_in[['genetree'] + s_in.columns[w_mean_cols].tolist()]
    s_in_w_std =  s_in[['genetree'] + s_in.columns[w_std_cols].tolist()]

    s_in.head()
    s_out.head()
    s_out_glob.head()
    s.columns


    # ## Load ages

    ages_file = "../ages/Simiiformes_m1w04_ages.subtreesCleanO2-um2-ci-grepoutSG.tsv"
    outbase = "Simiiformes_m1w04_ages.subtreesCleanO2-um2-ci-grepoutSG"

    ages = pd.read_table(ages_file, sep='\t', index_col=0)
    #if not set(('taxon', 'genetree')) & set(ages.columns):
    #    ages = splitname2taxongenetree(ages, "name")

    # ### Compute the number of duplications and speciations in the tree

    #add_robust_info(ages, ts)

    Ndup = ages.groupby('subgenetree').type.agg(lambda v: sum(v == "dup"))
    Ndup.name = 'Ndup'
    Ndup.describe()

    ages = merge_criterion_in_ages(Ndup, ages)

    Nspe = ages.groupby('subgenetree').type.agg(lambda v: sum(v == "spe"))
    Nspe.name = 'Nspe'
    Nspe.describe()

    ages = merge_criterion_in_ages(Nspe, ages)

    robust_info = pd.concat((ts, Ndup, Nspe), join='outer', axis=1, sort=False)

    robust_info.shape

    robust_info[~robust_info.robust & (robust_info.Ndup == 0) & (robust_info.Nspe == 7)]

    # First row has a _split gene_ (removed by me, so shows up as correct in terms of Ndup and Nspe).
    bad_robusts = robust_info[robust_info.robust & ((robust_info.Ndup > 0) | (robust_info.Nspe != 7))]
    print(bad_robusts.shape)
    bad_robusts.head()

    robust_info.root_location.unique()

    ages = merge_criterion_in_ages(robust_info.robust, ages)

    print(ages.columns)
    ages.head()


    # ### Fetch parent node info

    # LEFT JOIN to keep 'Simiiformes', or INNER JOIN to discard it.
    ages_p = pd.merge(ages,
                      ages[['taxon', 'type', 'age_t', 'age_dS', 'age_dN',
                            'age_dist', 'calibrated']],
                      how="left", left_on="parent", right_index=True,
                      suffixes=('', '_parent'))

    # Select only branches without duplications (from a speciation to another)
    ages_spe2spe = ages_p[(ages_p.type.isin(('spe', 'leaf'))) & (ages_p.type_parent == 'spe')]


    # ### Subset node ages data (robusts gene trees)

    ages_nodup = ages_p[ages_p.Ndup == 0].drop("Ndup", axis=1)

    ages_robust = ages_p[ages_p.robust & (ages_p.Ndup == 0) & (ages_p.Nspe == 7)].drop(["robust", "Ndup", "Nspe"], axis=1)

    ages_robust.shape

    # Robust AND with valid parent node
    ages_best = ages_spe2spe[ages_spe2spe.robust & (ages_spe2spe.Ndup == 0) & (ages_spe2spe.Nspe == 7)].drop(['Ndup', 'Nspe'], axis=1)
    print(ages_best.shape)


    # ### Merge control dates

    #add_control_dates_lengths(ages_robust, phyltreefile, timetree_ages_CI=None)

    median_taxon_ages = ages_robust[ages_robust.type.isin(("spe", "leaf"))]\
                                   .groupby("taxon").age_dS.median()
                                   #& (ages_robust.taxon != 'Simiiformes')]\
    median_taxon_ages.name = 'median_taxon_age'

    # Comparison with TimeTree dates

    phyltree = myPhylTree.PhylogeneticTree("/users/ldog/glouvel/ws2/DUPLI_data85/PhylTree.TimeTree2018.Ensembl-like.nwk")

    timetree_ages = median_taxon_ages.index.to_series().apply(phyltree.ages.get)
    timetree_ages.name = 'timetree_age'

    # Confidence intervals from TimeTree.org
    timetree_ages_CI = pd.DataFrame([[41,   46],
                                     [27.6, 31.3],
                                     [18.6, 21.8],
                                     [14.7, 16.8],
                                     [ 9.0, 14.9],
                                     [ 8.4, 9.7],
                                     [6.23, 7.07]],
                                    columns=['timetree_CI_inf', 'timetree_CI_sup'],
                                    index=['Simiiformes', 'Catarrhini',
                                        'Hominoidea', 'Hominidae',
                                        'Cercopithecinae', 'Homininae', 'HomoPan'])

    control_ages = pd.concat((median_taxon_ages, timetree_ages, timetree_ages_CI), axis=1, sort=False)

    print(control_ages.sort_values('timetree_age', ascending=False).head(10))

    # Adding the median age into `ages_best`

    #ages_median_taxon_ages =
    ages_controled = pd.merge(ages_best, control_ages,
                              left_on="taxon", right_index=True, validate="many_to_one")


    # ### Merge control branch lengths

    # Control branch length in million years.
    ages_controled['median_brlen'] = \
            ages_controled.taxon_parent.apply(control_ages.median_taxon_age.get) \
            - ages_controled.median_taxon_age

    # Resulting branch lengths

    median_brlen = ages_controled[
                        ~ages_controled.duplicated(["taxon_parent", "taxon"])
                        ][
                            ["taxon_parent", "taxon", "median_brlen",
                             "median_taxon_age"]
                        ].sort_values("taxon_parent", ascending=False)

    branch_info = ["taxon_parent", "taxon"]
    median_brlen.index = pd.MultiIndex.from_arrays(
                                            median_brlen[branch_info].values.T,
                                            names=branch_info)
    median_brlen.drop(branch_info, axis=1, inplace=True)
    median_brlen

    control_treelen = median_brlen.median_brlen.sum()
    print("Control tree length (robust) =", control_treelen, "My")

    real_control_treelen = median_brlen.loc[[("Simiiformes", "Catarrhini"),
                                         ("Simiiformes", "Callithrix jacchus"),
                                         ("Catarrhini", "Cercopithecinae"),
                                         ("Cercopithecinae", "Macaca mulatta"), 
                                         ("Cercopithecinae", "Chlorocebus sabaeus"),
                                         ("Cercopithecinae", "Papio anubis"),
                                         ("Catarrhini", "Hominoidea"),
                                         ("Hominoidea", "Nomascus leucogenys"),
                                         ("Hominoidea", "Hominidae"),
                                         ("Hominidae", "Pongo abelii"),
                                         ("Hominidae", "Homininae"),
                                         ("Homininae", "Gorilla gorilla gorilla"),
                                         ("Homininae", "HomoPan"),
                                         ("HomoPan", "Homo sapiens"),
                                         ("HomoPan", "Pan troglodytes")]].median_brlen.sum()
    real_control_treelen


    # #### Checks
    #check_control_brlen(ages_controled, phyltree)

    # Check out unexpected species branches for robust trees
    ages_best[(ages_best.taxon_parent == "Hominoidea") & (ages_best.taxon == "Homininae")\
        | (ages_best.taxon_parent == "Homininae") & (ages_best.taxon == "Homo sapiens")\
        | (ages_best.taxon_parent == "Catarrhini") & (ages_best.taxon == "Macaca mulatta")\
        | (ages_best.taxon_parent == "Catarrhini") & (ages_best.taxon == "HomoPan")]

    ages_p.loc["HomininaeENSGT00390000008575.b"]
    ages_p[ages_p.subgenetree == "SimiiformesENSGT00390000008575"]

    # Ignoring the source of the problem for now, just dropping the erroneous genetree.
    # 
    # **`TODO:`** Fix the detection of duplication VS speciation node in `generate_dNdStable`

    (ages_controled.subgenetree == "SimiiformesENSGT00390000008575").any()

    erroneous_nodes = ages_controled[
                        ages_controled.subgenetree=="SimiiformesENSGT00390000008575"].index
    ages_controled.drop(erroneous_nodes, inplace=True)

    ages_best.loc[erroneous_nodes]
    ages_best.drop(erroneous_nodes, inplace=True)

    ages_best.shape, ages_controled.shape


    # ### Quality measures

    #compute_dating_errors(ages_controled)

    ages_controled["abs_age_error"] = \
                    (ages_controled.age_dS - ages_controled.median_taxon_age).abs()
    ages_controled["signed_age_error"] = \
                    (ages_controled.age_dS - ages_controled.median_taxon_age)
    ages_controled["abs_brlen_error"] = \
                    (ages_controled.age_dS_parent - ages_controled.age_dS - \
                     ages_controled.median_brlen).abs()
    ages_controled["signed_brlen_error"] = \
                    (ages_controled.age_dS_parent - ages_controled.age_dS - \
                     ages_controled.median_brlen)

    mean_errors = ages_controled[
                     ['subgenetree', 'abs_age_error', 'signed_age_error',
                      'abs_brlen_error', 'signed_brlen_error']
                    ].groupby("subgenetree").mean()

    # #### Display

    print(ages_controled.subgenetree.unique().size, mean_errors.shape)
    mean_errors.tail()

    scatter_density("abs_age_error", "signed_age_error", mean_errors, alpha=0.3);

    _, (ax0, ax1) = plt.subplots(2)
    mean_errors.abs_age_error.hist(bins=50, ax=ax0)
    np.log10(mean_errors.abs_age_error).hist(bins=50, ax=ax1);

    scatter_density("abs_brlen_error", "signed_brlen_error", data=mean_errors, alpha=0.3);

    _, (ax0, ax1) = plt.subplots(2)
    mean_errors.abs_brlen_error.hist(bins=50, ax=ax0)
    np.log10(mean_errors.abs_brlen_error).hist(bins=50, ax=ax1);

    ax = mean_errors.plot.scatter("abs_age_error", "abs_brlen_error", alpha=0.3)
    ax.set_yscale('log')
    ax.set_xscale('log')


    # ## Correct codeml summary stats with branch length information

    # Aim: add mean tree rates, by taking theoretical branch length (My) into account.

    print(cl_params)


    dist_measures = ["branch_dist", "branch_t", "branch_dS", "branch_dN"]

    cs_rates, cs_wstds = compute_branchrate_std(ages_controled, dist_measures)
    # #### Checks

    test_g = sgg.get_group('SimiiformesENSGT00390000000002.a.a.a').drop("subgenetree", axis=1)
    #print(test_g.shape)
    test_g.sort_values(branch_info)

    print(test_g[dist_measures + ['median_brlen']].sum() / real_control_treelen)
    cs.brlen_mean['SimiiformesENSGT00390000000002.a.a.a']

    test_gm = pd.merge(test_g, median_brlen, how="outer", left_on=branch_info, right_index=True, indicator=True)
    test_gm[["_merge", "taxon_parent", "taxon", "median_brlen_x", "median_brlen_y", "median_taxon_age_x", "median_taxon_age_y"]]

    sgg.median_brlen.sum().describe()

    # Still some errors remaining. **MUSTFIX**. [Edit: seems ok now] 

    all_median_brlen = sgg.median_brlen.sum()
    real_control_treelen, all_median_brlen[0]

    epsilon = 1e-13
    ((all_median_brlen - real_control_treelen).abs() < epsilon).all()
    # OK!

    cs_wstds.head()


    # #### Checks

    test_g.sort_values(branch_info)

    test_rates = cs_rates.loc["SimiiformesENSGT00390000000002.a.a.a"]
    test_rates

    test_g.branch_t / test_g.median_brlen - test_rates.t_rate

    (test_g.branch_t / test_g.median_brlen - test_rates.t_rate)**2 * test_g.median_brlen

    ((test_g.branch_t / test_g.median_brlen - test_rates.t_rate)**2 * test_g.median_brlen).sum() / real_control_treelen

    ((test_g.branch_dS / test_g.median_brlen - test_rates.dS_rate)**2 * test_g.median_brlen).sum() / real_control_treelen

    # Ok.


    # ## Merge all statistics by subgenetree

    print("Shapes:", mean_errors.shape, s[al_params].shape, cs[cl_params].shape, cs_rates.shape, cs_wstds.shape, Ndup.shape, Nspe.shape)
    print("\nParameters:")
    print(" ".join(mean_errors.columns) + "\n")
    print(" ".join(al_params) + "\n\n" + " ".join(cl_params) + "\n")
    print(" ".join(cs_rates.columns.values) + "\n")
    print(" ".join(cs_wstds.columns.values) + "\n")
    print(Ndup.name, Nspe.name)

    cs = pd.concat((cs[cl_params], cs_rates, cs_wstds), axis=1, sort=False)
    cs.shape

    params = ["ls", "ns", "Nbranches", "NnonsynSites", "NsynSites", "kappa",
              "treelen", "dN_treelen", "dS_treelen", "brOmega_mean", "brOmega_std",
              "brOmega_med", "brOmega_skew", "dist_rate", "t_rate",
              "dS_rate", "dN_rate", "dist_rate_std", "t_rate_std", "dS_rate_std",
              "dN_rate_std", "lnL", "Niter", "seconds"]

    alls = pd.concat([mean_errors, Ndup, Nspe, s[al_params], cs[params]], axis=1, join='inner')


    # ### Checks

    print(alls.shape)
    print(alls.columns)
    alls.head()

    # Check that we got only real robust trees

    alls.ns.describe()  # Should always be 10 or 11 (9 ingroup, +1 or 2 outgroup)
    alls[~alls.ns.isin((10, 11))]
    alls[alls.ns == 11]
    "SimiiformesENSGT00390000000097" in alls.index

    ages_p[ages_p.subgenetree == "SimiiformesENSGT00390000000097"].sort_values("taxon_parent")

    # Why isn't there the `Macaca.mulattaENSGT...` entry???
    # 
    # **OK**: it's because it was a identified as a split gene but those where greped out...


    # # Overview of distributions

    # ## Alignment stats

    s_out_glob.hist(bins=20);

    axes = s[['glob_len', 'ingroup_glob_len']].hist(bins=1500, layout=(2,1), sharex=True)
    axes[0, 0].set_xlim(0, 10000)

    s.ingroup_glob_len.quantile([0.25, 0.75])

    s_in_mean.hist(bins=30);
    s_in_std.hist(bins=30);
    s_in_med.hist(bins=30);


    # ## Codeml output stats

    axes = cs[params].hist(bins=50)

    print(cs.dS_treelen.max())
    dS_treelen_heights, dS_treelen_bins, dS_treelen_patches = plt.hist(np.log10(cs.dS_treelen), bins=50)

    plt.hist(np.log10(cs.NsynSites), bins=50);
    plt.hist(np.log10(cs.brdS_mean), bins=50);


    # ## New codeml rate stats

    cs_rates.hist(bins=50)
    np.log10(cs_rates).hist(bins=50)
    cs_wstds.hist(bins=50)
    np.log10(cs_wstds).hist(bins=50);


    # # Subset data based on criterion

    # ## step by step example (ignore)

    # We can split the data based on a %GC threshold of 0.52

    s.ingroup_mean_GC[s.ingroup_mean_GC >= 0.52].to_csv("Simiiformes_topGC.tsv", sep='\t', header=False)
    s.ingroup_mean_GC[s.ingroup_mean_GC < 0.52].to_csv("Simiiformes_lowGC.tsv", sep='\t', header=False)

    #get_ipython().magic('pinfo pd.merge')
    s_in_glob.columns
    s_in_glob.head()

    # ## function definition


    # ## GC content

    #subset_on_criterion_tails(s.ingroup_glob_GC, ages, outbase=outbase, criterion_name="GC")


    # ### Step by step setup (ignore)

    ages_GC = merge_criterion_in_ages(s.ingroup_mean_GC, ages=ages, criterion_name="GC")
    q = 0.25
    low_lim, high_lim = s.ingroup_mean_GC.quantile([q, 1. - q])
    print(low_lim, high_lim)

    ages_lowGC = ages_GC[(ages_GC["GC"] <= low_lim) & (ages_GC.type == "spe")].copy()
    ages_highGC = ages_GC[(ages_GC["GC"] >= high_lim) & (ages_GC.type == "spe")].copy()

    groups_highGC = ages_highGC.groupby("taxon")
    groups_lowGC = ages_lowGC.groupby("taxon")

    groups_lowGC.age_dS.describe()
    groups_highGC.age_dS.describe()
    groups_lowGC.age_dS.median()
    groups_highGC.age_dS.median()

    ax = ages_lowGC.boxplot("age_dS", by="taxon", widths=0.4)
    ages_highGC.boxplot("age_dS", by="taxon", ax=ax, widths=0.4, positions=np.arange(len(groups_highGC.groups))+0.4)

    print(ages_lowGC.shape, ages_highGC.shape)

    ages_lowGC["tail"] = "lowGC"
    ages_highGC["tail"] = "highGC"
    #pd.concat([ages_lowGC, ages_highGC])

    ages_by_GC = pd.concat([ages_lowGC, ages_highGC])
    ages_by_GC.shape
    #ages_GC[ages_GC.t.group.median()

    sb.violinplot(x="taxon", y="age_dS", hue="tail", data=ages_by_GC, split=True)


    # ### Concise execution

    ages_GC = annot_quantiles_on_criterion(ages, s.ingroup_glob_GC, "GC", nquantiles=10)
    ages_GC.Q_GC.unique()
    ages_GC[["GC", "Q_GC"]].head()
    ax = _violin_spe_ages_vs_criterion_quantiles(ages_GC, "GC", split=True)

    ages_nodup_GC = annot_quantiles_on_criterion(ages_nodup, s.ingroup_glob_GC,
                                                 "GC", nquantiles=10)

    ax = _violin_spe_ages_vs_criterion_quantiles(ages_nodup_GC, "GC", split=True)

    ages_GC.plot.scatter("GC", "abs_error_mean");


    # ## Aligned length

    ages_len = annot_quantiles_on_criterion(ages, s.ingroup_glob_len, criterion_name="len")

    subset_on_criterion_tails(s.ingroup_glob_len, ages, outbase=outbase, criterion_name="length", save=False)

    get_ipython().run_cell_magic('bash', '', '''dSvisualizor.py tree \
    --sharescale \
    -t "Alignment length (without outgroup) <= 813" \
    -p ~/ws2/DUPLI_data85/PhylTree.TimeTree2018.Ensembl-like.nwk \
    Simiiformes_m1w04_ages.subtreesCleanO2-um2-ci-grepoutSG-lengthlowQ4.tsv \
    Simiiformes_m1w04_ages.subtreesCleanO2-um2-ci-grepoutSG-lengthlowQ4.svg''')

    get_ipython().run_cell_magic('bash', '', '''dSvisualizor.py tree \
    --sharescale \
    -t "Alignment length (without outgroup) >= 2094" \
    -p ~/ws2/DUPLI_data85/PhylTree.TimeTree2018.Ensembl-like.nwk \
    Simiiformes_m1w04_ages.subtreesCleanO2-um2-ci-grepoutSG-lengthhighQ4.tsv \
    Simiiformes_m1w04_ages.subtreesCleanO2-um2-ci-grepoutSG-lengthhighQ4.svg''')

    #_violin_spe_ages_vs_criterion_quantiles(ages_len, "len", split=False)
    violin_spe_ages_vs_criterion(ages_best, alls.ingroup_glob_len, split=True)

    ages_len.plot.scatter("len", "abs_error_mean", alpha=0.25, logy=True);


    # ## N content

    subset_on_criterion_tails(s.glob_N, ages, outbase=outbase, criterion_name="N")

    ages_N = annot_quantiles_on_criterion(ages, s.ingroup_mean_N, criterion_name="N", nquantiles=50)

    ages_N.Q_N.unique()

    _violin_spe_ages_vs_criterion_quantiles(ages_N, "N")


    # ## gap content

    subset_on_criterion_tails(s.ingroup_glob_gaps, ages, outbase=outbase, criterion_name="gaps")

    get_ipython().run_cell_magic('bash', '', '''dSvisualizor.py tree \
    --sharescale \
    -t "%Gaps (without outgroup) <= 0.0029" \
    -p ~/ws2/DUPLI_data85/PhylTree.TimeTree2018.Ensembl-like.nwk \
    Simiiformes_m1w04_ages.subtreesCleanO2-um2-ci-grepoutSG-gapslowQ4.tsv \
    Simiiformes_m1w04_ages.subtreesCleanO2-um2-ci-grepoutSG-gapslowQ4.svg''')

    get_ipython().run_cell_magic('bash', '', '''dSvisualizor.py tree \
    --sharescale \
    -t "%Gaps (without outgroup) >= 0.11" \
    -p ~/ws2/DUPLI_data85/PhylTree.TimeTree2018.Ensembl-like.nwk \
    Simiiformes_m1w04_ages.subtreesCleanO2-um2-ci-grepoutSG-gapshighQ4.tsv \
    Simiiformes_m1w04_ages.subtreesCleanO2-um2-ci-grepoutSG-gapshighQ4.svg''')

    ages_gaps = annot_quantiles_on_criterion(ages, s.ingroup_mean_gaps,
                                             criterion_name="gaps", nquantiles=10)

    _violin_spe_ages_vs_criterion_quantiles(ages_gaps, "gaps", split=False)


    # ## GC heterogeneity
    # 
    # (Intra-alignment GC standard-deviation) 

    ages_stdGC = annot_quantiles_on_criterion(ages, s.ingroup_std_GC,
                                              criterion_name="stdGC", nquantiles=4)

    _violin_spe_ages_vs_criterion_quantiles(ages_stdGC, "stdGC")


    # ## Number of duplications in the subtree

    sb.violinplot("taxon", "age_dS", data=ages[ages.type == "spe"], hue=(ages.Ndup > 0));


    # ## tree length in dS

    violin_spe_ages_vs_criterion(ages, cs.dS_treelen)


    # ## mean branch length in dS

    violin_spe_ages_vs_criterion(ages, cs.brdS_mean)
    violin_spe_ages_vs_criterion(ages_nodup, cs.brdS_mean)


    # ## mean branch length in dN

    violin_spe_ages_vs_criterion(ages, cs.brdN_mean, split=False)


    # ## branch length std

    violin_spe_ages_vs_criterion(ages, cs.brlen_std / cs.brlen_mean, criterion_name="Rbrlen_std", split=False)

    ages.subgenetree.unique().size
    ages_robust.subgenetree.unique().size

    violin_spe_ages_vs_criterion(ages_robust, cs.brdN_mean)


    # ## Mean Omega

    violin_spe_ages_vs_criterion(ages_robust, cs.brOmega_mean)
    violin_spe_ages_vs_criterion(ages_robust, cs.brOmega_med)


    # ## tree dS rate (ingroup)

    set(ages_best.subgenetree.unique()) - set(cs_rates.index)
    violin_spe_ages_vs_criterion(ages_best, cs_rates.dS_rate)

    cs_wstds.dS_rate_std.sort_values(ascending=False).head(20)

    violin_spe_ages_vs_criterion(ages_best, cs_wstds.dS_rate_std, split=True)


    # Try to plot a scatter plot for each speciation

    taxon_ages = merge_criterion_in_ages(cs.dS_rate_std, ages_best).groupby("taxon")
    #ages_best.groupby('taxon')
    #taxon_ages[["age_dS", ""]]

    taxon_ages.plot("dS_rate_std", "age_dS", kind='scatter')


    # ## Number of synonymous sites

    ages_Nsynsites = annot_quantiles_on_criterion(ages, cs.NsynSites, criterion_name="Nsynsites", nquantiles=10)

    ax = violin_spe_ages_vs_criterion(ages_Nsynsites[ages_Nsynsites.taxon != 'Simiiformes'], "Nsynsites", isin=list(range(10)), split=False)


# # Multiple linear regression

    # We'll fit our measured error to all of our parameters using a linear model

    responses = ['abs_age_error', 'signed_age_error', 'abs_brlen_error', 'signed_brlen_error']
    features = al_params + params
    # because using robust trees:
    features.remove('Nbranches')
    features.remove('ns')
    features.remove('Niter')

    # Suggested additional variables:
    # 
    # - position in genome
    # - position relative to centromeres/telomeres
    # - recombinaison
    # - sequencing error rate
    # - protein function
    # - number of alternative transcripts

    # ## Transform variables

    # Remove NaN values

    print(alls.shape)
    alls.isna().sum(axis=0).sort_values(ascending=False).head()

    alls_nona = alls[~alls.abs_brlen_error.isna()]

    print(alls_nona.shape)
    alls_nona.isna().sum(axis=0).head()

    # Check if variables should be log-transformed/square-root-transformed.
        
    all_test_transforms(alls, responses+features)

    # Transform variables

    # Choose function to apply based on the distributions above.

    totransform = {"abs_age_error":      log,
                   "signed_age_error":   make_logpostransform_inc(1.),
                   "abs_brlen_error":    log,
                   "signed_brlen_error": make_logpostransform_inc(1.),
                   "Ndup":              notransform,
                   "Nspe":              notransform,
                   "ingroup_glob_len":  log,
                   "ingroup_mean_GC":   notransform,
                   "ingroup_mean_N":    notransform,
                   "ingroup_mean_gaps": sqrt,
                   "ingroup_mean_CpG":  make_logtransform_inc(0.02), # ~ 10**(-1.71487)
                   "ingroup_std_len":   make_logtransform_inc(1.),
                   "ingroup_std_GC":    log,
                   "ingroup_std_N":     notransform,
                   "ingroup_std_gaps":  sqrt,
                   "ingroup_std_CpG":   make_logtransform_inc(0.00165), # ~ 10**(-2.78173)
                   "ls":                log,
                   "ns":                notransform,
                   "Nbranches":         notransform,
                   "treelen":           log,
                   "NnonsynSites":      log,
                   "NsynSites":         log,
                   "kappa":             log,
                   #"brlen_mean":        log,
                   #"brlen_std":         log,
                   #"brlen_med":         log,
                   #"brlen_skew":        make_logtransform_inc(2.1483),  # ~ 10**(0.332096)
                   "dN_treelen":        make_logtransform_inc(0.0254),  # ~ 10**(-1.59517)
                   "dS_treelen":        log,
                   #"brdS_mean":         log,
                   #"brdS_std":          log,
                   #"brdS_med":          make_logtransform_inc(0.005158),  # ~ 10**(-2.28751)
                   #"brdS_skew":         notransform,
                   #"brdN_mean":         log,
                   #"brdN_std":          log,
                   #"brdN_med":          make_logtransform_inc(0.000006),  
                   #"brdN_skew":         notransform,
                   "brOmega_mean":      log,
                   "brOmega_std":       sqrt, #make_logtransform(),  # ~ 10**(-0.838446)
                   "brOmega_med":       log,  # ~ 10**()
                   "brOmega_skew":      notransform,
                   "dist_rate":         log,
                   "t_rate":            log,
                   "dS_rate":           log,
                   "dN_rate":           log,
                   "dist_rate_std":     log,
                   "t_rate_std":        log,
                   "dS_rate_std":       log,
                   "dN_rate_std":       log,
                   "lnL":               logneg,
                   "Niter":             notransform,
                   "seconds":           log}

    new_feature_names = {ft: (func.__name__ + ("" if "%s" in func.__name__ else "(%s)")) % ft for ft, func in totransform.items() if func.__name__ != "<lambda>"}

    alls_transformed = alls.transform(totransform)
    assert not alls_transformed.isna().sum(axis=0).any()

    # Standardize/center variables

    forgotten = set(alls.columns) - set(totransform)
    assert not forgotten, forgotten

    zscore = lambda x: (x - x.mean()) / x.std()
    alls_normed = alls_transformed.transform({ft: zscore for ft in features})
    alls_normed[responses] = alls_transformed[responses]

    print(alls_normed.shape)
    print("NaNs:", alls_normed.columns.values[alls_normed.isna().sum(axis=0) > 0])
    alls_normed.head()

    alls_normed.hist(bins=50);

    Y = alls_normed.abs_brlen_error
    X = alls_normed[features]

    # ## PCA of parameters

    # colinearity?

    mpl.rcParams['figure.figsize'] = (16, 8) # width, height

    from sklearn.decomposition import PCA
    #import sklearn_panda

    detailed_pca(alls_normed, features)

    # ## Check variable colinearity

    # From the above covariance matrix, the following features seem to strongly covary:
    # 
    # - Ndup, ns, Nbranches, Niter
    # - ingroup_mean_CpG, ingroup_mean_GC
    # - ingroup_mean_N, ingroup_std_N
    # - treelen, brlen_mean, brlen_std, dN_treelen, dS_treelen, brdS_mean, brdS_std
    # - ingroup_glob_len, lnL

    _, (ax0, ax1) = plt.subplots(1,2)
    ax0.plot(alls.treelen, alls.brlen_mean*alls.Nbranches, '.', alpha=0.5)
    ax1.plot(X.treelen, X.brlen_mean + np.log10(X.Nbranches), '.', alpha=0.5);

    # CpG should be normalized by GC
    _, (ax0, ax1, ax2) = plt.subplots(1,3, sharex=True, sharey=False)
    ax0.plot(alls.ingroup_mean_GC, alls.ingroup_mean_CpG, '.', alpha=0.5);
    ax1.plot(alls.ingroup_mean_GC, alls.ingroup_mean_CpG/(alls.ingroup_mean_GC**2), '.', alpha=0.5);

    ax2.plot(s.ingroup_mean_GC, s.ingroup_mean_CpG/(s.ingroup_mean_G*s.ingroup_mean_C), '.', alpha=0.5);

    #a_i = alls_inde  not yet defined
    scatter_density(a.t_rate, a.dS_rate + a.dN_rate)

    _, (ax0, ax1) = plt.subplots(1, 2)
    ax0.plot(X.ns, X.Nbranches, '.', alpha=0.5)
    ax0.set_xlabel("N sequences")
    ax0.set_ylabel("N branches")
    ax1.plot(X.ns, X.Ndup, '.', alpha=0.5)
    ax1.set_xlabel("N sequences")
    ax1.set_ylabel("N dup");

    # _, (ax0, ax1) = plt.subplots(1,2)
    # ax0.plot("ingroup_glob_len", "lnL", '.', alpha=0.5, data=X)
    # ax1.plot("ls", "lnL", '.', alpha=0.5, data=X);
    check_decorrelate("lnL", "ingroup_glob_len", alls, X)

    a = alls
    plt.plot(a.NsynSites / (a.ls*3), a.NnonsynSites / (a.ls*3), '.', alpha=0.7); # Alright, Ok.

    check_decorrelate("NsynSites", "ls", alls, X)
    check_decorrelate("ingroup_std_N", "ingroup_mean_N", alls, logdata=None)
    check_decorrelate("brlen_std", "brlen_mean", alls, alls_transformed)
    check_decorrelate("brdS_std", "brdS_mean", alls, X)
    #plt.plot(X.brdS_mean, X.brdS_std / X.brdS_mean, '.', alpha=0.5);

    check_decorrelate("t_rate_std", "t_rate", alls, alls_transformed)
    check_decorrelate("dS_rate_std", "dS_rate", alls, X)
    #plt.plot(X.brdS_mean, X.brdS_std / X.brdS_mean, '.', alpha=0.5);

    check_decorrelate("dN_rate_std", "dN_rate", alls, X)
    #plt.plot(X.brdS_mean, X.brdS_std / X.brdS_mean, '.', alpha=0.5);

    check_decorrelate("dist_rate_std", "dist_rate", alls, X)
    #plt.plot(X.brdS_mean, X.brdS_std / X.brdS_mean, '.', alpha=0.5);

    # There is an outlier point with super high dS_rate and dN_rate and std for each
    alls.sort_values(["dS_rate", "dS_rate_std"], ascending=False).head(10)
    alls.sort_values(["dS_rate_std", "dS_rate"], ascending=False).head(10)

    "SimiiformesENSGT00760000119097.F" in ages_best.subgenetree.values
    ts.loc["SimiiformesENSGT00760000119097.F"]
    alls.sort_values(["dN_rate_std", "dN_rate"], ascending=False).head(10)
    alls.sort_values(["t_rate_std", "t_rate"], ascending=False).head(10)

    fig, axes = plt.subplots(2,2)
    xyvars = ["ingroup_mean_gaps", "ingroup_std_gaps"]

    scatter_density(*xyvars, data=alls, ax=axes[0,0], s=9, alpha=0.6)
    scatter_density(*xyvars, data=np.sqrt(alls[xyvars]), ax=axes[0,1], s=9, alpha=0.6)

    decorr_var = alls.ingroup_std_gaps / alls.ingroup_mean_gaps
    decorr_var[decorr_var.isna()] = 0
    scatter_density(alls.ingroup_mean_gaps, decorr_var, ax=axes[1,0], s=9, alpha=0.6)
    scatter_density(np.sqrt(alls.ingroup_mean_gaps), np.sqrt(decorr_var), ax=axes[1,1], s=9, alpha=0.6);

    plt.plot(a.treelen,
            (a.dS_treelen * a.NsynSites + a.dN_treelen * a.NnonsynSites)\
                /(a.NsynSites + a.NnonsynSites),
            '.', alpha=0.5);
    # Alright, OK.

    plt.plot(a.ls, (a.NnonsynSites + a.NsynSites)*3, '.', alpha=0.5); # Alright
    plt.plot(a.ls, a.ingroup_glob_len, '.', alpha=0.5);

    scatter_density(X.ingroup_glob_len, X.ls + X.ingroup_mean_gaps, s=9, alpha=0.5);


    # ## Create independent variables

    a_n = alls
    a_t = alls_transformed

    alls_inde = a_t[responses + ["seconds", #"ns", #"Ndup", "Nspe",
                                 "ingroup_glob_len", "ingroup_mean_gaps", "ingroup_mean_N", "kappa",
                                 #"brlen_mean", "brlen_med", "brlen_skew", "brdN_mean",
                                 "brOmega_mean", "brOmega_med", "brOmega_skew",
                                 "dist_rate", "t_rate", "dS_rate", "dN_rate"]].copy()

    # Dropped features:
    # - ns (~Ndup)
    # - ls, NsynSites (~ingroup_glob_len)
    # - treelen (~brlen_mean)

    #alls_inde["Rdup"] = a_n.Ndup / a_n.ns
    alls_inde["Ringroup_std_gaps"] = a_t.ingroup_std_gaps / a_t.ingroup_mean_gaps  # (sqrt transformed)
    alls_inde.Ringroup_std_gaps[alls_inde.Ringroup_std_gaps.isna()] = 0

    #alls_inde["RbrdS_mean"] = a_t.brdS_mean - a_t.brlen_mean
    #alls_inde["RbrdN_mean"] = a_t.brdN_mean - a_t.brlen_mean
    alls_inde["sitelnL"] = a_t.lnL / a_t.ingroup_glob_len
    alls_inde["RsynSites"] = a_t.NsynSites + np.log10(3) - a_t.ls
    alls_inde["Ringroup_std_len"] = a_t.ingroup_std_len - a_t.ingroup_glob_len
    alls_inde["Ringroup_std_N"] = a_t.ingroup_std_N - a_n.ingroup_mean_N
    #alls_inde["RbrdS_std"] = a_t.brdS_mean - a_t.brdS_std
    #alls_inde["Rbrlen_std"] = a_t.brlen_std - a_t.brlen_mean
    #alls_inde["RbrdS_std"] = a_t.brdS_std - a_t.brdS_mean
    alls_inde["R_t_rate_std"]    = a_t.t_rate_std    - a_t.t_rate
    alls_inde["R_dS_rate_std"]   = a_t.dS_rate_std   - a_t.dS_rate
    alls_inde["R_dN_rate_std"]   = a_t.dN_rate_std   - a_t.dN_rate
    alls_inde["R_dist_rate_std"] = a_t.dist_rate_std - a_t.dist_rate

    alls_inde["RbrOmega_std"] = a_t.brOmega_std - a_t.brOmega_mean
    #alls_inde["RbrdN_std"] = a_t.brdN_std - a_t.brdN_mean
    alls_inde["CpG_odds"] = a_n.ingroup_mean_CpG / (a_n.ingroup_mean_GC**2) # Assuming %C == %G

    inde_features = ["seconds", #"ns", "Ndup", "Nspe",
                     "ingroup_glob_len", "ingroup_mean_gaps", "ingroup_mean_N", "kappa",
                     #"brlen_mean", "brdS_mean", "brdN_mean", "brlen_med", "brlen_skew",
                     #"Rdup", "RbrdS_mean", "RbrdN_mean",
                     "t_rate", "dN_rate", "dS_rate", "dist_rate",
                     "brOmega_mean", "brOmega_med", "brOmega_skew",
                     "sitelnL", "RsynSites", "Ringroup_std_len", "Ringroup_std_N", "Ringroup_std_gaps",
                     #"Rbrlen_std", "RbrdS_std",
                     "R_t_rate_std", "R_dS_rate_std", "R_dN_rate_std", "R_dist_rate_std",
                     "RbrOmega_std",
                     "CpG_odds"]

    print(set(alls_inde.columns) - set(inde_features))
    print(set(inde_features) - set(alls_inde.columns))

    alls_inde_normed = alls_inde.drop(responses, axis=1).transform(zscore)
    alls_inde_normed[responses] = alls_inde[responses]

    print(set(alls_inde_normed.columns) - set(inde_features))
    print(set(inde_features) - set(alls_inde_normed.columns))

    alls_inde_normed.head()
    alls_inde.ingroup_mean_N.describe()
    alls.ingroup_mean_N.isna().head()
    s.ingroup_mean_N.describe()

    alls_inde_normed.hist(bins=50);


    # ## PCA on the new supposedly independent variables

    ft_pca_inde = PCA(n_components=12)
    ft_pca_inde_components = ft_pca_inde.fit_transform(alls_inde_normed[inde_features])
    plot_cov(ft_pca_inde.get_covariance(), inde_features)
    ft_cov_inde = ft_pca_inde.get_covariance()

    print(inde_features.index("brlen_mean"))
    print(inde_features.index("RbrdS_mean"))
    ft_cov_inde[5,6]
    ft_cov_inde[9,6]

    print(ft_cov_inde.shape, len(inde_features))

    pd.DataFrame(ft_cov_inde, index=inde_features, columns=inde_features)\
            .style.apply(centered_background_gradient, extend=0.1, axis=None)\
            .set_properties(**{'max-width': '80px', 'font-size': '1pt'})\
            .set_table_styles(magnify())


    # ## Linear Regression with Scikit-learn

    from sklearn.linear_model import LinearRegression
    from sklearn.preprocessing import StandardScaler

    lm = LinearRegression()
    lm.fit(X, y)
    lm_summary(lm, features, "abs_error", alls_normed)

    # Only the non-covariating variables

    lm2 = LinearRegression()
    lm2.fit(alls_inde_normed[inde_features], alls_inde_normed.abs_error)
    lm_summary(lm2, inde_features, "abs_error", alls_inde_normed)


    # ## Regression with StatsModels

    # ### With all non-colinear transformed features

    import statsmodels.api as sm
    #import statsmodels.formula.api as smf


    # Add intercept
    olsfit0 = sm.OLS(alls_inde_normed.abs_brlen_error, sm.add_constant(alls_inde_normed[inde_features])).fit()
    olsfit0.params
    olsfit0.summary()
    sm_ols_summary(olsfit0)

    # #### robust trees

    alls_inde_normed.head()

    # Add intercept
    data = alls_inde_normed[(alls_inde_normed.Ndup == 0) & (alls_inde_normed.Nspe == 7)]
    olsfitr = sm.OLS(data.abs_error, sm.add_constant(data[inde_features])).fit()
    olsfitr.params
    olsfitr.summary()
    sm_ols_summary(olsfitr)


    # ### Same with the formula syntax

    formula = 'abs_error ~ ' + ' + '.join(inde_features)
    print(formula)
    ols = smf.ols(formula, data=alls_inde_normed)
    results = ols.fit()

    results.params
    r_summary = results.summary()
    r_summary

    sm_ols_summary(results)


    # ### Add square effects

    # Test for an _optimum_ branch length:
    # 
    # - branches _too long_: saturation leads to bad dating?
    # - branches _too short_ (not enough substitutions): lack of information/too high stochasticity leads to bad dating?

    fig, ((ax0, ax1), (ax2, ax3)) = plt.subplots(2,2)
    ax0.hist(alls_inde_normed.brdS_mean_sq, bins=50);
    ax1.hist(np.log10(alls_inde_normed.brdS_mean_sq), bins=50);
    ax2.hist(alls_inde_normed.brlen_mean_sq, bins=50);
    ax3.hist(np.log10(alls_inde_normed.brlen_mean_sq), bins=50);

    alls_inde_normed['RbrdS_mean_sq'] = alls_inde_normed.RbrdS_mean ** 2
    # TODO: square the untransformed data
    alls_inde_normed['brlen_mean_sq'] = alls_inde_normed.brlen_mean ** 2

    # Check with the log of the squared variable (which was already logged and
    # centered normalized)

    olsfit_sq = sm.OLS(alls_inde_normed.abs_error,
                       sm.add_constant(alls_inde_normed[inde_features + \
                                            ["RbrdS_mean_sq", "brlen_mean_sq"]]
                                       )
                       ).fit()

    sm_ols_summary(olsfit_sq)

    # There does not seem to be a squared relation for the branch length.

    # ### With only the 2 best parameters

    ols_2params = smf.ols("abs_error ~ brlen_mean + ingroup_glob_len",
                          data=alls_inde_normed).fit()
    ols_2params.summary()

    # $R^2 = 0.217$

    fig, (ax0, ax1) = plt.subplots(2, 1, figsize=(15,18))

    ax0.plot("brlen_mean", "abs_error", '.', data=alls_inde_normed, alpha=0.5)
    ax0.set_title("Mean branch length (substitutions per site)")
    x0 = np.array([alls_inde_normed.brlen_mean.min(), alls_inde_normed.brlen_mean.max()])
    y0 = 0.6454 + 0.0714 * x0
    ax0.plot(x0, y0)

    ax1.plot("ingroup_glob_len", "abs_error", '.', data=alls_inde_normed, alpha=0.5)
    ax1.set_title("Length of the alignment")
    x1 = np.array([alls_inde_normed.ingroup_glob_len.min(), alls_inde_normed.ingroup_glob_len.max()])
    y1 = 0.6454 - 0.1021 * x1
    ax1.plot(x1, y1);

    # Observation : some outliers might drive the relation with `brlen_mean`.


    # ### Excluding gene trees with duplications

    print(alls_inde_normed.shape, alls.shape)
    print(alls_inde_normed[responses].head())
    print(alls[responses].head())

    alls_inde_normed_nodup = alls_inde_normed[alls.Ndup == 0].drop(["Ndup", "Rdup"], axis=1)
    inde_features_nodup = ["seconds", "ns", "ingroup_glob_len", "kappa",
                           "brlen_mean", "brdS_mean", "sitelnL", "RsynSites",
                           "Ringroup_std_len", "Ringroup_std_N", "Rbrlen_std",
                           "CpG_odds"]

    alls_inde_normed_nodup.shape
    alls_inde_normed_nodup.head()

    olsfit_nodup_sq = sm.OLS(alls_inde_normed_nodup.abs_error,
                       sm.add_constant(
                           alls_inde_normed_nodup[inde_features_nodup + \
                                            ["brdS_mean_sq", "brlen_mean_sq"]]
                           )
                       ).fit()

    sm_ols_summary(olsfit_nodup_sq)
    olsfit_nodup_sq.summary()

    fig, (ax0, ax1) = plt.subplots(2, 1, figsize=(15,18))

    ax0.plot("brlen_mean", "abs_error", '.', data=alls_inde[alls.Ndup == 0], alpha=0.5)
    ax0.set_title("Log10(Mean branch length) (substitutions per site)")
    #x0 = np.array([alls_inde_normed.brlen_mean.min(), alls_inde_normed.brlen_mean.max()])
    #y0 = 0.6454 + 0.0714 * x0
    #ax0.plot(x0, y0)

    ax1.plot("ingroup_glob_len", "abs_error", '.', data=alls_inde[alls.Ndup == 0], alpha=0.5)
    ax1.set_title("Log10(Length of the alignment)");
    #x1 = np.array([alls_inde_normed.ingroup_glob_len.min(), alls_inde_normed.ingroup_glob_len.max()])
    #y1 = 0.6454 - 0.1021 * x1
    #ax1.plot(x1, y1);


    # ## Regression in R

    #get_ipython().magic('reload_ext rpy2.ipython')
    #get_ipython().run_cell_magic('R', '-i alls_inde_normed', 'library(nlme)')

    # # Extract the worst trees

    alls_normed.sort_values('abs_error', ascending=False).head(30)

    # 5 best trees

    alls_inde_normed[responses + inde_features].sort_values('abs_brlen_error').head(5)    .style.highlight_min(subset="abs_brlen_error")

    alls[responses + features].sort_values('brOmega_mean', ascending=False).head(5)    .style.highlight_max(subset="brOmega_mean")

    ax = sb.violinplot(x="taxon", y="age_dS", data=ages_best[ages_best.type == "spe"])
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45);


# # Chronos dating (PL)

    ages_PL1 = pd.read_table("/users/ldog/glouvel/ws2/DUPLI_data85/alignments_analysis/ages/Simiiformes_m1w04_ages.subtreesCleanO2-um2-withSG-PL1.tsv",
                             sep='\t', index_col=0)

    # ### Compute the number of duplications and speciations in the tree

    Ndup_PL1 = ages_PL1.groupby('subgenetree').type.agg(lambda v: sum(v == "dup"))
    Ndup_PL1.name = 'Ndup'
    Ndup_PL1.describe()

    ages_PL1 = merge_criterion_in_ages(Ndup_PL1, ages_PL1)

    Nspe_PL1 = ages_PL1.groupby('subgenetree').type.agg(lambda v: sum(v == "spe"))
    Nspe_PL1.name = 'Nspe'
    Nspe_PL1.describe()

    ages_PL1 = merge_criterion_in_ages(Nspe_PL1, ages_PL1)

    robust_info_PL1 = pd.concat((ts, Ndup_PL1, Nspe_PL1), join='outer', axis=1, sort=False)
    robust_info_PL1.shape

    robust_info_PL1.head()

    robust_info_PL1[ (~robust_info_PL1.robust.astype(bool)) & (robust_info_PL1.Ndup == 0) & (robust_info_PL1.Nspe == 7)]

    bad_robusts_PL1 = robust_info_PL1[robust_info_PL1.robust.astype(bool) & ((robust_info_PL1.Ndup > 0) | (robust_info_PL1.Nspe != 7))]
    print(bad_robusts_PL1.shape)
    bad_robusts_PL1.head()

    robust_info_PL1.root_location.unique()

    ages_PL1 = merge_criterion_in_ages(robust_info_PL1.robust, ages_PL1)

    print(ages_PL1.columns)
    ages_PL1.head()


    # ### Fetch parent node info

    # LEFT JOIN to keep 'Simiiformes', or INNER JOIN to discard it.
    ages_PL1_p = pd.merge(ages_PL1, ages_PL1[['taxon', 'type', 'age', 'calibrated']],
                      how="left", left_on="parent", right_index=True, suffixes=('', '_parent'))

    # Select only branches without duplications (from a speciation to another)
    ages_PL1_spe2spe = ages_PL1_p[(ages_PL1_p.type.isin(('spe', 'leaf'))) & (ages_PL1_p.type_parent == 'spe')]

    # ### Subset node ages data (robusts gene trees)

    ages_nodup = ages_p[ages_p.Ndup == 0].drop("Ndup", axis=1)
    ages_robust = ages_p[ages_p.robust & (ages_p.Ndup == 0) & (ages_p.Nspe == 7)].drop(["robust", "Ndup", "Nspe"], axis=1)

    ages_robust.shape
    # Robust AND with valid parent node

    ages_PL1_best = ages_PL1_spe2spe[ages_PL1_spe2spe.robust & (ages_PL1_spe2spe.Ndup == 0) & (ages_PL1_spe2spe.Nspe == 7)].drop(['Ndup', 'Nspe'], axis=1)
    print(ages_PL1_best.shape)



