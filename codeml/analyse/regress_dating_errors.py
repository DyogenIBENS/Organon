#!/usr/bin/python3
# coding: utf-8

# Source:
# jupyter nbconvert --to python \
#   ~/ws2/DUPLI_data85/alignments_analysis/subtrees_stats/subtrees_stats.ipynb

import sys
import os.path as op

import warnings
from copy import copy, deepcopy
import numpy as np
import pandas as pd
import matplotlib as mpl
#mpl.use('TkAgg', warn=False)
#try:
#    get_ipython().magic('matplotlib inline')
#except UnknownBackend:
import matplotlib.pyplot as plt
import seaborn as sb

from seqtools.compo_freq import weighted_std
from dendro.bates import dfw_pairs_generalized, dfw_pairs
from datasci.graphs import scatter_density, \
                           plot_cov, \
                           heatmap_cov, \
                           plot_loadings, \
                           plot_features_radar, \
                           plottree, \
                           stackedbar, \
                           dodged_violin
from datasci.compare import pairwise_intersections, align_sorted
from datasci.stats import r_squared, adj_r_squared, multicol_test, multi_vartest,\
                          rescale_groups, iqr, mad, trimstd, trimmean, \
                          mean_absdevmed, f_test
from datasci.routines import *
from datasci.dataframe_recipees import *

from dendro.any import myPhylTree as phyltree_methods, ete3 as ete3_methods
import ete3

from LibsDyogen import myPhylTree

import scipy.stats as stats  # stats.levene
from sklearn.decomposition import PCA, FactorAnalysis
from sklearn.linear_model import LinearRegression, Lasso  # Lasso doesn't work.
from sklearn.preprocessing import StandardScaler
import statsmodels.api as sm
#import statsmodels.formula.api as smf
import statsmodels.stats.api as sms
import statsmodels.stats.multitest as smm
import statsmodels.graphics as smg  # smg.gofplots.qqplot

from IPython.display import display_html
import logging


logfmt = "%(levelname)-7s:l.%(lineno)3s:%(funcName)-20s:%(message)s"
logf = logging.Formatter(logfmt)

try:
    from UItools import colorlog
    #clogfmt = "$LVL%(levelname)-7s:$RESET${white}l.%(lineno)3s:%(funcName)-20s:$RESET%(message)s"
    clogfmt = "$LVL%(levelname)-7s:$RESET${white}l.%(lineno)2s:$RESET%(message)s"
    colorlogf = colorlog.ColoredFormatter(clogfmt)
except ImportError:
    colorlogf = logf

# Notebook setup
logger = logging.getLogger(__name__)
if not logger.hasHandlers():
    sh = logging.StreamHandler(sys.stdout)
    sh.setFormatter(colorlogf)
    logger.setLevel(logging.INFO)
    logger.handlers = []
    logger.addHandler(sh)
    #logging.basicConfig(handlers=[sh])
    logging.basicConfig()
# to reset:
#logging.shutdown()

mpl.style.use("softer")
pd.set_option("display.max_columns", 50)
pd.set_option("display.width", 115)
pd.set_option("display.max_colwidth", 50)
pd.set_option("display.show_dimensions", True)  # even on non truncated dataframes

# wide_screen_style
mpl.rcParams['figure.figsize'] = (14, 10) # width, height


# Useful variables: the name of statistics of interest:
## Piping of the metrics:
## measures: '{m}'

## chronos ages: 'age.{m}'     with m in ('PL1', 'PL100', 'R1', ..., 'MPL')

## age_measures -> 'age_{m}    m in ('dS', 't', 'dN', 'dist') (represents substitutions/site)   

## Original branch lengths in substitutions
## dist_measures/branch_measures: 'branch_{m}'
## same but with the parent branch: += '_parent'

## rate_measures -> 'rate_{m}' computed by dividing the branch_m by the real length (My)
## rate_std_measures -> 'rate_{m}_std' std dev of the rate in each tree.

## consecutive_null_{m}
## freq_null_{m}
## null_{m}_after
## null_{m}_before

## median_measures

MEASURES = ['dS', 'dN', 't', 'dist']
DIST_MEASURES = ['branch_%s' % m for m in MEASURES]
RATE_MEASURES = ['%s_rate' % m for m in MEASURES]
RATE_STD_MEASURES = [r + '_std' for r in RATE_MEASURES]


# Convert "time used" into seconds.  # ~~> numbertools? timetools? converters?
def time2seconds(time_str):
    factors = [1, 60, 3600, 3600*24]
    s = 0
    for factor, n_units in zip(factors, reversed(time_str.split(':'))):
        s += factor * int(n_units)
    return s


def load_stats_al(alfile):
    return pd.read_csv(alfile, index_col=0, sep='\t')

load_any_stats = load_stats_al


def load_stats_tree(treefile):
    ts = pd.read_csv(treefile, index_col=0, sep='\t',
                     dtype={'leaves_robust':      bool,
                            'single_child_nodes': bool,
                            'nodes_robust':       bool,
                            'only_treebest_spe':  int,
                            'aberrant_dists':     int})
                            #'root2tip_var': float
    ts['really_robust'] = ts.leaves_robust & ts.nodes_robust & ~ts.only_treebest_spe
    try:
        ts['unlike_clock'] = ts.root2tip_sd / ts.root2tip_mean
    except AttributeError:
        logger.warning('"root2tip_mean" or "root2tip_sd" not in treestats columns.')
    return ts

def load_stats_codeml(codemlfile):
    cs = pd.read_csv(codemlfile, index_col=0, sep='\t')
    cs['seconds'] = cs['time used'].apply(time2seconds)
    return cs


def load_stats_chronosruns(runsfile):
    rs = pd.read_csv(runsfile, sep='\t', index_col=0)
    rs.set_axis([c.replace('.', '_') for c in rs.columns], axis=1, inplace=True)
    for c in rs.columns:
        if c.endswith('_message'):
            rs[c] = rs[c].astype('category')
    return rs

def load_stats_chronoslogs(logsfile):
    ls = pd.read_csv(logsfile, sep='\t', index_col=0,
                     true_values='T', false_values='F',
                     dtype={'action': 'category'})
    return ls

stat_loaders = {'al':     load_stats_al,
                'tree':   load_stats_tree,
                'codeml': load_stats_codeml,
                'cleaning': load_any_stats,
                'chronos-runs': load_stats_chronosruns,
                'chronos-logs': load_stats_chronoslogs}

stat_loaders['treeI'] = stat_loaders['tree']
stat_loaders['codemlfsahmmc'] = stat_loaders['codemlfsa'] = stat_loaders['codemlI'] = stat_loaders['codeml']
stat_loaders['alfsahmmc'] = stat_loaders['alfsa'] = stat_loaders['alI'] = stat_loaders['al']
stat_loaders['cleaningfsa'] = stat_loaders['cleaning']


def load_subtree_stats(template, stattypes=('al', 'tree', 'codeml')):
    """
    param: template to the csv/tsv files containing al/tree/codeml stats.

    Example template: 'subtreesRawTreeBestO2_{stattype}stats_Simiiformes.tsv'
    """
    
    # ## Load tree/alignment statistics
    output_tables = []
    for stattype in stattypes:
        filename = template.format(stattype=stattype)
        logger.info('Load %s', filename)
        
        try:
            loader = stat_loaders[stattype]
        except KeyError:
            loader = load_any_stats
            logger.warning('No specific parameters for loading %r',
                    stattype)
        output_tables.append(loader(filename))

    return tuple(output_tables)


# ~~> datasci.graphs
def intersection_plot(**named_sets):
    labels, *values = pairwise_intersections(**named_sets)

    intersections = np.array(values)
    stackedbar(np.arange(len(labels)), intersections, zero=1, orientation='horizontal');
    ax = plt.gca()
    ax.set_yticks(np.arange(len(labels)))
    ax.set_yticklabels(['%s ^ %s' % lab for lab in labels], size='large')
    return ax

#def check_load_subtree_stats(aS, ts, cs, clS=None, figsize=(10,3)):
def check_subtree_stats(subtrees_stats, figsize=(10,3)):
    
    #logger.info("shapes: aS %s, ts %s, cs %s" % (aS.shape, ts.shape, cs.shape))
    #logger.info("aS has dup: %s", aS.index.has_duplicates)
    #logger.info("ts has dup: %s", ts.index.has_duplicates)
    #logger.info("cs has dup: %s", cs.index.has_duplicates)
    #common_subtrees = set(aS.index) & set(ts.index) & set(cs.index)
    #logger.info("%d common subtrees" % len(common_subtrees))
    #only_al = aS.index.difference(ts.index.union(cs.index))
    #only_tr = ts.index.difference(aS.index.union(cs.index))
    #only_co = cs.index.difference(aS.index.union(ts.index))
    #l_al = len(only_al)
    #l_tr = len(only_tr)
    #l_co = len(only_co)
    #logger.warning("%d only in al stats: %s" % (l_al, list(only_al)[:min(5, l_al)]))
    #logger.warning("%d only in tree stats: %s" % (l_tr, list(only_tr)[:min(5, l_tr)]))
    #logger.warning("%d only in codeml stats: %s" % (l_co, list(only_co)[:min(5, l_co)]))
    #
    #intersect_kwargs = dict(aS=aS.index, ts=ts.index, cs=cs.index)
    #if clS is not None:
    #    logger.info("shape clS %s; clS has dup: %s; %d common subtrees.",
    #                clS.shape, clS.index.has_duplicates,
    #                len(common_subtrees & set(clS.index)))
    #    intersect_kwargs['clS'] = clS.index

    for stype, ss in subtrees_stats:
        logger.info("%s shape: %s. Has dup: %s", stype, ss.shape, ss.index.has_duplicates)
        only_stype = ss.index.difference(set.union(*(set(oss.index)
                                            for ostype,oss in subtrees_stats
                                            if ostype != stype))
                                        )
        l_only = len(only_stype)
        if l_only:
            logger.warning("%d only in %s stats: %s", l_only, stype,
                           list(only_stype)[:min(5, l_only)])
    
    common_subtrees = set.intersection(*(set(ss.index) for _,ss in subtrees_stats))
    logger.info("%d common subtrees", len(common_subtrees))
    intersect_kwargs = {stype: ss.index for stype,ss in subtrees_stats}

    ax = intersection_plot(**intersect_kwargs)
    ax.get_figure().set_size_inches(figsize)


# ## Function to merge additional subgenetree information into `ages`
def merge_criterion_in_ages(criterion_serie, ages=None, ages_file=None,
                            criterion_name=None):
    """Merge a column into the *node ages* table: the common field is the *subgenetree*."""
    
    assert (ages is not None or ages_file) and not (ages_file and ages is not None), "At least `ages` (dataframe) or `ages_file` (filename) must be given."
    if ages is None:
        ages = pd.read_csv(ages_file, sep='\t')
    
    logger.info("Input shape: %s", ages.shape)
    criterion = criterion_serie.name if not criterion_name else criterion_name
    ages_subgenetrees = ages.subgenetree.unique()
    
    if len(criterion_serie) < ages_subgenetrees.size \
            or set(criterion_serie.index) < set(ages_subgenetrees):
        logger.error("Not all genetrees have a criterion value: %d < %d" % (len(criterion_serie), ages_subgenetrees.size))
    
    #criterion_df = pd.DataFrame({criterion: criterion_serie})
    criterion_df = criterion_serie.to_frame(criterion)
    assert (criterion_df.index == criterion_serie.index).all()
    
    ages_c = pd.merge(ages, criterion_df, how='left', left_on='subgenetree', right_index=True)
    # If all subgenetrees have a criterion, equivalent to:
    #ages[criterion] = criterion[ages.subgenetree]
    return ages_c



def load_prepare_ages(ages_file, ts, measures=['dist', 'dS', 'dN', 't']):
    """Load ages dataframe, join with parent information, and compute the
    'robust' info"""
    
    
    beast_renames = {'height': 'age_beastS',
                     'height_median': 'age_beastSmedian',
                     'length': 'branch_beastS',
                     'length_median': 'branch_beastSmedian'}
    
    ages = pd.read_csv(ages_file, sep='\t', index_col=0)\
            .rename(columns=lambda n: n.replace('.', '_'))\
            .rename(columns=beast_renames)\
            .rename_axis('name')  # If coming from date_dup.R
    if measures is None:
        measures = [c[4:] for c in ages.columns if c.startswith('age_')] 

    # Type checks:
    dtypes = {} if measures is None else {'%s_%s' % (s,m): float
                                            for s in ('age', 'branch')
                                            for m in measures}
    dtypes.update(calibrated=bool, parent=str, taxon=str, is_outgroup=bool,
                  type=str, root=str, subgenetree=str)
    try:
        ages = ages.astype(dtypes, copy=False, errors='raise')
    except KeyError as err:
        err.args += ('%r in `dtypes` not in columns.' % set(dtypes).difference(ages.columns),)
        raise

    # Convert beast string ranges to floats, in 2 columns:
    rangevars = ['height_95%_HPD', 'height_range', 'length_95%_HPD', 'length_range', 'rate_95%_HPD', 'rate_range']
    if len(ages.columns.intersection(rangevars)):
        logger.debug('ages[rangevars].dtypes = %s', ages[rangevars].dtypes)
        for var in rangevars:
            parsedrange = ages[var].replace("None", "NaN").astype(str)\
                          .str.strip('{}').str.split(r'\.\.', n=2, expand=True)\
                          .astype(float)
            logger.debug('%s: parsedrange.shape = %s; parsedrange.columns = %s',
                         var, parsedrange.shape, parsedrange.columns.tolist())
            ages[var+'_low'] = parsedrange[0]
            ages[var+'_up'] = parsedrange[1]
        ages.drop(columns=rangevars, inplace=True)

    logger.info("Shape ages: %s; has dup: %s" % (ages.shape, ages.index.has_duplicates))
    n_nodes = ages.shape[0]
    ages_int = ages[ages.type != 'leaf']
    logger.info("Shape ages internal nodes: %s; has dup: %s",
                ages_int.shape,
                ages_int.index.has_duplicates)
    n_nodes_int = (ages.type != 'leaf').sum()

    # Fetch parent node info  ## FIXME: handle missing column 'is_outgroup'
    ages_p = pd.merge(ages.reset_index(),
                      ages.loc[ages.type != 'leaf',
                            ['taxon', 'type', 'calibrated', 'is_outgroup',
                             'subgenetree']
                            + ['%s_%s' % (s,m) for s in ('age', 'branch')
                               for m in measures]],
                      how="left", left_on=["parent", "subgenetree"],
                      right_on=["name", "subgenetree"],
                      suffixes=('', '_parent'), indicator=True,
                      validate='many_to_one')\
                .set_index('name')
    logger.info("Shape ages with parent info: %s" % (ages_p.shape,))
    n_nodes_p = ages_p.shape[0]

    if n_nodes_p < n_nodes:
        logger.warning('%d nodes were lost when fetching parent information.',
                       n_nodes - n_nodes_p)

    orphans = ages_p._merge=='left_only'
    ages_orphans = ages_p[orphans]

    logger.info("\nOrphans: %d\n" % (orphans.sum()))
    #display_html(ages_orphans.head())

    #expected_orphans = ages.is_outgroup_parent.isna()
    recognized_outgroups = (ages_orphans.is_outgroup | ages_orphans.is_outgroup_parent)
                            #| ages_orphans.is_outgroup_parent == 1)
    orphan_taxa = ages_orphans[~recognized_outgroups].taxon.unique()

    #print("All orphans are expected (leaf, or child of the root).")
    logger.info('Computed ages of orphan taxa:\n%s\n%s',
                ages_orphans.loc[~recognized_outgroups,
                                 ['age_%s' %m for m in measures]].head(),
                ages_orphans.loc[~recognized_outgroups,
                                 ['age_%s' %m for m in measures]].isna().sum()\
                            .rename('Sum of NaNs').to_frame().T)
    logger.warning('CHECK those orphan taxa (should all be outgroups sequences (could be duplicates from the ingroup taxa): %s.',
                   ', '.join(orphan_taxa))
    
    e_nochild = set(ages.index) - set(ages.parent)  # expected
    parent_nodata = set(ages[ages.type!='leaf'].parent) - set(ages_p.index)
    n_nochild = len(e_nochild)
    logger.info("\nExpected nodes without children (leaves): %d" % (n_nochild,))
    logger.info("Observed nodes not found as parents: %d" % \
            len(set(ages.index) - set(ages_p[ages_p._merge=='both'].index)))
    logger.info("Parent nodes without data: %d" % (len(parent_nodata),))
    #assert len(nochild) == n_nochild, \
    #    "Found %d unexpected nodes without child:\n%s" % \
    #        (len(nochild) - n_nochild,
    #         ages.loc[nochild - e_nochild])

    # Drop outgroups! (they create duplicated index values and might confuse stuff
    # HOWEVER DON'T drop the ingroup root.
    # (doesn't have an is_outgroup_parent info)

    child_of_root = ages_p.root == ages_p.parent
    ingroup_root = ages_p.index.to_series() == ages_p.subgenetree
    #recognized_outgroups = (ages_p.is_outgroup == 1 | ages_orphans.is_outgroup_parent == 1)
    #to_remove = (~child_of_root & orphans
    #             | (ages_p.is_outgroup == 1)
    #             | ages_p.is_outgroup_parent == 1)
    #to_remove = (orphans
    #             | (ages_p.is_outgroup == 1)
    #             | ages_p.is_outgroup_parent == 1)
    to_remove = (~ingroup_root & orphans
                 | (ages_p.is_outgroup == 1)
                 | ages_p.is_outgroup_parent == 1)

    ages_p = ages_p.loc[~to_remove].copy(deep=False)
    if ages_p.index.has_duplicates:
        debug_columns = ['parent', 'subgenetree', 'taxon', 'taxon_parent',
                         'median_age_'+measures[0], 'branch_dS', 'age_dS']
        logger.error("Failed to remove index duplicates in 'ages_p':\n%s\n...",
                     ages_p.loc[ages_p.index.duplicated(), debug_columns].head(10))

    ages_spe2spe = ages_p[(ages_p.type.isin(('spe', 'leaf')))
                          & (ages_p.type_parent == 'spe')]
    
    logger.info("\nShape ages speciation to speciation branches (no dup): %s",
                ages_spe2spe.shape)

    #ages_robust = ages_treestats[ages_treestats.really_robust & \
    #                             (ages_treestats.aberrant_dists == 0)]\
    #                        .drop(['really_robust', 'aberrant_dists'], axis=1)
    return add_robust_info(ages_p, ts, measures)


def add_robust_info(ages_p, ts, measures=['dist', 'dS', 'dN', 't']):
    """
    Compute Number of duplications/speciation per tree,
    and additional tree specific statistics.
    """
    branch_measures = ['branch_'+m for m in measures]  # global dist_measures
    #NOTE: already computed by `subtrees_stats` in current version.
    for m in measures:
        logger.debug('measure: %s' %m)
        ages_p['consecutive_null_%s' %m] = \
                ~ages_p[['branch_%s_parent' %m, 'branch_%s' %m]].any(axis=1,
                                                                  skipna=False)
        #NOTE: there's a pitfall with this method as the `branch_dist` of the ingroup node is not documented by `generate_dNdS.py`!
        #-> with skipna=False, the NaN are True.

    sgg = subgenetree_groups = ages_p.groupby('subgenetree', sort=False)

    logger.info('Aggregating `ns` ("new stats" specific to this dataset)...')
    # This is a slow operation. Needs optimization.
    
    # tree statistics making use of each node.
    ns = pd.concat((sgg.type.agg(lambda v: sum(v == "dup")),
                    sgg.type.agg(lambda v: sum(v == "spe"))),
                   axis=1, keys=['Ndup', 'Nspe'])

    #print(ns.Ndup.describe())
    #print(ns.Nspe.describe())

    def freq_of_null(v):
        return (v==0).mean()

    # Already computed by codeml.subtrees_stats.get_codeml_stats
    agg_funcs = {'consecutive_null_%s' %m: 'mean' for m in measures}

    #
    agg_funcs.update({'branch_%s' %m: freq_of_null for m in measures})
    ns = ns.join(sgg.agg(agg_funcs))\
           .rename(columns={'branch_%s' %m: 'freq_null_%s' %m for m in measures})

    # Now, count null branches **around** dated nodes (before/after)

    sgg_before = ages_p[ages_p.calibrated==0].groupby('subgenetree', sort=False)
    sgg_after = ages_p[ages_p.calibrated_parent==0].groupby('subgenetree', sort=False)

    ns = ns.join(sgg_before[branch_measures].agg(freq_of_null))\
           .rename(columns={'branch_%s' %m: 'null_%s_before' %m
                            for m in measures})\
           .join(sgg_after[branch_measures].agg(freq_of_null))\
           .rename(columns={'branch_%s' %m: 'null_%s_after' %m
                            for m in measures})\
           .join(ns.index.to_series()\
                 .str.extract(r'^[A-Z][a-zA-Z_]+ENSGT00(\d\d)0.*').astype(int)\
                 .set_axis(['ensembl_version'], axis=1, inplace=False))
           #.join(sgg_after[['branch_dN', 'branch_dS']]\
           #      .apply(lambda r: (r==0).all(axis=1))

    # merge tree stats to select robust trees
    ages_treestats = pd.merge(ages_p.drop('_merge', axis=1),
                              ns[['Ndup', 'Nspe']].join(
                                  ts[['really_robust', 'aberrant_dists',
                                      'rebuilt_topo']]),  # Is this useful?
                              how='left',
                              left_on='subgenetree',
                              right_index=True,
                              indicator=True, validate='many_to_one')
    logger.info("Ndup %s; Nspe %s; ages_treestats %s",
                (~ns.Ndup.isna()).sum(), (~ns.Nspe.isna()).sum(),
                ages_treestats.shape)
    logger.info('merge types:\n%s',
                ages_treestats._merge.value_counts())

    return ages_treestats, ns
    #return ages_p, ns


# for averaging by taking into account branch length: with Omega.
# NOTE: columns will be reordered following `var`.
def old_group_average(g, var, weight_var="median_brlen"):
    g = g[~g[weight_var].isna()]
    #g = g.dropna(subset=[weight_var])
    #if not g.shape[0]:
    #    return pd.Series([np.NaN]*len(var))
    return pd.Series(np.average(g[var], axis=0, weights=g[weight_var]))
    # TODO: append the count of removed rows (NA in weight_var)

def raw_group_average(g, var, weight_var):
    gv = g[var].values
    gw = g[weight_var].values
    keep = ~np.isnan(gw)
    return np.average(gv[keep], axis=0, weights=gw[keep])

def npraw_group_average(gvalues, var_idx, weight_idx):
    #values = np.ma.array(g[var].values, mask=g[var].isna().values)
    gv = gvalues[:,var_idx]
    gw = gvalues[:,weight_idx]
    keep = ~np.isnan(gw)
    return np.average(gv[keep], axis=0, weights=gw[keep])

def npraw_group_average_w0(gvalues):
    """Vectorized weighted average from numpy input (2D) to numpy output (1D).
    
    From an array `gvalues` of dim (N, M) (lines, columns):
    - take the row weights in column 0;
    - ignore rows with NaN weight;
    - compute the weighted average of the other columns.
    """
    #values = np.ma.array(g[var].values, mask=g[var].isna().values)
    gv = gvalues[:,1:]
    gw = gvalues[:,0]
    keep = ~np.isnan(gw)
    return np.average(gv[keep], axis=0, weights=gw[keep])

def group_average(g, var, weight_var):
    #values = np.ma.array(g[var].values, mask=g[var].isna().values)
    gv = g[var].values
    gw = g[weight_var].values
    keep = ~np.isnan(gw)
    return pd.Series(np.average(gv[keep], axis=0, weights=gw[keep]))


def group_weighted_std(g, var, weight_var="median_brlen"):
    return pd.Series(weighted_std(g[var], axis=0, weights=g[weight_var]))

# for calculating rates: with dS, dN, t
def tree_dist_2_rate(g, dist_var, norm_var="median_brlen"):
    # in pandas, sum is done per columns.
    return g[dist_var].sum() / g[norm_var].sum()


def add_control_dates_lengths(ages, phyltree, control_ages_CI=None,
                              measures=['dist', 'dS', 'dN', 't'],
                              control_condition='really_robust & aberrant_dists == 0',
                              control_names=['timetree']  # 'dosreis'
                              ):
    # 'consecutive_zeros_dS==0'
    # Merge control dates
    ages_forcontrol = ages.query(control_condition)
    logger.info("%d nodes from robust trees", ages_forcontrol.shape[0])

    age_measures = ['age_'+m for m in measures]
    median_age_measures = ['median_age_'+m for m in measures]
    control_ages = ages_forcontrol[ages_forcontrol.type\
                                                       .isin(("spe", "leaf"))]\
                                   .groupby("taxon", sort=False)[age_measures]\
                                   .median()\
                                   .rename(columns=dict(zip(age_measures,
                                                            median_age_measures)))
    logger.debug('median_taxon_ages.columns = [%s]',
                 ' '.join(control_ages.columns))
    # FIXME: WHY does `median_taxon_ages` only have the 1st age_measure?
    #TODO: change for each measure as in `join_extra_ages`
    
    if control_ages_CI is not None:
        toconcat = (control_ages,)
        if 'timetree_age' not in control_ages_CI.columns:
            timetree_age = control_ages.index.to_series()\
                                             .apply(phyltree.ages.get)\
                                             .rename('timetree_age')  # age_timetree
            toconcat += (timetree_age,)

        toconcat += (control_ages_CI,)
        control_ages = pd.concat(toconcat, axis=1, sort=False)
    logger.debug('Control ages (columns): ' + ' '.join(control_ages.columns))

    # Add potentially missing ages (0 for extant species).
    for ctl in control_names:
        na_ages = control_ages['%s_age' % ctl].isna()
        if na_ages.any():
            for taxon in control_ages.index.values[na_ages]:
                # PhylTree ages are those from TimeTree
                if taxon in phyltree.listSpecies or ctl=='timetree':
                    control_ages.loc[taxon, '%s_age' % ctl] = phyltree.ages.get(taxon, np.NaN)

    #print(control_ages.sort_values('timetree_age', ascending=False))

    ages_controled = pd.merge(ages, control_ages,
                              left_on="taxon", right_index=True,
                              validate="many_to_one")

    # Merge control branch lengths
    invalid_taxon_parent = ages_controled.taxon_parent.isna()
    invalid_measures  = ages_controled[age_measures].isna().all(1)
    invalid_nodes = invalid_taxon_parent & ~invalid_measures
    # calibrated | branch_dS.isna()
    # Should be nodes whose parent node is the root.
    debug_columns = ['parent', 'subgenetree', 'taxon', 'taxon_parent'] \
                    + median_age_measures \
                    + ['branch_dS', 'age_dS']
    if invalid_nodes.any():
        #assert (ages_controled[invalid_taxon_parent].parent == \
        #        ages_controled[invalid_taxon_parent].root).all()
        logger.warning("%d invalid 'taxon_parent':head:\n%s\n"
                     "The following taxa have no parent taxa information, "
                     "please check:\n%s\n**DROPPING** this data!",
                       invalid_nodes.sum(),
                       ages_controled[invalid_nodes][
                           debug_columns
                           ].head(),
                       ', '.join(ages_controled[invalid_nodes].taxon.unique()))
        ages_controled.dropna(subset=['taxon_parent'], inplace=True)
    

    # Do not include length from/to duplication.
    ages_controled_spe2spe = ages_controled.dropna(subset=['taxon_parent']).query('type != "dup" & type_parent != "dup"')
    same_sp = ages_controled_spe2spe.taxon == ages_controled_spe2spe.taxon_parent
    if same_sp.any():
        logger.error("Failed to filter out duplications:\n%s",
                     ages_controled_spe2spe[same_sp].head(15))

    for ctl in control_names:
        ages_controled[ctl+'_brlen'] = \
            ages_controled_spe2spe.taxon_parent\
                                  .apply(control_ages[ctl+'_age'].get) \
            - ages_controled_spe2spe[ctl+'_age']

    #TODO: NaN for data not in ages_forcontrol?
    for m,am,mm in zip(measures, age_measures, median_age_measures):
        ages_controled['median_brlen_'+m] = \
            ages_controled_spe2spe.taxon_parent.apply(control_ages[mm].get) \
            - ages_controled_spe2spe[mm]
    #control_ages.reindex(ages_controled.taxon_parent)\
    #        .set_axis(ages_controled.index, inplace=False)
    # would be more 'Pandas-like'.

    # An approximation of the real branch length
    for m,am in zip(measures, age_measures):
        ages_controled['approx_%s' %m] = (ages_controled[am+'_parent']
                                            - ages_controled[am])

    # Resulting branch lengths
    branch_info = ["taxon_parent", "taxon"]
    #control_brlen = ages_controled.loc[
    #                    ~ages_controled.duplicated(branch_info),
    control_brlen = ages_controled.query('type != "dup" & type_parent != "dup"')\
                    .groupby(branch_info, sort=False)\
                    [['median_%s_%s' % (typ, m) for typ in ('age', 'brlen')
                        for m in measures]
                     + ["%s_%s" %(ctl, typ) for typ in ('age', 'brlen')
                         for ctl in control_names]]\
                    .first()#agg(lambda s: s[0])
                #.sort_values('timetree_age')
                    #].sort_values("taxon_parent", ascending=False)
                    #.reset_index(branch_info, drop=True)
                    #.set_axis(pd.MultiIndex.from_arrays(
                    #                control_brlen[branch_info].values.T,
                    #                names=branch_info),
                    #          inplace=False)\
                    #.drop(branch_info, axis=1)
    ###FIXME: Still some branches that are not in PhylTree (jump over taxa)
    display_html(control_brlen)
    return ages_controled, control_ages, control_brlen


def check_control_dates_lengths(control_brlen, phyltree, root,
                                measures=MEASURES):
    """Check NA values and if the branches fit the phylogenetic tree.
    
    Return: (unexpected_branches, lost_branches)"""
    get_phylchildren = lambda phyltree, ancdist: phyltree.items.get(ancdist[0], [])

    expected_branches, expected_dists = zip(*(((p[0],ch[0]),ch[1]) for p,ch in
                             dfw_pairs_generalized(phyltree,
                                                   get_phylchildren,
                                                   queue=[(None, (root,0))])))

    #logger.debug(expected_branches)
    #logger.debug(expected_dists)

    na_brlen = control_brlen.isna().any(axis=1)
    if na_brlen.any():
        print("MISSING branch lengths:\n" + str(control_brlen[na_brlen]))

    median_measures = ['median_brlen_%s' % m for m in measures]
    median_brlen_sum = control_brlen[median_measures].sum()
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

    median_treelen_phyltree = control_brlen.reindex(list(expected_branches))[median_measures].sum()
    timetree_treelen_phyltree = control_brlen.reindex(list(expected_branches)).timetree_brlen.sum()
    print("Sum of median branch lengths for branches found in phyltree =\n",
          str(median_treelen_phyltree).replace('\n', '\t\n'))
    print("Sum of timetree branch lengths for branches found in phyltree =",
          timetree_treelen_phyltree)
    return unexpected_branches, lost_branches


def compute_dating_errors(ages_controled, control='median', measures=['dS'],
                          rescale=None):
    """:param: `control` in median/timetree"""
    if measures is None:
        age_vars = ['age']
        brlen_vars = ['brlen']
    else:
        age_vars = ['age_'+m for m in measures]
        brlen_vars = ['brlen_'+m for m in measures]

    control += '_'
    prefix = '' if control=='median_' else control

    ### Inplace
    for age_var, brlen_var in zip(age_vars, brlen_vars):
        #ages_controled = ages_controled.assign(**{
        
        if control == 'median_':
            ctl_age, ctl_brlen = control+age_var, control+brlen_var
        else:
            # control in ('timetree_', 'dosreis_')
            ctl_age, ctl_brlen = control+'age', control+'brlen'

        logger.debug('Control=%r age_var=%r => subtract %r', control, age_var, ctl_age)

        ages_controled[prefix+"signed_dev_"+age_var] = (ages_controled[age_var]
                                                   - ages_controled[ctl_age])
        ages_controled[prefix+"signed_dev_"+brlen_var] = (ages_controled[age_var+'_parent']
                                                     - ages_controled[age_var]
                                                     - ages_controled[ctl_brlen])
        if rescale is None:
            ages_controled[prefix+"abs_dev_"+age_var] = ages_controled[prefix+"signed_dev_"+age_var].abs()
            ages_controled[prefix+"abs_dev_"+brlen_var] = ages_controled[prefix+"signed_dev_"+brlen_var].abs()
        elif rescale=='sym':
            # 1. Symetrical method
            raise NotImplementedError
        elif rescale=='asym':
            # 2. Asymetrical method. Let c be the real age, c0 and c1 the calibrations before and after, and x the estimation:
            # if x - c > 0 (the estimation is older) => normalize by c0 - c
            # if x - c < 0 (the estimation is younger) => normalize by c1 - c
            #ages_controled["sym_dev_" + age_var] = (ages_controled["signed_dev_" + age_var]
            #                                        /( - ages_controled[control+age_var]))
            raise NotImplementedError
        else:
            raise ValueError('`rescale` not in (None, "sym", "asym")')

    # Compute the mean only for nodes that were not calibrated.
    sgg = ages_controled.groupby("subgenetree", sort=False)
    #dev_measures = ['abs_age_dev', 'signed_age_dev', 'abs_brlen_dev', 'signed_brlen_dev']
    mean_errors = pd.concat((
                    sgg[[prefix+dev+age_var for age_var in age_vars
                            for dev in ('abs_dev_', 'signed_dev_')]].sum().div(
                        sgg['calibrated'].agg(lambda v: (1-v).sum()), axis=0),
                    sgg[[prefix+dev+brlen_var for brlen_var in brlen_vars
                            for dev in ('abs_dev_', 'signed_dev_')]].sum().div(
                        sgg[['calibrated', 'calibrated_parent']]\
                           .apply(lambda df: ((1-df.calibrated) | (1-df.calibrated_parent)).sum()), axis=0)
                    ),
                    axis=1)
    #TODO: median_errors

    return mean_errors


def join_extra_ages(new_ages_file, ages_data,
                    control_condition='really_robust & aberrant_dists == 0',
                    suffixes=('_extra', ''),
                    control_names=['timetree']):
                    #subset=None):
    """Re-use outputs of `analyse_age_errors`/`load_prepare_ages`/
    `add_control_dates_lengths`
    """
    extra_ages = pd.read_csv(new_ages_file, sep='\t', index_col=0)\
                   .rename(columns=lambda n: n.replace('.', '_'))\
                   .rename_axis(index='name')\
                   .reset_index()\
                   .assign(subgenetree=lambda df: df.subgenetree.fillna(method='ffill'))
    age_measures = [colname for colname in extra_ages.columns
                        if colname.startswith('age_')]

    logger.info('Joining new age measures: ' + ' '.join(age_measures))

    # 1. Needs to concat the parent data!! (for compute_mean_errors of 'brlen')
    extra_ages = pd.merge(extra_ages,
                          extra_ages.loc[extra_ages.type!='leaf',
                                         ['name', 'subgenetree'] + age_measures],
                      how="left",
                      left_on=["parent", "subgenetree"],
                      right_on=["name", "subgenetree"],
                      suffixes=('', '_parent'),
                      indicator=True, validate='many_to_one')


    # 2. Join to the previous data.
    ages = pd.merge(
               extra_ages.drop(columns=['calibrated', 'type', 'parent', 'taxon']),
               ages_data.ages_controled_withnonrobust.reset_index(),
               'inner',  # to remove what has been removed (outgroups)
               on=['name', 'subgenetree'],
               suffixes=suffixes,
               sort=False)\
           .set_index('name')
    logger.debug('New ages shape: %s', ages.shape)
    logger.debug('ages.columns = %s', ' '.join(ages.columns))

    if set(ages_data.ages_controled.columns).intersection(age_measures):
        logger.warning('New ages already in the data: %s (joined with suffix %r)',
                       age_measures, suffixes[0])
        ages.rename(columns={am+'_parent'+suffixes[0]: am+suffixes[0]+'_parent'
                             for am in age_measures
                             if am in ages_data.ages_controled.columns},
                    inplace=True)
        #logger.debug('extra_ages.columns = %s', ' '.join(extra_ages.columns))
        age_measures = [am + suffixes[0]
                        if am in ages_data.ages_controled.columns 
                        else am
                        for am in age_measures]
    measures=[am[4:] for am in age_measures]

    # 3. Rerun a small part of `add_control_dates_lengths`
    ages_forcontrol = ages.query(control_condition)
    logger.info("%d nodes from robust trees", ages_forcontrol.shape[0])

    # 3.1. Recompute the median ages.
    median_taxon_ages = ages_forcontrol[ages_forcontrol.type\
                                                       .isin(("spe", "leaf"))]\
                                   .groupby("taxon", sort=False)[age_measures].median()\
                                   .rename(columns={am: 'median_'+am
                                                    for am in age_measures})
    control_ages = pd.concat((median_taxon_ages, ages_data.control_ages),
                             axis=1, sort=False)

    #print(control_ages.sort_values('timetree_age', ascending=False))

    ages_controled = pd.merge(ages,
                              control_ages[['median_'+am for am in age_measures]],
                              left_on="taxon", right_index=True,
                              validate="many_to_one")

    # 3.2. Check duplication nodes (optional?)
    ages_controled_spe2spe = ages_controled.dropna(subset=['taxon_parent']).query('type != "dup" & type_parent != "dup"')
    same_sp = ages_controled_spe2spe.taxon == ages_controled_spe2spe.taxon_parent
    if same_sp.any():
        logger.error("Failed to filter out duplications:\n%s",
                     ages_controled_spe2spe[same_sp].head(15))

    # 3.3. Median branch lengths from the given measures.
    ages_controled = ages_controled.assign(**{
        ('median_brlen_%s' % m):
            (ages_controled_spe2spe.taxon_parent\
                        .apply(control_ages['median_%s' % am].get)
             - ages_controled_spe2spe['median_%s' % am])
        for m,am in zip(measures, age_measures)})

    branch_info = ["taxon_parent", "taxon"]
    #control_brlen = ages_controled.loc[
    #                    ~ages_controled.duplicated(branch_info),
    control_brlen = ages_data.control_brlen.join(
                        ages_controled.query('type != "dup" & type_parent != "dup"')\
                        .groupby(branch_info, sort=False)\
                        [['median_%s' % am for am in age_measures]]\
                        .first())
    #check_control_dates_lengths(control_brlen, phyltree, root, measures)

    ages_controled_cond = ages_controled.query(control_condition).copy(deep=False)
    # copy to silence the `SettingWithCopyWarning`. Alternatively, set `_is_copy=None`
    mean_errors = pd.concat(
                    [compute_dating_errors(ages_controled_cond, ctl, measures)
                     for ctl in ['median']+control_names],
                    axis=1, sort=False, copy=False)

    return age_analysis_data(ages_controled_cond,
                             ages_controled,
                             None,  # ns  #TODO: add a dict of functions to aggregate columns.
                             control_ages,
                             control_brlen,
                             mean_errors)


def lineage_evolutionary_rates(anc, ages_file,
                               phyltree,
                               control_ages_CI,
                               stats_tmpl='subtreesGoodQualO2_%sstats-%s.tsv',
                               control_condition='really_robust & aberrant_dists == 0',
                               measures=MEASURES,
                               control_names=['timetree'],  # 'dosreis'
                               saveas=None):
    # Also see `analyse_age_errors`
    ts = load_stats_tree(stats_tmpl % ('tree', anc))
    aS = load_stats_al(stats_tmpl % ('al', anc))

    dist_measures = ['branch_%s' % m for m in measures]

    ages, ns = load_prepare_ages(ages_file, ts, measures)
    ages = pd.merge(ages, aS.join(ns, how='outer'),
                    left_on='subgenetree', right_index=True, sort=False,
                    validate='many_to_one')
                     
    ages_controled_withnonrobust, control_ages, control_brlen =\
        add_control_dates_lengths(ages, phyltree, control_ages_CI, measures,
                                  control_condition, control_names)

    ages_controled = ages_controled_withnonrobust.query(control_condition).copy(deep=False)
    unexpected_branches, lost_branches = check_control_dates_lengths(
                                                control_brlen, phyltree, anc,
                                                measures)
    control_brlen.drop(unexpected_branches, inplace=True)
    age_analysis = age_analysis_data(ages_controled,
                            ages_controled_withnonrobust,
                            ns, control_ages, control_brlen, mean_errors=None)
    return (age_analysis,) + lineage_evolutionary_rates_fromdata(age_analysis, dist_measures, control_names)


def lineage_evolutionary_rates_fromdata(age_analysis, dist_measures,
                                        control_names=['timetree']):
    def wmean(g, wkey='ingroup_glob_len'):
        gdata = g.drop(columns=wkey)
        return pd.Series(np.average(gdata, axis=0, weights=g[wkey]), index=gdata.columns)
        #return np.average(gdata, axis=0, weights=g[wkey])

    #def func(g, wkey='ingroup_glob_len'):
    #    gdata = g.drop(columns=wkey)
    #    return pd.DataFrame(
    #                np.average(gdata, axis=0, weights=g[wkey]),
    #                     index=gdata.columns)

    lineage_groups = age_analysis.ages_controled.groupby(['taxon_parent', 'taxon'], sort=False)
    lineage_brlen = lineage_groups[dist_measures]\
            .agg(['median', 'mean', 'std'])\
            .join(lineage_groups[dist_measures + ['ingroup_glob_len']]\
                        .apply(wmean)\
                        .set_axis(pd.MultiIndex.from_product([dist_measures, ['wmean']]),
                            axis=1, inplace=False),
                  sort=False)\
            .sort_index(1)
    ctl = control_names[0]

    print('lineage_brlen =\n', lineage_brlen.head().iloc[:,:4])
    print('control_brlen =\n', age_analysis.control_brlen.head().iloc[:,:4])
    lineage_rates = lineage_brlen.div(age_analysis.control_brlen[ctl+'_brlen'],
                                      axis='index')

    branch_na_rates = lineage_rates.isna().any(axis=1)
    if branch_na_rates.any():
        logger.info('DROP NA rates rows:\n%s', lineage_rates.index[branch_na_rates].values)
        lineage_rates.dropna(inplace=True)

    return lineage_rates, lineage_brlen


def display_lineage_evolutionary_rates(lineage_rates, lineage_brlen,
                                       age_analysis, phyltree, anc,
                                       measures=MEASURES,
                                       control='timetree',
                                       figsize=None,
                                       extra_text='',
                                       cmap_minmax=None):

    ages_controled = age_analysis.ages_controled
    control_brlen = age_analysis.control_brlen
    dist_measures = ['branch_%s' % m for m in measures]
    rate_measures = ['%s_rate' % m for m in measures]

    ordered_branches = ['%s--%s' % br for br in dfw_pairs(phyltree, queue=[(None, anc)], closest_first=True)]

    ordered_branches_bylen = ['%s--%s' % v
                              for v in control_brlen.sort_values(control+'_brlen').index]

    # Table of summary values
    styled_rates = lineage_rates\
            .style.background_gradient(cmap='PRGn')#,
                                       #subset=[(d, stat)
                                       #    for d in dist_measures
                                       #    for stat in ('median', 'mean', 'std')])
    display_html(styled_rates)
    outputs = [styled_rates]

    m = measures[0]
    br_m = 'branch_'+m
    rates = lineage_rates[(br_m, 'median')]
    print('Min of median %s rate: %s (%s)' %(m, rates.min(), rates.idxmin()))
    print('Max of median %s rate: %s (%s)'%(m, rates.max(), rates.idxmax()))
    
    # Violin plot of rates
    ages_data = ages_controled\
                    .assign(
                        branchtaxa=(ages_controled.taxon_parent + '--' + ages_controled.taxon)
                    )\
                    .join(ages_controled[dist_measures].div(
                                ages_controled.timetree_brlen, axis=0),
                          how='outer', rsuffix='_rate')

    ax = sb.violinplot('branchtaxa', br_m+'_rate',
                       data=ages_data,
                       width=1, order=ordered_branches_bylen)
    
    values = ages_data.groupby('branchtaxa', sort=False)[br_m+'_rate']
    #rate_q = values.quantile([0.01, 0.99]).values
    #rate_range = rate_q[1] - rate_q[0]
    #ylim = (rate_q[0] - 0.01*rate_range, rate_q[0] + 0.01*rate_range)
    ndev = 2
    get_ylow = lambda v: max(v.min(), v.mean() - ndev*v.std())
    get_yup = lambda v: min(v.max(), v.mean() + ndev*v.std())

    ylim = (values.apply(get_ylow).min(),
            values.apply(get_yup).max())

    logger.debug('Setting new ylim: %s', ylim)
    ax.set_ylim(ylim)
    ax.set_ylabel(m+' rate')
    plt.setp(ax.get_xticklabels(), rotation=45, va='top', ha='right')
    ax.set_title('Variation of synonymous substitutions/site/My' +
                 (('. '+extra_text) if extra_text else ''))
    outputs.append(ax.figure)

    # Tree with branch length in median number of substitutions
    fig, ax = plt.subplots(figsize=figsize)
    plottree(phyltree,
             get_items=(lambda t,nd: [(ch,
                                       lineage_brlen.loc[(nd[0], ch),
                                                         (br_m, 'median')])
                                     for ch,d in t.items.get(nd[0], [])]),
             get_label=(lambda t,n: n),
             root=anc, ax=ax)
    ax.set_xlabel('Median %s branch length' % m)
    ax.set_title(extra_text)
    outputs.append(fig)

    # Tree with colored branches
    fig, (cax,ax) = plt.subplots(1, 2, figsize=figsize,
                                 gridspec_kw={'width_ratios': [1,20]})
    cax.set_position(cax.get_position().anchored('SE').shrunk(1, 0.5).shrunk_to_aspect(20))
    
    lines, coords, subaxes = plottree(phyltree,
                 lambda t,nd: [(ch, control_brlen.loc[(nd[0], ch), control+'_brlen'])
                                 for ch,d in t.items.get(nd[0], [])],
                 phyltree_methods.get_label,
                 root=anc,
                 edge_colors=lineage_rates.reset_index('taxon_parent', drop=True)[
                             (br_m, 'median')],
                 edge_cmap='afmhot', add_edge_axes=None, style='squared',
                 ax=ax)
    #cax = ax.inset_axes((0, 0.5, 0.05, 0.5))
    cbar = fig.colorbar(lines, cax=cax)
    cbar.set_label('Median '+br_m)
    cax.yaxis.set_ticks_position('left')
    cax.yaxis.set_label_position('left')
    ax.set_title('Median evolutionary %s rate across %s descendants' % (measures[0], anc)
                 + (('. '+extra_text) if extra_text else ''))
    ax.set_xlabel('Age (My)') #+ control
    outputs.append(fig)
    #fig.tight_layout()
    return outputs


def compute_branchrate_std(ages_controled, dist_measures,
                           branchtime='median_brlen_dS', taxon_age=None,
                           mean_condition=None,
                           std_condition=None,
                           weighted=False, poisson=True):
    """
    Example filter_condition for approximated dS (nonlocal):
    '(calibrated==1) & (calibrated_parent==1)'
    """

    groupby_cols = ["subgenetree", "taxon_parent", "taxon",
                    branchtime] + dist_measures
    if taxon_age is not None:
        groupby_cols.append(taxon_age)  ## "median_age"/"median_age_dS"
    #ages_controled["omega"] = ages_controled.branch_dN / ages_controled.branch_dS

    if mean_condition:
        ages_controled = ages_controled.query(mean_condition)
        if not ages_controled.shape[0]:
            logger.warning('Selected 0 rows with `mean_condition`.')
    elif not ages_controled.shape[0]:
        logger.warning('0 rows in data')
    sgg = subgenetree_groups = ages_controled[groupby_cols].groupby('subgenetree', sort=False)

    ### Average (substitution) rates over the tree:
    #     sum of all branch values / sum of branch lengths

    # Sum aggregation + division broadcasted on columns
    # This is the maximum likelihood estimate for a constant Poisson rate.
    if poisson:
        cs_rates = sgg[dist_measures].sum().div(sgg[branchtime].sum(), axis=0)
    #    Problem: directly proportional to the tree length. Because branch times same across genes
    elif not weighted:
        cs_rates = ages_controled[dist_measures]\
                    .div(ages_controled[branchtime], axis=0)\
                    .join(ages_controled[['subgenetree']], sort=False)\
                    .groupby('subgenetree', sort=False)\
                    .mean()
    else:
        # Will be slow.
        cs_rates = sgg.apply(lambda g: g[dist_measures].mul(g[branchtime], axis=0))\
                        .div(sgg[branchtime].sum()**2, axis=0)
    #cs_rates["omega"] = (sgg.branch_dN / sgg.branch_dS).apply()

    # This weird operation fills the subgenetree-mean-rate into each node name.
    aligned_csrates = ages_controled[['subgenetree']]\
                            .join(cs_rates, on='subgenetree', sort=False)\
                            .drop(columns='subgenetree')

    ### Weighted standard deviation of substitution rates among branches
    # the squared deviations from the mean rate:
    rate_dev = (ages_controled[dist_measures]\
                    .div(ages_controled[branchtime], axis=0)\
                    .sub(aligned_csrates, axis=0)**2)\
                .join(ages_controled[['subgenetree', branchtime]], sort=False)
    #logger.debug('rate_dev.shape = %s; columns = %s', rate_dev.shape, rate_dev.columns)

    cs_rates.rename(columns=lambda d: d.replace('branch_', '')+'_rate', inplace=True)
    rate_measures = cs_rates.columns.tolist()  #[(m.replace('branch_', '') + '_rate') for m in dist_measures]

    #tmp = pd.merge(ages_controled[["subgenetree", branchtime] + dist_measures],
    #               cs_rates, left_on="subgenetree", right_index=True, sort=False)

    ### subtract branch rate with mean rate, then square.
    #rate_dev_dict = {d: (tmp[d] / tmp[branchtime] - tmp[r])**2
    #                 for d,r in zip(dist_measures, rate_measures)}
    
    #rate_dev = pd.DataFrame(rate_dev_dict)\
    #                .join(tmp[['subgenetree', branchtime]], sort=False)

    # Caching the column indices to apply raw numpy computations (faster).
    #dist_measures_idx = rate_dev.columns.get_indexer(dist_measures)
    #branchtime_idx = rate_dev.columns.tolist().index(branchtime)
    rsgg = rate_dev.groupby("subgenetree", sort=False)
    #rsgg_g0name = list(rsgg.groups.keys())[0]
    #rsgg_g0 = rsgg.get_group(rsgg_g0name)
    #logger.debug('group 0 "%s" shape = %s; columns = %s; name = %s',
    #             rsgg_g0name, rsgg_g0.shape, rsgg_g0.columns, getattr(rsgg_g0, 'name', None))
    #logger.debug('group 0 mean rate =\n%s;\nvalues = %s',
    #             cs_rates.loc[rsgg_g0name], rsgg_g0.values)
    #logger.debug('group 0 group_average = %s', npraw_group_average_w0(rsgg_g0.values))

    if weighted:
        cs_stds = rsgg[[branchtime] + dist_measures]\
                        .apply(lambda g: pd.Series(np.sqrt(
                                                    npraw_group_average_w0(g.values)),
                                                  name=g.name)
                        )
    else:
        cs_stds = np.sqrt(rsgg[dist_measures].mean())
    cs_stds.set_axis([r+'_std' for r in rate_measures], axis=1, inplace=True)

    #cs_stds = pd.DataFrame(
    #            np.stack(
    #                rsgg.apply(
    #                (lambda x:
    #                        sqrt(
    #                            npraw_group_average_w0(x.values)
    #                            )
    #                        )
    #                )
    #            ),
    #            columns=[r+'_std' for r in rate_measures])

    cs_rates = pd.concat((cs_rates, cs_stds), axis=1, sort=False, verify_integrity=True)
    # Checks
    inf_cols = np.isinf(cs_rates).any(axis=0)
    inf_rows = np.isinf(cs_rates).any(axis=1)
    if inf_rows.any():
        logger.warning('%d Inf values in columns %s. DROPPING rows!',
                       inf_rows.sum(),
                       cs_rates.columns[inf_cols].tolist())
        cs_rates = cs_rates[~inf_rows].copy(deep=False)
    return cs_rates


def triplet_aggfunc(func, func_args, func_kwargs, parentgrouped, row):
    """To use within `apply` only! (because of .name attribute)"""
    try:
        children_rows = parentgrouped.get_group(row.name)
    except KeyError:
        #logger.debug('parent %r not found', row.name)
        return pd.Series(np.NaN, index=row.index)
    return pd.concat((row, children_rows), sort=False).agg(func, *func_args, **func_kwargs)


def raw_triplet_std(row, parentgrouped):
    try:
        children_rows = parentgrouped.get_group(row.name)
    except KeyError:
        #logger.debug('parent %r not found', row.name)
        return np.full(row.shape[0], np.NaN)
    return np.nanstd(np.vstack((row, children_rows)), axis=0, ddof=0)


def compute_correlated_rate(ages_controled, dist_measures,
                            branchtime='median_brlen_dS', taxon_age=None,
                            mean_condition=None,
                            std_condition=None):
    rates = ages_controled[dist_measures].div(ages_controled[branchtime], axis=0)
    sister_rates = rates.groupby(ages_controled.parent, sort=False)
    #compute_triplet_std = partial(triplet_aggfunc, 'std', (), {'ddof': 0}, sister_rates)
    #compute_triplet_std = lambda row: pd.Series(raw_triplet_std(row, sister_rates))
    #logger.debug('rates: type=%s; shape=%s; columns=%s', type(rates), rates.shape, rates.columns)
    #logger.debug('rates.head() =\n%s', rates.head(20))

    #TODO: subset ages_controled where 'type!="leaf"'
    triplet_rate_corrstds = rates.apply(raw_triplet_std, axis=1, raw=False,
                                    result_type='expand', args=(sister_rates,))

    groupby_cols = ["subgenetree", "taxon_parent", "taxon",
                    branchtime] + dist_measures
    if taxon_age is not None:
        groupby_cols.append(taxon_age)  ## "median_age"/"median_age_dS"

    if mean_condition:
        triplet_rate_corrstds = triplet_rate_corrstds.loc[ages_controled.query(mean_condition).index]
        if not triplet_rate_corrstds.shape[0]:
            logger.warning('Selected 0 rows with `mean_condition`.')
    elif not triplet_rate_corrstds.shape[0]:
        logger.warning('0 rows in data')
    
    #return triplet_rate_corrstds
    return triplet_rate_corrstds.groupby(ages['subgenetree'], sort=False).mean()


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
    # WARNING: unused?
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
                                           split=True, order=None, cut=0,
                                           points=False,
                                           **kwargs):
    Q_col = "Q_" + criterion_name
    if isin is None:
        # Look at extreme quantiles only
        isin = (annot_df[Q_col].min(), annot_df[Q_col].max())
        logger.info(isin)
        
    ax = sb.violinplot(x="taxon", y="age_dS", hue=Q_col,
                       data=annot_df[(annot_df.type == "spe")
                                     & annot_df[Q_col].isin(isin)],
                       order=order,
                       split=split,
                       cut=cut,
                       **kwargs)
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, va='top', ha='right')
    if points:
        #violins = ax.get_children()
        orig_legend_args = ax.get_legend_handles_labels()
        orig_legend_title = ax.legend_.get_title().get_text()

        sb.pointplot(x="taxon", y="age_dS", hue=Q_col,
                     data=annot_df[(annot_df.type == "spe") & annot_df[Q_col].isin(isin)],
                     order=order,
                     estimator=np.nanmedian,
                     #ci='sd',
                     join=False,
                     dodge=(1./(len(isin)+1)),
                     palette=['#b8b8b8'],
                     alpha=0.8,
                     ax=ax,
                     legend=False)
        # Reset the legend.
        #print(ax.get_children())
        #plt.setp(set(ax.get_children()).difference(violins), label='')
        ax.legend_ = None
        ax.legend(*orig_legend_args, title=orig_legend_title)
    return ax


def violin_spe_ages_vs_criterion(ages, criterion_serie, criterion_name=None,
                                 nquantiles=10, split=False, order=None, **kwargs):
    criterion_name = criterion_name or criterion_serie.name
    annot_ages = annot_quantiles_on_criterion(ages, criterion_serie,
                                              criterion_name, nquantiles)
    isin = None if split else list(range(nquantiles))
    return _violin_spe_ages_vs_criterion_quantiles(
                                annot_ages[annot_ages.taxon != 'Simiiformes'],
                                criterion_name, isin, split, order, **kwargs)


class age_analysis_data(object):
    """Just a namespace to hold several objects of one analysis."""
    def __init__(self,
                 ages_controled,
                 ages_controled_withnonrobust,
                 ns,
                 control_ages,
                 control_brlen,
                 mean_errors):
        for k,v in locals().items():
            if k != 'self':
                setattr(self, k, v)


def analyse_age_errors(ages_file, root, phyltree, control_ages_CI, ts,
                       measures=['dS', 'dist', 'dN', 't'],
                       control_names=['timetree']):  # 'dosreis'
    #, aS, cs, clS=None
    ages, ns = load_prepare_ages(ages_file, ts, measures)
    ages_controled_withnonrobust, control_ages, control_brlen =\
        add_control_dates_lengths(ages, phyltree, control_ages_CI, measures,
                                  control_names=control_names)
    ages_controled = ages_controled_withnonrobust.query('really_robust & aberrant_dists == 0').copy(deep=False)
    unexpected_branches, lost_branches = check_control_dates_lengths(
                                                control_brlen, phyltree, root,
                                                measures)
    control_brlen.drop(unexpected_branches, inplace=True)
    # Rates
    #cs_rates = compute_branchrate_std(ages_controled, dist_measures)
    #cs_rates_approx
    #cs_rates_onlyAnc
    #cs_rates_withoutAnc

    mean_errors = pd.concat(
                [compute_dating_errors(ages_controled, ctl, measures)
                 for ctl in ['median']+control_names],
                axis=1, sort=False, copy=False)
    
    return age_analysis_data(ages_controled, ages_controled_withnonrobust, ns,
                             control_ages, control_brlen, mean_errors)


def lassoselect_and_refit(a, y, features, atol=1e-3, method='elastic_net',
                          L1_wt=1, cov_type='HC1', **psummary_kw):
    exog = sm.add_constant(a[features])
    if 'const' not in exog:
        logger.warning('No constant added to `exog`: some features are already constant')
    ols = sm.OLS(a[y], exog)
    fitlasso = ols.fit_regularized(method=method,
                                   L1_wt=L1_wt,  # lasso)
                                   start_params=None)

    NaN_coefs = fitlasso.params.isna()
    if NaN_coefs.all():
        raise RuntimeError("ALL coefs estimated as NaN.")
    elif NaN_coefs.any():
        logger.warning("Some coefs estimated as NaN.")

    print('Non zeros params: %d/%d (from %d input features).\n'
          'Almost zeros params (<%g): %d.' % (
            (fitlasso.params.drop('const') != 0).sum(),
             fitlasso.params.drop('const').shape[0],
             len(features),
             atol,
             np.isclose(fitlasso.params.drop('const'), 0, atol=atol).sum()))

    sorted_coefs = fitlasso.params.drop('const').abs().sort_values(ascending=False)
    sorted_features = sorted_coefs.index.tolist()
    multicol_test_cumul = pd.DataFrame(columns=['cumul_colinearity',
                                                'cumul_colinearity_noconst'],
                               index=['const'] + sorted_features,
                               dtype=float)
    multicol_test_cumul.loc['const', 'cumul_colinearity'] = 1
    for i, ft in enumerate(sorted_features, start=1):
        multicol_test_cumul.loc[ft] = [multicol_test(exog[['const'] + sorted_features[:i]]),
                                       multicol_test(exog[sorted_features[:i]])]

    pslopes = sm_pretty_summary(fitlasso, multicol_test_cumul)
    display_html(pslopes)

    #scatter_density('ingroup_glob_len', y, data=a, alpha=0.5);
    #sb.violinplot('triplet_zeros_dS', 'abs_age_dev', data=a);

    print('\n#### OLS refit of Lasso-selected variables (%g)' % atol)
    #TODO: delete but make a refit.
    selected_features = sorted_coefs[sorted_coefs > atol].index.tolist()
    fit = sm.OLS(a[y], exog[['const'] + selected_features]).fit(cov_type='HC1')

    bars_to_show = ['coef', 'Lasso coef', 'Simple regression coef']
    bars_to_show += psummary_kw.pop('bars', [])
    info_to_join = psummary_kw.pop('join', None)
    param_info = pd.concat((
                    pslopes.data.coef.rename('Lasso coef'),
                    pd.Series({ft: sm.OLS(a[y],
                                          sm.add_constant(a[ft])
                                          ).fit().params[ft]
                               for ft in features
                              }, name='Simple regression coef'),
                    ) + (() if info_to_join is None else (info_to_join,)),
                    axis=1, sort=False)

    preslopes = sm_pretty_summary(fit, param_info, bars=bars_to_show, **psummary_kw)
    display_html(preslopes)

    assert set(selected_features) == set(preslopes.data.drop('const').index)
    # Test of homoscedasticity
    #sms.linear_harvey_collier(fit)
    return fitlasso, fit, pslopes, preslopes


_must_transform = dict(
        ingroup_nucl_entropy_median=binarize,
        #ingroup_nucl_parsimony_median=binarize, #constant
        ingroup_codon_entropy_median=binarize,
        ingroup_codon_parsimony_median=binarize,
        rebuilt_topo=binarize,
        consecutive_zeros=binarize, # - triplet_zeros
        sister_zeros=binarize,      # - triplet_zeros
        triplet_zeros=binarize,
        bootstrap_min=notransform,
        prop_splitseq=binarize,
        convergence_warning=binarize,  # notransform
        consecutive_zeros_dS=binarize,
        sister_zeros_dS=binarize,
        triplet_zeros_dS=binarize,
        consecutive_zeros_dN=binarize,
        sister_zeros_dN=binarize,
        triplet_zeros_dN=binarize,
        r2t_dN_mean=make_best_logtransform,
        gb_Nblocks=notransform,
        hmmc_propseqs=notransform,
        freq_null_dS=binarize,
        null_dist_before=binarize,
        null_dS_before=binarize,
        null_dN_before=binarize,
        null_dist_after=binarize,
        null_dS_after=binarize,
        null_dN_after=binarize,
        #null_dN_before=sqrt
        **{'dN_rate'+('_'+setting if setting else ''): make_best_logtransform
           for setting in ('', 'global', 'local', 'nonlocal', 'global_approx',
                           'local_approx', 'global_beastS')#self.rate_settings
          }
       )

# Example parametrisation (seemed reasonable for most regressions so far)
_must_drop_features = ["ls", "seconds",  # ~ ingroup_glob_len
                       "ingroup_std_gaps", # ~ ingroup_std_len
                       "dS_treelen",       # ~dS_rate
                       "dN_treelen",
                       "treelen",
                       "ingroup_codon_entropy_mean", # ~ ingroup_codon_parsimony_mean
                       "ingroup_codon_entropy_std",  # 
                       "ingroup_nucl_entropy_std",   # ~ ingroup_codon_parsimony_std
                       "ingroup_nucl_entropy_mean",
                       "ingroup_nucl_entropy_median",  # ~ ingroup_codon_entropy_median
                       "ingroup_nucl_parsimony_mean",
                       "ingroup_nucl_parsimony_std",
                       "r2t_t_mean", "r2t_dS_mean", "r2t_dN_mean",
                       "r2t_t_std",  "r2t_dS_std",  "r2t_dN_std",
                       "bootstrap_mean",  # ~ bootstrap_min
                       "brOmega_skew",
                       "ingroup_std_N",
                       "ingroup_std_CpG",  # ~ ingroup_std_GC
                       "NnonsynSites",  # ~ ls/3 - NsynSites
                       # decorrelated:
                       "ingroup_mean_CpG",
                       # "ingroup_codon_parsimony_std",
                       # "NnonsynSites", "Nsynsites", "brOmega_std",
                       # "ingroup_mean_CpG", "ingroup_std_N", "lnL",
                       # "dS_rate_std", "t_rate_std", "dN_rate_std", "dist_rate_std"
                       # BeastS
                       "treeL_12_med", #~mean
                       "birthRateY_med",
                       "ucldMean_12_med",
                       "ucldStdev_12_med",
                       "ucldMean_3_med",
                       "TreeHeight_med",
                       "rate_12_mean_med",
                       "gammaShape_med",
                       "rateAG_med"]

# FIXME: only if the variable was log-transformed!
_must_renormlogdecorr = (
    [('brOmega_std', 'brOmega_mean')]
     #('ingroup_std_N', 'ingroup_mean_N')]
      # Normalise the rate deviations by the rate mean.
    + [('%s_rate_std%s' %(m, ('_'+setting if setting else '')),
        '%s_rate%s' %(m, ('_'+setting if setting else '')))
       for setting in ('', 'global', 'local', 'nonlocal', 'global_approx',
                       'local_approx', 'global_beastS')#self.rate_settings
       for m in MEASURES]
    + [('treeL_3_stddev', 'treeL_3_med'),
       ('treeL_12_stddev', 'treeL_12_med')],
    dict(RsynSites=('NsynSites',     'ls'),
         sitelnL=('lnL', 'ingroup_glob_len')))

_must_renormdecorr = (
    [('ingroup_codon_parsimony_std', 'ingroup_codon_parsimony_mean')],
    {})

_must_logdecorr = (
        #[('%s_zeros%s' %(how, ('' if m=='dist' else '_'+m)),
        #  'triplet_zeros'+('' if m=='dist' else '_'+m))
        # for m in MEASURES
        # for how in ('sister', 'consecutive')],
        [],  # No need anymore since this is now handled in treestats.
        {})

# variable name, variable values. Dropped by .isin()
# Must be the decorr variable name.
_must_drop_data = dict(prop_splitseq=(1,),
                       **{'null_%s_%s' % (m,where): (1,)
                          for m in MEASURES
                          for where in ('before', 'after')},
                       **{'%s_zeros%s' %(what,('' if m=='dist' else '_'+m)): (1,)
                          for m in MEASURES
                          for what in ('triplet', 'Rsister', 'sister',
                                       'Rconsecutive', 'consecutive')})

_protected_features = ['RdS_rate_std', 'dS_rate_std',
                       'RdS_rate_std_local', 'dS_rate_std_local',
                       'RdS_rate_std_nonlocal', 'dS_rate_std_nonlocal',
                       'RbeastS_rate_std', 'beastS_rate_std',
                       'ingroup_glob_len']
# anc = 'Catarrhini'
# param = 'um1.new'
# nsCata = age_analyses[anc][param].ns
class full_dating_regression(object):

    init_vars = ['data',
                 'same_alls',
                 'dataset_params',
                 'responses',
                 'features',
                 'measures',  # measures of branch lengths and ages.
                 'ref_suggested_transform',
                 'impose_transform',
                 'to_renormdecorr',
                 'to_renormlogdecorr',
                 'to_subtract',
                 'must_drop_features',
                 'protected_features',
                 'must_drop_data']

    init_defaults = {'ref_suggested_transform': dict,
                     'impose_transform': dict,
                     'must_drop_features': list,
                     'to_renormdecorr': lambda: ([], {}),
                     'to_renormlogdecorr': lambda: ([], {}),
                     'to_subtract': lambda: ([], {}),
                     'protected_features': list,
                     'must_drop_data': dict}

    def __init__(self, data, same_alls, dataset_params, responses, features,
                 measures=MEASURES, ref_suggested_transform=None,
                 impose_transform=None, must_drop_features=None,
                 to_renormdecorr=None, to_renormlogdecorr=None, to_subtract=None,
                 protected_features=None, must_drop_data=None):
        for k,v in locals().items():
            if k != 'self':
                if k in self.init_defaults and v is None:
                    v = self.init_defaults[k]()  # initialize to the proper type.
                setattr(self, k, v)
        self.displayed = []  # Figures and styled dfs

    @classmethod
    def from_other(cls, other_regression):
        """Instantiate by copying all attributes from another instance."""
        self = cls(**{ivar: getattr(other_regression, ivar)
                               for ivar in cls.init_vars})
        #vars() doesn't contain properties nor private attributes
        for k,v in vars(other_regression).items():
            if k not in cls.init_vars:
                try:
                    setattr(self, k, v)
                except AttributeError as err:
                    err.args += (k,)
                    raise
        return self

    def do_rates(self, unnamed_rate_setting=None, **named_rate_settings):
        # Compute cs_rates
        dist_measures = ['branch_'+m for m in self.measures]
        if unnamed_rate_setting is not None:
            named_rate_settings[''] = unnamed_rate_setting
        self.rate_settings = named_rate_settings

        toconcat = []
        for setting, rate_args in named_rate_settings.items():
            print('### Compute rates with setting %r and measures %s.'
                    % (setting, ','.join(self.measures)))
            kwargs = {'branchtime': 'median_brlen_dS',
                      'taxon_age': 'median_age_dS',  # defaults kwargs
                      **rate_args}
        #display_html('<h3>Compute rates</h3>', raw=True)
            cs_rates = compute_branchrate_std(self.data.ages_controled,
                                              dist_measures, **kwargs)
            if setting:
                cs_rates.rename(columns=lambda c: c+'_'+setting, inplace=True)
            toconcat.append(cs_rates)

        self.cs_rates = pd.concat(toconcat, axis=1,
                                  sort=False, verify_integrity=True)
        #self.features.extend(self.cs_rates.columns)

    @property  #@dependency
    def suggested_transform(self):
        try:
            return self._suggested_transform
        except AttributeError:
            responses, features = self.responses, self.features
            ref_suggested_transform = self.ref_suggested_transform
            impose_transform = self.impose_transform
            alls = self.alls

            suggested_transform = test_transforms(alls,
                                    [ft for ft in responses+features
                                        if ft not in impose_transform]) #ages_features + rate_features
            #TODO: add to self.displayed

            suggested_transform.update(**dict(
                                        (ft, t(alls[ft]))
                                        if (t is make_best_logtransform and ft in alls)
                                        else (ft,t)
                                        for ft,t in impose_transform.items())
                                      )

            if ref_suggested_transform:
                onlyref = []
                onlynew = []
                diff_funcs = []
                diff_args = []
                for k in set(ref_suggested_transform).union(suggested_transform):
                    if k not in ref_suggested_transform:
                        onlynew.append(k)
                    elif k not in suggested_transform:
                        onlyref.append(k)
                    else:
                        reffunc = ref_suggested_transform[k].__name__
                        newfunc = suggested_transform[k].__name__
                        if reffunc != newfunc:
                            reffunc0, *refargs = reffunc.split('(')
                            newfunc0, *newargs = newfunc.split('(')
                            if reffunc0 == newfunc0:
                                newargs = newargs[0] if newargs else ''
                                refargs = refargs[0] if refargs else ''
                                diff_args.append((k, '%s VS %s' % (newargs.rstrip(')'),
                                    refargs.rstrip(')'))))
                            else:
                                diff_funcs.append((k, '%s VS %s : *updating*.' % (reffunc, newfunc)))
                                # *Updating*
                                suggested_transform[k] = ref_suggested_transform[k]
                if onlynew:
                    print('Transforms only in new:', ', '.join(onlynew))
                if onlyref:
                    print('Transforms only in ref:', ', '.join(onlyref))
                if diff_funcs:
                    print('Different transforms:\n', '\n'.join('%35s\t%s' % t for t in diff_funcs))
                if diff_args:
                    print('Different transform args:\n', '\n'.join('%35s\t%s' % t for t in diff_args))

            for ft in list(suggested_transform.keys()):
                if ft not in alls.columns:
                    logger.warning('Hardcoded feature %s not available: delete.', ft)
                    suggested_transform.pop(ft)

            self._suggested_transform = suggested_transform
            return suggested_transform


    def do(self):
        """Run self.do_rates() first."""
        cs_rates = self.cs_rates

        data = self.data
        same_alls = self.same_alls
        dataset_params = self.dataset_params
        responses = self.responses
        features = self.features
        ref_suggested_transform = self.ref_suggested_transform
        impose_transform = self.impose_transform
        must_drop_features = self.must_drop_features

        ages_controled = data.ages_controled
        mean_errors = data.mean_errors

        print('\n# Merge features')
        self.alls = alls = pd.concat((mean_errors,
                                      data.ns[dataset_params],
                                      cs_rates,
                                      same_alls),
                                     axis=1, join='inner', sort=False,
                                     verify_integrity=True)

        y = responses[0]
        print('Variable Y :', y)
        print('%d observation × %d features' % alls.shape)
        if alls.shape[0] == 0:
            logger.error('NULL inner join from:\nmean_errors\n------\n'
                         + str(mean_errors.iloc[:5, :5]) + '\n'
                         '\ndata.ns[dataset_params]\n------\n'
                         + str(data.ns[dataset_params].iloc[:5, :5]) + '\n'
                         '\ncs_rates\n------\n'
                         + str(cs_rates.iloc[:5, :5]) + '\n'
                         '\nsame_alls\n------\n'
                         + str(same_alls.iloc[:5, :5]) + '\n'
                         )
            sys.exit(1)
        na_amount = alls.isna().sum(axis=0).sort_values(ascending=False)
        print('Amount of NA:', na_amount.head(10), sep='\n')
        print('Amount of Inf:', np.isinf(alls.select_dtypes(np.number)).sum(axis=0).sort_values(ascending=False).head(10), sep='\n')
        
        too_many_na = na_amount.index[(na_amount >= 0.9*alls.shape[0])].tolist()
        many_na = na_amount.index[(na_amount >= 0.5*alls.shape[0])].tolist()
        if too_many_na:
            logger.warning('**DROP** these columns with >90%% of NA: %s',
                           ' '.join(too_many_na))
            alls.drop(columns=too_many_na, inplace=True)
            self.features = features = [ft for ft in features if ft not in too_many_na]
            many_na = [ft for ft in many_na if ft not in too_many_na]
        if many_na:
            logger.warning('Columns with >50%% of NA (kept as is): %s',
                           ' '.join(many_na))

    #def do_transforms(self):
        # Calling the cached property: depends on self.alls and self.features.
        suggested_transform = self.suggested_transform

        # Remove constant features
        self.features = features = [ft for ft in features if ft in suggested_transform]
        # All binary variables should **NOT** be z-scored!
        bin_features = [ft for ft in features if suggested_transform[ft].__name__ == 'binarize']

        logger.warning('Replace Inf,-Inf by NaN')
        alls_transformed = alls.replace([-np.Inf, np.Inf], np.NaN).transform(suggested_transform)

        print('transformed -> Any NA:', alls_transformed.isna().sum(axis=0).any())
        print('transformed -> *All* NA:', alls_transformed.columns[alls_transformed.isna().all()])

        print('transformed shape:', alls_transformed.shape)

        self.a_t = a_t = alls_transformed

    #def do_norm(self):
        logger.debug('Will Z-score: %s', [ft for ft in responses+features if ft not in bin_features])
            #{ft: zscore for ft in responses+features if ft not in bin_features})
        logger.debug('Will join binary features: %s', bin_features)
        a_n = a_t.transform({ft: zscore for ft in responses+features
                             if ft not in bin_features})\
                 .join(alls_transformed[bin_features])
        print('\nResponse transforms and normalisations:\n'
                + '\n'.join('%s -> %s : Mean=%.4f Std=%.4f'
                            %(r, suggested_transform[r].__name__,
                              a_t[r].mean(), a_t[r].std())
                          for r in responses))

        print(a_n.shape)
        print('normed -> Any NA:', a_n.columns.values[a_n.isna().any(axis=0)])
        print('Check a_n.head():')
        display_html(a_n.head())
        self.na_rows = na_rows = a_n.isna().any(axis=1)
        if na_rows.any():
            print('Drop %d NA rows' % na_rows.sum())
            a_n.dropna(inplace=True)

        a_nn = a_n.select_dtypes(np.number)  # Numeric columns only.
        inf_rows = np.isinf(a_nn).any(axis=1)
        if inf_rows.any():
            print('Drop %d Inf rows' % inf_rows.sum())
            print('Data with Inf: column %s:\n%s' % (
                        a_nn.columns[np.isinf(a_nn).any()],
                        alls[~na_rows][inf_rows].head(10)))
            a_n = a_n[~inf_rows].copy(deep=False)
        
        constant_vars = a_n.columns[(a_n.agg(np.ptp) == 0)].tolist()
        if constant_vars:
            logger.warning('Introduced constant features (**DROPPING NOW**): %s', 
                           ' '.join(constant_vars))
            a_n.drop(columns=constant_vars, inplace=True)
            self.features = features = [ft for ft in features if ft not in constant_vars]
        self.a_n = a_n


    #def do_fitall(self):
        print('\n### Fit of all features (lasso)')

        #fitlasso, fit, pslopes, displayed =
        *_, pslopes0, preslopes0 = lassoselect_and_refit(a_n, y, features)
        self.displayed.extend((pslopes0, preslopes0))

        #sb.violinplot('null_dS_before', 'abs_age_dev', data=a_n, cut=0);
        #scatter_density('r2t_dS_mean', 'abs_brlen_dev', data=a_n, alpha=0.5)
        #scatter_density('ingroup_glob_len', 'abs_age_dev', data=a_n, alpha=0.5)
        #scatter_density('dS_rate_std', 'abs_age_dev', data=a_n, alpha=0.5)
        #scatter_density(alls.null_dist_before, alls.null_dS_before)

    #def do_pca(self):
        print('\n### Dimension reduction of features\n#### PCA')

        ft_pca, ft_pca_outputs = detailed_pca(a_n, features, abs_cov=True,
                                              make_corr=True, heat_dendro=True,
                                              dendro_pad=0.2)
        self.ft_pca = ft_pca
        self.displayed.extend(ft_pca_outputs)

        print('\n#### Factor analysis')

        # To take into account continuous and categorical variables

        FA, fa_outputs = detailed_pca(a_n, features, FA=True)
        self.FA = FA
        self.displayed.extend(fa_outputs)

        # Outlier study (can be separated by a line in PC1-PC2 space)
        # (example from Catarrhini exact m1w04_fsa)
        #scatter_density(transformed0[:,0], transformed0[:,1], alpha=0.4)
        #ax = plt.gca()
        #ax.plot([-3, 3], [0.25, -0.5])

        #xsep = [-3, 3]; ysep = [0.25, -0.5]
        #eq_sep = lambda x: -0.125 - 0.75/6*x
        ##x = np.array([-2, -1, 0, 1, 2])
        ##plt.plot(x, eq_sep(x), 'r.')

        #transformed0_out = transformed0[:,1] > eq_sep(transformed0[:,0])
        #plt.plot(transformed0[transformed0_out, 0], transformed0[transformed0_out, 1], 'r.')

        print('\n### Feature decorrelation')
        #must_drop_features
        try:
            a_n_inde = a_n.drop(must_drop_features, axis=1)
        except KeyError as err:
            logger.warning("Not in `self.a_n`:" + err.args[0])
            a_n_inde = a_n.drop(must_drop_features, axis=1, errors='ignore')
        print('%d Independent features (%d rows)' % a_n_inde.shape[::-1])

        self.missing_todecorr = [pair for pair in
                         (self.to_renormlogdecorr[0] + list(self.to_renormlogdecorr[1].values())
                          + self.to_renormdecorr[0] + list(self.to_renormdecorr[1].values()))
                         if pair[0] not in a_t or pair[1] not in a_t]
        if self.missing_todecorr:
            logger.warning('Unexpected missing pairs to renorm(log)decorr: %s',
                           self.missing_todecorr)
        to_renormlogdecorr = ([pair for pair in self.to_renormlogdecorr[0]
                               if pair[0] in a_t and pair[1] in a_t],
                              {key:pair for key,pair in self.to_renormlogdecorr[1].items()
                               if pair[0] in a_t and pair[1] in a_t})
        to_renormdecorr = ([pair for pair in self.to_renormdecorr[0]
                             if pair[0] in a_t and pair[1] in a_t],
                            {key:pair for key,pair in self.to_renormdecorr[1].items()
                             if pair[0] in a_t and pair[1] in a_t})
        # Check if the proposed decorrelation is in line with the transforms.
        # (subtract for logs, divide for raw/sqrt)
        def iter_renormlogdecorr():
            yield from enumerate(to_renormlogdecorr[0])
            yield from to_renormlogdecorr[1].items()
        for k, (var, corrvar) in iter_renormlogdecorr():
            msg = []
            var_transform = suggested_transform[var].__name__
            corrvar_transform = suggested_transform[corrvar].__name__
            if 'log' not in var_transform and 'binary' not in var_transform:
                msg.append('suggested_transform[%r] = %s' % (var, var_transform))
            if 'log' not in corrvar_transform and 'binary' not in corrvar_transform:
                msg.append('/ suggested_transform[%r] = %s' % (corrvar, corrvar_transform))
            if msg:
                logger.warning(' '.join(msg) + ': automatic switch from "logdecorr" to "decorr"')
                if isinstance(k, int):
                    to_renormlogdecorr[0].remove((var, corrvar))
                    to_renormdecorr[0].append((var, corrvar))
                else:
                    to_renormlogdecorr[1].pop(k)
                    to_renormdecorr[1][k] = (var, corrvar)

        # TODO: verbose change in correlation coef
        if a_n_inde.shape[0] != a_t.shape[0]:
            logger.error('a_n_inde and a_t have incompatible shapes: %s %s (NA rows: %d/%d)',
                         a_n_inde.shape, a_t.shape, na_rows.sum(), na_rows.shape[0])

        a_n_inde = renorm_logdecorrelate(a_n_inde, a_t[~na_rows],
                                         *to_renormlogdecorr[0],
                                         **to_renormlogdecorr[1])
        a_n_inde = renorm_decorrelate(a_n_inde, a_t[~na_rows],
                                      *to_renormdecorr[0],
                                      **to_renormdecorr[1])

        self.missing_zeros_todecorr = [pair for pair in self.to_subtract#[0]
                                  if pair[0] not in a_t or pair[1] not in a_t]
        if self.missing_zeros_todecorr:
            logger.warning('Unexpected missing pairs of zero measures to decorr: %s',
                           self.missing_zeros_todecorr)
        zeros_to_decorrelate = [pair for pair in self.to_subtract
                                if pair[0] in a_t and pair[1] in a_t]
        
        a_n_inde = logdecorrelate(a_n_inde, a_t[~na_rows], *zeros_to_decorrelate)

        # special_decorr
        a_n_inde['CpG_odds'] = zscore(log(alls.ingroup_mean_CpG / (alls.ingroup_mean_GC**2)))[~na_rows]
        print('%d independent features (%d rows)' % a_n_inde.shape[::-1])

        inde_features = [ft for ft in features if ft in a_n_inde]
        #print('inde_features', len(inde_features))
        inde_features += [colname for colname in a_n_inde.columns.difference(a_n.columns)]
        print('inde_features', len(inde_features))

        new_inde_features = set(inde_features) - set(features)
        # Description of the transformation for the user.
        self.decorr = {'R'+var1: '- '+var2 for var1, var2 
                        in to_renormlogdecorr[0] + zeros_to_decorrelate}
        self.decorr.update(CpG_odds='ingroup_mean_CpG / ingroup_mean_GC^2',
                **{'R'+var1: '/ '+var2 for var1, var2 in to_renormdecorr[0]},
                **{key: '%s / %s' % (v1, v2) for key, (v1,v2)
                              in to_renormlogdecorr[1].items()},
                           **{key: '%s / %s' % (v1, v2) for key, (v1,v2)
                              in to_renormdecorr[1].items()})

        self.decorr_suggested_transform = {'R'+var1: self.suggested_transform[var1]
                for var1,_ in to_renormlogdecorr[0] + to_renormdecorr[0] + zeros_to_decorrelate}
        self.decorr_suggested_transform.update(**{key: self.suggested_transform[var1]
            for key, (var1,_) in to_renormlogdecorr[1].items()},
            **{key: self.suggested_transform[var1]
            for key, (var1,_) in to_renormdecorr[1].items()})

        self.a_n_inde = a_n_inde
        self.inde_features = inde_features

        #ft_pca_inde, ft_pca_inde_outputs = detailed_pca(a_n_inde, inde_features)
        self.ft_pca_inde = ft_pca_inde = PCA(n_components=15)
        # The above line is necessary to get all attributes (e.g. covariance)
        ft_pca_inde.fit_transform(a_n_inde[inde_features]) # -> transformed data

        heatmap_cov(np.abs(ft_pca_inde.get_covariance()), inde_features,
                    make_corr=True).suptitle('Inde features absolute correlation (PCA)')
        self.displayed.append(plt.gcf())
        plt.show()

        #FA_inde, FA_inde_outputs = detailed_pca(a_n_inde, inde_features, FA=True)
        self.FA_inde = FA_inde = FactorAnalysis(n_components=15)
        self.transformed_inde = transformed_inde = FA_inde.fit_transform(a_n_inde[inde_features])

        heatmap_cov(np.abs(FA_inde.get_covariance()), inde_features, make_corr=True)
        self.displayed.append(plt.gcf())
        self.displayed[-1].suptitle('Inde features absolute correlation (FA)')
        plt.show()
        fig, (ax1, ax2) = plt.subplots(1, 2, sharey=True)
        scatter_density(transformed_inde[:,1], transformed_inde[:,0], alpha=0.4,
                        ax=ax1)
        ax1.set_ylabel('PC1')
        ax1.set_xlabel('PC2')
        scatter_density(transformed_inde[:,2], transformed_inde[:,0], alpha=0.4,
                        ax=ax2)
        ax2.set_xlabel('PC3')
        plt.show()
        self.displayed.append(fig)


    #def do_fit(self):
        print('\n### Fit of less colinear features')

        self.fitlasso, self.fit, pslopes, preslopes = lassoselect_and_refit(
                                a_n_inde, y, inde_features, atol=1e-3,
                                method='elastic_net', L1_wt=1, cov_type='HC1') #**psummary_kw
        self.slopes = pslopes.data
        self.reslopes = preslopes.data
        self.displayed.extend((pslopes, preslopes))

    #def do_randomforest(self):
        print('\n#### Random Forest Regression')
        RFcrossval_r2 = randomforest_regression(a_n_inde[inde_features], a_n_inde[y])

        reslopes2 = self.do_bestfit()
        self.do_worsttrees()

        return reslopes2


    @property
    def suggest_multicolin(self):
        try:
            return self._suggest_multicolin
        except AttributeError:
            print('\n#### Finding the most colinear features')
            #multicol_test(a_n[features]), multicol_test(a_n_inde[inde_features])
            protected_features = [ft for ft in self.protected_features
                                    if ft in self.inde_features]
            self._suggest_multicolin, outputs = loop_drop_eval(self.inde_features,
                                          lambda x: multicol_test(self.a_n_inde[x]),
                                          stop_criterion=20,
                                          protected=protected_features)
            self.displayed.extend(outputs)
            return self._suggest_multicolin


    def do_bestfit(self):

    #def do_dropcolinear(self):
        a_n_inde = self.a_n_inde
        inde_features = self.inde_features
        y = self.responses[0]

        print(self.suggest_multicolin)

        print('\n#### Refit without most colinear features')
        
        print('Drop trees with bad properties: null_{dS,dN,t,dist}_{after,before}, prop_splitseq')
        bad_props = {}
        for badp, badval in self.must_drop_data.items():
            if badp not in a_n_inde.columns:
                continue
            vcounts = a_n_inde[badp].value_counts()
            if vcounts.shape[0] > 2:
                logger.warning('%r not binary.', badp)
            if (a_n_inde[badp].isin(badval).sum()) > 0.05 * vcounts.sum():
                logger.warning('Discarding %s.isin(%s) trees will remove >5%% of the subtrees',
                               badp, badval)
            print(vcounts, '\n')
            bad_props[badp] = badval

        #a_n_inde2 = a_n_inde.query(' & '.join('(%s==0)' % p for p in bad_props))\
        #                .drop(bad_props, axis=1, errors='ignore')
        a_n_inde2 = self.a_n_inde2 = a_n_inde.loc[a_n_inde.isin(bad_props).any(axis=1)]
        print("New shape after removing bad data:", a_n_inde2.shape)
        self.inde_features2 = inde_features2 = [ft for ft in inde_features
                                                if ft not in list(bad_props)+self.suggest_multicolin]

        constant_vars = a_n_inde2[inde_features2].agg(np.ptp)
        print('Check for newly introduced constant variables:\n',
              constant_vars.sort_values())
        self.inde_features2 = inde_features2 = [ft for ft in inde_features2
                                                if constant_vars[ft]>0]
        #print('\n##### OLS fit (2)')
        #self.fit2 = fit2 = ols.fit()
        #fit2.summary()
        ##display_html(sm_pretty_slopes(fit2))
        #display_html(fit2.summary())
        
        print('\n##### LASSO fit (2)')
        param_info = pd.concat((
                        pd.Series({ft: func.__name__ for ft,func
                                   in self.suggested_transform.items()},
                                  name='transform'),
                        self.a_t[[ft for ft in inde_features2 if ft not in
                                  self.decorr]].agg(['mean', 'std'])\
                                .rename(lambda c: 'transformed '+c).T,
                        pd.Series({ft: self.decorr.get(ft, '')
                                    for ft in inde_features2}, name='decorr')
                        ),
                        axis=1, sort=False)
        fitlasso2, refitlasso2, slopes2_styled, reslopes2_styled = \
                lassoselect_and_refit(a_n_inde2, y, inde_features2, join=param_info)
        self.displayed.extend((slopes2_styled, reslopes2_styled))
        self.fitlasso2 = fitlasso2
        self.refitlasso2 = refitlasso2
        self.slopes2_styled = slopes2_styled
        self.reslopes2_styled = reslopes2_styled

        self.slopes2 = slopes2_styled.data
        reslopes2 = self.reslopes2 = self.reslopes2_styled.data
        self.selected_features2 = reslopes2.drop('const').index.tolist()
        self.F_pval2 = refitlasso2.f_pvalue
        self.lL2 = refitlasso2.llf
        #print('fstat = %g\nP(F > fstat) = %g\nlog-likelihood = %g' %(
        #        refitlasso2.fvalue, self.F_pval2, self.lL2))

        figsize = mpl.rcParams['figure.figsize']
        fig, axes = plt.subplots(3, figsize=(figsize[0], figsize[1]*2))
        # Plot Y~X with largest coef.
        ax = axes[0]
        scatter_density(a_n_inde2[y], a_n_inde2[self.slopes2.index[1]], alpha=0.4, ax=ax)
        x = np.array(ax.get_xlim())
        a, b = reslopes2.loc[['const', self.slopes2.index[1]], 'Simple regression coef']
        ax.plot(x, 0 + b*x, '--', label='Simple Regression line (a=%g b=%g)' % (0,b))
        a, b = reslopes2.loc[['const', self.slopes2.index[1]], 'Lasso coef']
        ax.plot(x, a + b*x, '--', label='Lasso regression line (a=%g b=%g)' % (a,b))
        ax.legend()
        ax.set_title('Response against the feature with the largest Lasso coef.')
        ax.set_ylabel(y)
        ax.set_xlabel(self.slopes2.index[1])

        # Residual plot
        ax = axes[1]
        scatter_density(refitlasso2.fittedvalues,
                        a_n_inde2[y] - refitlasso2.fittedvalues,
                        alpha=0.4, ax=ax)
        ax.set_title('Residuals plot')
        ax.set_ylabel('Residual error')
        ax.set_xlabel('Predicted response')

        smg.gofplots.qqplot(refitlasso2.resid, line='r', ax=axes[2])
        plt.show()
        self.displayed.append(fig)

        if hasattr(refitlasso2, 'cov_HC0'):
            fig = heatmap_cov(refitlasso2.cov_HC0, ['const']+self.selected_features2,
                        cmap='seismic', make_corr=True)
            fig.suptitle('cov_HC0')
            plt.show()
            self.displayed.append(fig)
        else:
            logger.warning("No attribute 'cov_HC0' in `refitlasso2`.")
        if hasattr(refitlasso2, 'cov_HC1'):
            fig = heatmap_cov(refitlasso2.cov_HC1, ['const']+self.selected_features2,
                        cmap='seismic', make_corr=True)
            fig.suptitle('cov_HC1')
            plt.show()  # Not needed because this is the last figure before return.
            self.displayed.append(fig)
        else:
            logger.warning("No attribute 'cov_HC1' in `refitlasso2`.")

        # See this page for all possible accessible attributes/methods:
        # https://www.statsmodels.org/stable/generated/statsmodels.regression.linear_model.OLSResults.html#statsmodels.regression.linear_model.OLSResults

        return reslopes2

    def do_worsttrees(self):
        print('\n### Investigate the worst trees')
        display_html(self.alls.sort_values(self.responses[0], ascending=False).head(50))

    def show(self):
        for out in self.outputs:
            if isinstance(out, (pd.Series, pd.DataFrame, pd.io.formats.style.Styler)):
                display_html(out)
            elif isinstance(out, mpl.figure.Figure):
                display(out)
            else:
                print(out)


from siphon import dependency_func, dependency, auto_internalmethod
from functools import wraps


class analysis(object):
    """Analysis pipeline and environment"""
    
    @property
    def measures(self):
        return self._measures
    
    @property
    def dist_measures(self):
        return self._dist_measures

    @property
    def age_measures(self):
        return self._age_measures

    @property
    def rate_measures(self):
        return self._rate_measures

    @measures.setter
    def set_measures(self, value):
        self._measures = value
        self._dist_measures = ['branch_%s' %m for m in self._measures]
        self._rate_measures = ['rate_%s' %m for m in self._measures]
        self._age_measures =  ['age_%s' %m for m in self._measures]


    def __init__(self, age_file, stats_template,
                 stattypes=('al', 'tree', 'codeml'),
                 measures = ['dist', 'dS', 'dN', 't'],
                 control_condition='really_robust && aberrant_dists == 0',
                 phyltreefile=None,
                 control_ages_CI=None,
                 common_info=None,
                 al_params=None,
                 tree_params=None,
                 cl_params_restricted=None):
        # Add a settings attribute that stores the keys
        self.settings = set()
        for k,v in locals():
            if k != 'self':
                setattr(self, k, v)
                self.settings.add(k)

    @dependency
    def phyltree(self):
        return myPhylTree.PhylogeneticTree(phyltreefile)
    
    @dependency
    def ts(self):
        return load_stats_tree(stats_template.format(stattype=stattypes[1]))

    @dependency
    def aS(self):
        return load_stats_al(stats_template.format(stattype=stattypes[0]))

    @dependency
    def cs(self):
        return load_stats_codeml(stats_template.format(stattype=stattypes[2]))

    @dependency
    def clS(self):
        return load_any_stats(stats_template.format(stattype=stattypes[3]))

    @dependency_func
    @wraps(load_prepare_ages)
    def load_prepare_ages(self):
        return load_prepare_ages(self.ages_file, self.ts)

    @dependency
    def ages_treestats(self):
        return self.load_prepare_ages()[0]

    @dependency
    def ns(self):
        return self.load_prepare_ages()[1]

    @auto_internalmethod
    @wraps(add_control_dates_lengths)
    def get_control_dates_length(self, *args, **kwargs): 
        return add_control_dates_lengths(*args, **kwargs)

    @dependency_func
    def getcheck_control_dates_lengths(self, root):
        a_ctrled, ctrl_ages, ctrl_brlen = self.get_control_dates_length()
        check_control_dates_lengths(ctrl_brlen, self.phyltree, root)
        return a_ctrled, ctrl_ages, ctrl_brlen

    @property
    def ages_controled(self):
        return self.get_control_dates_length[0]

    @property
    def control_ages(self):
        return self.get_control_dates_length[1]
    
    @property
    def control_brlen(self):
        return self.get_control_dates_length[2]

    @dependency_func
    def mean_errors(self, control='median'):
        return compute_dating_errors(self.ages_controled, control=control)

    @dependency_func
    def mean_speciation_errors(self, control='median'):
        pass

    def compare_error_with(other_analysis):
        pass

    compute_branchrate_std = dependency(auto_internalmethod(compute_branchrate_std))

    @dependency
    def csrates(self):
        return self.compute_branchrate_std

    @dependency
    def alls(self):
        return pd.concat((self.mean_errors,
                          self.aS[self.common_info + self.al_params],
                          self.ts[self.tree_params],
                          self.cs[self.cl_params_restricted],
                          self.cs_rates),
                         axis=1, join='inner')


if __name__ == '__main__':

    # # Load data


    # ## Load ages

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

    new_feature_names = {ft: (func.__name__ + ("" if "%s" in func.__name__ else "(%s)")) % ft for ft, func in totransform.items() if func.__name__ != "notransform"}

    alls_transformed = alls.transform(totransform)
    assert not alls_transformed.isna().sum(axis=0).any()

    # Standardize/center variables

    forgotten = set(alls.columns) - set(totransform)
    assert not forgotten, forgotten

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

    #mpl.rcParams['figure.figsize'] = (16, 8) # width, height

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

    a = alls
    plt.plot(a.NsynSites / (a.ls*3), a.NnonsynSites / (a.ls*3), '.', alpha=0.7); # Alright, Ok.

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
    sm_pretty_summary(olsfit0)

    # #### robust trees

    alls_inde_normed.head()

    # Add intercept
    data = alls_inde_normed[(alls_inde_normed.Ndup == 0) & (alls_inde_normed.Nspe == 7)]
    olsfitr = sm.OLS(data.abs_error, sm.add_constant(data[inde_features])).fit()
    sm_pretty_summary(olsfitr)


    # ### Same with the formula syntax

    formula = 'abs_error ~ ' + ' + '.join(inde_features)
    print(formula)
    ols = smf.ols(formula, data=alls_inde_normed)
    results = ols.fit()
    sm_pretty_summary(results)


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

    sm_pretty_summary(olsfit_sq)

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

    sm_pretty_summary(olsfit_nodup_sq)
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

    ages_PL1 = pd.read_csv("/users/ldog/glouvel/ws2/DUPLI_data85/alignments_analysis/ages/Simiiformes_m1w04_ages.subtreesCleanO2-um2-withSG-PL1.tsv",
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




