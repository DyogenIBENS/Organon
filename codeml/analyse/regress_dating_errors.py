#!/usr/bin/python3
# coding: utf-8

# Source:
# jupyter nbconvert --to python \
#   ~/ws2/DUPLI_data85/alignments_analysis/subtrees_stats/subtrees_stats.ipynb

import sys
import warnings
from copy import copy, deepcopy
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
from dendro.bates import dfw_pairs_generalized, dfw_pairs
from datasci.graphs import scatter_density, \
                           plot_cov, \
                           heatmap_cov, \
                           plot_loadings, \
                           plot_features_radar, \
                           plottree
from datasci.stats import normal_fit, cov2cor
from datasci.dataframe_recipees import centered_background_gradient, magnify
from datasci.compare import pairwise_intersections, align_sorted

from dendro.any import myPhylTree as phyltree_methods, ete3 as ete3_methods
import ete3

from scipy import stats
import scipy.cluster.hierarchy as hclust
import scipy.spatial.distance as spdist
#stats.skew, stats.kurtosis

from LibsDyogen import myPhylTree

from sklearn.decomposition import PCA
from sklearn.linear_model import LinearRegression, Lasso
from sklearn.preprocessing import StandardScaler
import statsmodels.api as sm
#import statsmodels.formula.api as smf
import statsmodels.stats.api as sms

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
sh = logging.StreamHandler(sys.stdout)
sh.setFormatter(colorlogf)
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
logger.handlers = []
logger.addHandler(sh)
#logging.basicConfig(handlers=[sh])

mpl.style.use("softer")
pd.set_option("display.max_columns", 50)
pd.set_option("display.width", 115)
pd.set_option("display.max_colwidth", 50)
pd.set_option("display.show_dimensions", True)  # even on non truncated dataframes

# wide_screen_style
mpl.rcParams['figure.figsize'] = (14, 10) # width, height


# Convert "time used" into seconds.  # ~~> numbertools? timetools? converters?
def time2seconds(time_str):
    factors = [1, 60, 3600, 3600*24]
    s = 0
    for factor, n_units in zip(factors, reversed(time_str.split(':'))):
        s += factor * int(n_units)
    return s


def load_stats_al(alfile):
    return pd.read_table(alfile, index_col=0)

def load_stats_tree(treefile):
    ts = pd.read_table(treefile, index_col=0,
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
    cs = pd.read_table(codemlfile, index_col=0)
    cs['seconds'] = cs['time used'].apply(time2seconds)
    return cs

stat_loaders = {'al':     load_stats_al,
                'tree':   load_stats_tree,
                'codeml': load_stats_codeml}

stat_loaders['codemlI'] = stat_loaders['codeml']
stat_loaders['treeI'] = stat_loaders['tree']


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
        output_tables.append(stat_loaders[stattype](filename))

    return tuple(output_tables)


def check_load_subtree_stats(aS, ts, cs):
    logger.info("shapes: aS %s, ts %s, cs %s" % (aS.shape, ts.shape, cs.shape))
    logger.info("aS has dup: %s", aS.index.has_duplicates)
    logger.info("ts has dup: %s", ts.index.has_duplicates)
    logger.info("cs has dup: %s", cs.index.has_duplicates)
    common_subtrees = set(aS.index) & set(ts.index) & set(cs.index)
    logger.info("%d common subtrees" % len(common_subtrees))
    only_al = aS.index.difference(ts.index.union(cs.index))
    only_tr = ts.index.difference(aS.index.union(cs.index))
    only_co = cs.index.difference(aS.index.union(ts.index))
    l_al = len(only_al)
    l_tr = len(only_tr)
    l_co = len(only_co)
    logger.warning("%d only in al stats: %s" % (l_al, list(only_al)[:min(5, l_al)]))
    logger.warning("%d only in tree stats: %s" % (l_tr, list(only_tr)[:min(5, l_tr)]))
    logger.warning("%d only in codeml stats: %s" % (l_co, list(only_co)[:min(5, l_co)]))
    # Todo: pyupset plot



# ## Function to merge additional subgenetree information into `ages`
def merge_criterion_in_ages(criterion_serie, ages=None, ages_file=None,
                            criterion_name=None):
    """Merge a column into the *node ages* table: the common field is the *subgenetree*."""
    
    assert (ages is not None or ages_file) and not (ages_file and ages is not None), "At least `ages` (dataframe) or `ages_file` (filename) must be given."
    if ages is None:
        ages = pd.read_table(ages_file, sep='\t')
    
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


def add_robust_info(ages_p, ts, measures=['dist', 'dS', 'dN', 't']):
    # Compute Number of duplications/speciation per tree.
    sgg = subgenetree_groups = ages_p.groupby('subgenetree')
    ts['Ndup'] = sgg.type.agg(lambda v: sum(v == "dup"))
    print(ts.Ndup.describe())
    ts['Nspe'] = sgg.type.agg(lambda v: sum(v == "spe"))
    print(ts.Nspe.describe())
    agg_funcs = {'consecutive_null_%s' %m: 'sum' for m in measures}
    agg_funcs.update({'branch_%s' %m: lambda v: (v==0).mean()
                      for m in measures})
    ts.join(sgg.agg(agg_funcs)).rename(columns={'branch_%s' %m: 'freq_null_%s' %m
                                                for m in measures})
    # merge tree stats to select robust trees
    ages_treestats = pd.merge(ages_p.drop('_merge', axis=1),
                              ts[['Ndup', 'Nspe', 'really_robust',
                                  'aberrant_dists', 'rebuilt_topo']],
                              how='left',
                              left_on='subgenetree',
                              right_index=True,
                              indicator=True, validate='many_to_one')
    logger.info("Ndup %s; Nspe %s; ages_treestats %s",
                (~ts.Ndup.isna()).sum(), (~ts.Nspe.isna()).sum(),
                ages_treestats.shape)
    logger.info('merge types:\n%s',
                ages_treestats.groupby('_merge')['_merge'].count())
    return ages_treestats
    

def load_prepare_ages(ages_file, ts, measures=['dist', 'dS', 'dN', 't']):
    """Load ages dataframe, join with parent information, and compute the
    'robust' info"""
    ages = pd.read_table(ages_file, sep='\t', index_col=0)

    logger.info("Shape ages: %s; has dup: %s" % (ages.shape, ages.index.has_duplicates))
    n_nodes = ages.shape[0]
    ages_int = ages[ages.type != 'leaf']
    logger.info("Shape ages internal nodes: %s; has dup: %s",
                ages_int.shape,
                ages_int.index.has_duplicates)
    n_nodes_int = (ages.type != 'leaf').sum()

    # Fetch parent node info
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
    recognized_outgroups = (ages_orphans.is_outgroup == 1)
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
    # However don't drop the ingroup root
    # (doesn't have an is_outgroup_parent info)

    child_of_root = ages_p.root == ages_p.parent
    #recognized_outgroups = (ages_p.is_outgroup == 1 | ages_orphans.is_outgroup_parent == 1)
    to_remove = (~child_of_root & orphans
                 | (ages_p.is_outgroup == 1)
                 | ages_p.is_outgroup_parent == 1)
    # If some nodes to remove appear to be valid ingroup nodes, it's because their
    # ancestors could not be dated (absence of the node required for calibration,
    # e.g Simiiformes)
    ages_p = ages_p.loc[~to_remove].copy(deep=False)
    if ages_p.index.has_duplicates:
        debug_columns = ['parent', 'subgenetree', 'taxon', 'taxon_parent', 'median_taxon_age', 'branch_dS', 'age_dS']
        logger.error("Failed to remove index duplicates in 'ages_p':\n%s\n...",
                     ages_p.loc[ages_p.index.duplicated(), debug_columns].head(10))
    ages_spe2spe = ages_p[(ages_p.type.isin(('spe', 'leaf')))
                          & (ages_p.type_parent == 'spe')]
    logger.info("\nShape ages speciation to speciation branches (no dup): %s",
                ages_spe2spe.shape)
    for m in measures:
        ages_p['consecutive_null_%s' %m] = \
                ~ages_p[['branch_%s_parent' %m, 'branch_%s' %m]].any(axis=1)
    ages_treestats = add_robust_info(ages_p, ts)

    #ages_robust = ages_treestats[ages_treestats.really_robust & \
    #                             (ages_treestats.aberrant_dists == 0)]\
    #                        .drop(['really_robust', 'aberrant_dists'], axis=1)
    return ages_treestats


# for averaging by taking into account branch length: with Omega.
# NOTE: columns will be reordered following `var`.
def group_average(g, var, weight_var="median_brlen"):
    g = g[~g[weight_var].isna()]
    if not g.shape[0]:
        return pd.Series([np.NaN]*len(var))
    values = np.ma.array(g[var], mask=g[var].isna())
    return pd.Series(np.average(g[var], axis=0, weights=g[weight_var]))
    # TODO: append the count of removed rows (NA in weight_var)

def group_weighted_std(g, var, weight_var="median_brlen"):
    return pd.Series(weighted_std(g[var], axis=0, weights=g[weight_var]))

# for calculating rates: with dS, dN, t
def tree_dist_2_rate(g, dist_var, norm_var="median_brlen"):
    # in pandas, sum is done per columns.
    return g[dist_var].sum() / g[norm_var].sum()


def add_control_dates_lengths(ages, phyltree, timetree_ages_CI=None,
                              measures=['dist', 'dS', 'dN', 't'],
                              control_condition='really_robust & aberrant_dists == 0'):
    # Merge control dates
    ages_forcontrol = ages.query(control_condition)
    logger.info("%d nodes from robust trees", ages_forcontrol.shape[0])
    
    median_taxon_ages = ages_forcontrol[ages_forcontrol.type\
                                                       .isin(("spe", "leaf"))]\
                                       .groupby("taxon").age_dS.median()\
                                       .rename('median_taxon_age')

    timetree_ages = median_taxon_ages.index.to_series()\
                                     .apply(phyltree.ages.get)\
                                     .rename('timetree_age')

    control_ages = pd.concat((median_taxon_ages, timetree_ages,
                              timetree_ages_CI), axis=1, sort=False)

    #print(control_ages.sort_values('timetree_age', ascending=False))

    ages_controled = pd.merge(ages, control_ages,
                              left_on="taxon", right_index=True,
                              validate="many_to_one")

    # Merge control branch lengths
    invalid_taxon_parent = ages_controled.taxon_parent.isna()
    invalid_measures  = ages_controled[['age_%s' %m for m in measures]].isna().all(1)
    invalid_nodes = invalid_taxon_parent & ~invalid_measures
    # calibrated | branch_dS.isna()
    # Should be nodes whose parent node is the root.
    debug_columns = ['parent', 'subgenetree', 'taxon', 'taxon_parent', 'median_taxon_age', 'branch_dS', 'age_dS']
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

    #TODO: NaN for data not in ages_forcontrol?
    ages_controled['median_brlen'] = \
        ages_controled_spe2spe.taxon_parent.apply(control_ages.median_taxon_age.get) \
        - ages_controled_spe2spe.median_taxon_age
    #control_ages.reindex(ages_controled.taxon_parent)\
    #        .set_axis(ages_controled.index, inplace=False)
    # would be more 'Pandas-like'.

    ages_controled['timetree_brlen'] = \
        ages_controled_spe2spe.taxon_parent.apply(control_ages.timetree_age.get) \
        - ages_controled_spe2spe.timetree_age

    # Resulting branch lengths
    branch_info = ["taxon_parent", "taxon"]
    #control_brlen = ages_controled.loc[
    #                    ~ages_controled.duplicated(branch_info),
    control_brlen = ages_controled.query('type != "dup" & type_parent != "dup"')\
                    .groupby(branch_info)\
                    [["median_brlen", "median_taxon_age", "timetree_brlen",
                      "timetree_age"]]\
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
    return unexpected_branches, lost_branches


def compute_dating_errors(ages_controled, control='median'):
    age_var = control + '_taxon_age'  # Control value
    brlen_var = control + '_brlen'    # Control value

    ages_controled["abs_age_error"] = \
                (ages_controled.age_dS - ages_controled[age_var]).abs()
    ages_controled["signed_age_error"] = \
                (ages_controled.age_dS - ages_controled[age_var])
    ages_controled["abs_brlen_error"] = \
                (ages_controled.age_dS_parent - ages_controled.age_dS - \
                 ages_controled[brlen_var]).abs()
    ages_controled["signed_brlen_error"] = \
                (ages_controled.age_dS_parent - ages_controled.age_dS - \
                 ages_controled[brlen_var])

    mean_errors = ages_controled[
                     ['subgenetree', 'abs_age_error', 'signed_age_error',
                      'abs_brlen_error', 'signed_brlen_error']
                    ].groupby("subgenetree").mean()
    return mean_errors



def display_evolutionary_rates():
    Glires_ts = load_stats_tree('subtreesGoodQualO2_treestats-Glires.tsv')

    measures = ['dist', 'dS', 'dN', 't']
    dist_measures = ['branch_%s' % m for m in measures]
    rate_measures = ['%s_rate' % m for m in measures]

    agesGlires_treestats = load_prepare_ages('../ages/Glires_m1w04_ages.subtreesGoodQualO2-um2-ci.tsv', Glires_ts)

    agesGlires_controled, control_agesGlires, control_brlenGlires =\
        add_control_dates_lengths(agesGlires_treestats, phyltree)
    agesGlires_controled_robust = agesGlires_controled.query('really_robust & aberrant_dists == 0') 

    agesGlires_controled_robust['branchtaxa'] = agesGlires_controled_robust.taxon_parent + '--' + agesGlires_controled_robust.taxon

    agesGlires_controled_robust[rate_measures] = agesGlires_controled_robust[dist_measures]\
                                          .div(agesGlires_controled_robust.timetree_brlen, axis=0)

    agesGlires_controled_robust.groupby(['taxon_parent', 'taxon'])[rate_measures].agg(['median', 'mean', 'std'])\
        .style.background_gradient(cmap='PRGn', subset=[(r, stat) for r in rate_measures for stat in ('median', 'mean')])

    ordered_Glires_branches = ['%s--%s' % br for br in dfw_pairs(phyltree, queue=[(None, 'Glires')], closest_first=True)]

    ordered_Glires_branches_bylen = ['%s--%s' % v for v in control_brlenGlires.sort_values('timetree_brlen').index.values]

    sb.violinplot('branchtaxa', 'dS_rate', data=agesGlires_controled_robust, width=1, order=ordered_Glires_branches_bylen)
    ax = plt.gca()
    ax.set_ylim(-0.005, 0.03);
    ax.set_xticklabels([xt.get_text() for xt in ax.get_xticklabels()], rotation=45, va='top', ha='right');


def compute_branchrate_std(ages_controled, dist_measures,
                           branchtime='median_brlen', taxon_age=None):
    
    groupby_cols = ["subgenetree", "taxon_parent", "taxon",
                    branchtime] + dist_measures
    if taxon_age is not None:
        groupby_cols.append(taxon_age)  ## "median_taxon_age"
    #ages_controled["omega"] = ages_controled.branch_dN / ages_controled.branch_dS

    sgg = subgenetree_groups = ages_controled[groupby_cols].groupby('subgenetree')

    # ### Average (substitution) rates over the tree:
    #     sum of all branch values / sum of branch lengths

    # Sum aggregation + division broadcasted on columns
    # FIXME: actually wouldn't it be simpler dividing before groupby?
    cs_rates = sgg[dist_measures].sum().div(sgg[branchtime].sum(), axis=0)
    #cs_rates["omega"] = (sgg.branch_dN / sgg.branch_dS).apply() 
    rate_measures = [(m.replace('branch_', '') + '_rate') for m in dist_measures]
    cs_rates.columns = rate_measures

    # ### Weighted standard deviation of substitution rates among branches

    tmp = pd.merge(ages_controled[["subgenetree", branchtime] + dist_measures],
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
    rate_dev_dict = {d: (tmp[d] / tmp[branchtime] - tmp[r])**2
                     for d,r in zip(dist_measures, rate_measures)}
    
    rate_dev = pd.DataFrame(rate_dev_dict)\
                    .join(tmp[['subgenetree', branchtime]])

    cs_wstds = rate_dev.groupby("subgenetree").apply(
                (lambda x, var, weight_var:
                                    sqrt(group_average(x, var, weight_var))),
                dist_measures, branchtime)
                
    cs_wstds.columns = [(r + '_std') for r in cs_rates.columns]  # Or multiindex

    return pd.concat((cs_rates, cs_wstds), axis=1)


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
                                           split=True, order=None, **kwargs):
    Q_col = "Q_" + criterion_name
    if isin is None:
        # Look at extreme quantiles only
        isin = (annot_df[Q_col].min(), annot_df[Q_col].max())
        print(isin)
        
    ax = sb.violinplot(x="taxon", y="age_dS", hue=Q_col,
                       data=annot_df[(annot_df.type == "spe")
                                     & annot_df[Q_col].isin(isin)],
                       split=split,
                       order=order,
                       **kwargs)
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45)
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


# Functions to check if variables need transformation


def all_test_transforms(alls, variables, figsize=(14, 5)):
    #fig, axes = plt.subplots(len(variables),3, figsize=(22, 5*len(variables)))
    nbins = 50

    suggested_transform = {}
    for i, ft in enumerate(variables):
        
        var = alls[ft]
        
        transform_skews = {}
        # Plot original variable distribution 
        if var.dtype != float:
            print("Variable %r not continuous: %s" % (ft, var.dtype))
            if var.dtype != int:
                print("Variable %r does not seem to be numeric. Skipping" % ft)
                continue
        
        fig, axes = plt.subplots(1, 3, figsize=figsize)
        #axes[i, 0].set_ylabel(ft)
        axes[0].set_ylabel(ft)

        axes[0].hist(var, bins=nbins, density=True)
        axes[0].plot(*normal_fit(var), '-')
        axes[0].set_title("original")
        _, xmax0 = axes[0].get_xlim()
        _, ymax0 = axes[0].get_ylim()
        
        varskew = var.skew()
        text = "Skew: %g\nKurtosis: %g\n" % (var.skew(), var.kurt())
        transform_skews[notransform] = varskew

        if (var < 0).any():
            if (var > 0).any():
                print("Variable %r has negative and positive values. Shifting to positive." % ft)
                text += "Negative and positive values. Shifting to positive.\n"
                var -= var.min()
                #best_transform = make_logpostransform_inc()
            else:
                print("Variable %r converted to positive values" % ft)
                text += "Converted to positive values.\n"
                var = -var
                #best_transform = logneg
        
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
        #transform_skews[logtransform] = logskew

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
            #transform_skews[logtransform_inc(10**suggested_increment)] = logtransformed_inc_var.skew()
            
            logtext_inc1 = ("Skew     (+1): %g\nKurtosis (+1): %g\n"
                            % (logtransformed_inc1_var.skew(),
                               logtransformed_inc1_var.kurt()))
            #transform_skews[logtransform_inc(1)] = logtransformed_inc1_var.skew()

        xmin1, xmax1 = axes[1].get_xlim()
        _, ymax1 = axes[1].get_ylim()

        axes[1].set_title("log10 transformed")
        
        sqrttransformed_var = np.sqrt(var)
        sqrttext = "Skew: %g\nKurtosis: %g\n" % (sqrttransformed_var.skew(),
                                                 sqrttransformed_var.kurt())
        transform_skews[sqrt] = sqrttransformed_var.skew()

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
        suggested_transform[ft] = sorted(transform_skews.items(),
                                         key=lambda x: abs(x[1]))[0][0]

    return suggested_transform


# Variable transformation

def notransform(x):
    return x

#notransform.__name__ = "%s"
sqrt = np.sqrt
log = np.log10

def logneg(x):
    return np.log10(-x)
logneg.__name__ = "log10(-%s)"

def make_logtransform_inc(inc=0):
    loginc = lambda x: np.log10(x + inc)
    loginc.__name__ = "log10(%g+%%s)" % inc
    return loginc

def make_logpostransform_inc(inc=0):
    loginc = lambda x: np.log10(x + x.min() + inc)
    loginc.__name__ = "log10(%g+min+%%s)" % (inc)
    return loginc


# ~~> datasci.routines ?
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

    # Just so that the features with similar component values are grouped together.
    # TO give more weight in the clustering to the first PC, multiply the PC
    # by their eigenvalue.
    ft_dendro = hclust.dendrogram(
                            hclust.linkage(components*ft_pca.explained_variance_,
                                           'average'),
                            labels=features,
                            count_sort='descending',
                            no_plot=True)
    ordered_ft = ft_dendro['ivl']
    ft_order = ft_dendro['leaves']

    #.sort_values(["PC1", "PC2"])
    styled_components = components.loc[ordered_ft].style.\
            apply(centered_background_gradient, cmap="PRGn", extend=0.15).\
            set_caption("Principal Components loadings").\
            set_properties(**{'max-width': '80px', 'font-size': '1pt'}).\
            set_table_styles(magnify())
    print("Rendered_components:", type(styled_components), styled_components)
    display_html(styled_components)

    # ### Feature covariance

    ft_cov = ft_pca.get_covariance()
    print("Covariance dimensions:", ft_cov.shape)
    plot_cov(ft_cov[ft_order,:][:,ft_order], ordered_ft, cmap='seismic')
    plt.show()

    # Plot feature vectors in the PC space
    fig, (ax0, ax1) = plt.subplots(2, 1, figsize=(10, 20),
                                   subplot_kw={'projection': 'polar'})

    plot_features_radar(components, features, ax=ax0)
    plot_features_radar(components, features, PCs=["PC1", "PC3"], ax=ax1)
    fig.suptitle("Features in Principal Component space"); 
    plt.show()
    return ft_pca


# Functions for checking colinearity between variables

# ~~> datasci.routine ?
def check_decorrelate(var, correlated_var, data, logdata=None):
    _, axes = plt.subplots(2, 1 if logdata is None else 2)
    axes = axes.flat
    
    if logdata is None:
        ax0, ax2 = axes
    else:
        ax0, ax1, ax2, ax3 = axes
        
    scatter_density(correlated_var, var, alpha=0.5, s=9, data=data, ax=ax0)
    
    # Plot the uncorrelated variables
    uncor_var = data[var] / data[correlated_var]
    uncor_null = (data[var] == 0) & (data[correlated_var] == 0)
    uncor_var[uncor_null] = 0
    scatter_density(data[correlated_var], uncor_var, s=9, alpha=0.5, ax=ax2)
    
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

# ~~> datasci.stattools ?
def lm_summary(lm, features, response, data):
    """Display the summary statistics of the sklearn multiple linear regression."""
    print("R^2       =", lm.score(data[features], data[response]))
    print("Intercept =", lm.intercept_)
    print("\nSlopes\n------")

    features_by_coef = sorted(zip(features, lm.coef_), key=lambda x: np.abs(x[1]), reverse=True) 

    for ft, coef in features_by_coef:
        print("%-17s: %10.6f" % (ft, coef))


def sm_ols_summary(olsfit, renames=None):
    r_coefs = pd.read_csv(StringIO(olsfit.summary().tables[1].as_csv()),
                          sep='\s*,\s*', index_col=0, engine='python')

    # Sort
    r_coefs['abs_coef'] = r_coefs.coef.abs()
    r_coefs.sort_values("abs_coef", ascending=False, inplace=True)
    r_coefs.drop("abs_coef", axis=1, inplace=True)
    #r_coefs_styled = r_coefs.style.apply(centered_background_gradient, axis=0, subset="coef")
    renames = {} if renames is None else renames

    param_names = [renames.get(n, n) for n in olsfit.model.exog_names
                   if n not in ('Intercept', 'const')]
    r_coefs_styled = r_coefs.rename(renames).style.bar(
                        subset=pd.IndexSlice[param_names, "coef"],
                        axis=0,
                        align="zero")
    return r_coefs_styled



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
                 timetree_ages_CI=None,
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
    def ages_treestats(self):
        return load_prepare_ages(self.ages_file, self.ts)

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
        return compute_dating_error(self.ages_controled, control=control)

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

    #add_control_dates_lengths(ages_treestats, ages_robust, phyltree, timetree_ages_CI=None)

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
    #check_control_dates_lengths(ages_controled, phyltree, root)

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

    new_feature_names = {ft: (func.__name__ + ("" if "%s" in func.__name__ else "(%s)")) % ft for ft, func in totransform.items() if func.__name__ != "notransform"}

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

    #mpl.rcParams['figure.figsize'] = (16, 8) # width, height

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




