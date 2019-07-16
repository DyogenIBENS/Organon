#!/usr/bin/env python3
# -*- coding: utf-8 -*-


"""Various statistical and graph routines.

Warning: some of them can only be run in a notebook (because using display_html
or Pandas dataframe styling.)"""

from io import StringIO
from functools import partial
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from IPython.display import display_html
import warnings

import scipy.cluster.hierarchy as hclust

from sklearn.decomposition import PCA
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import cross_val_score

from datasci.graphs import scatter_density, \
                           plot_cov, \
                           plot_features_radar
from datasci.stats import normal_fit, cov2cor  #, r_squared, adj_r_squared
from datasci.dataframe_recipees import centered_background_gradient, magnify


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

def make_lognegtransform_inc(inc=0):
    logneginc = lambda x: np.log10(-x + inc)
    logneginc.__name__ = "log10(%g-%%s)" % inc
    return logneginc

def make_logpostransform_inc(inc=0):
    loginc = lambda x: np.log10(x + x.min() + inc)
    loginc.__name__ = "log10(%g+min+%%s)" % (inc)
    return loginc


def make_best_logtransform(var, quantile=0.05):
    inc = 0
    if (var < 0).any():
        if (var > 0).any():
            # Both negative and positive values
            current_logtransform = make_logpostransform_inc
        else:
            current_logtransform = make_lognegtransform_inc
            if (var==0).any():
                inc = -var[var!=0].quantile(1-quantile)
    else:
        current_logtransform = make_logtransform_inc
        if (var==0).any():
            inc = var[var!=0].quantile(quantile)

    return current_logtransform(inc)


def discretizer(bins=[0]):
    """"""
    def discretize(x):
        return np.digitize(x, bins, right=True)
    return discretize


def binarize(x):
    return np.digitize(x, [0], right=True)


def zscore(x):
    return (x - x.mean()) / x.std() # if x.std() else 1) )


def decorrelator(decorrfunc, data, src_data=None, *args, **kwargs):
    if src_data is None:
        src_data = data

    kwargs.update({'R%s' % v[0]: v for v in args})

    return data.assign(**{newvar: decorrfunc(src_data[var],src_data[corrvar])
                          for newvar, (var, corrvar) in kwargs.items()})\
               .drop([var for var,_ in kwargs.values()],
                     axis=1, errors='ignore')

decorrelate = partial(decorrelator, np.divide)
logdecorrelate = partial(decorrelator, np.subtract)

renorm_decorrelate = partial(decorrelator,
                             lambda v,cv: zscore(np.divide(v,cv)))
renorm_logdecorrelate = partial(decorrelator,
                                lambda v,cv: zscore(np.subtract(v,cv)))


def test_transforms(alls, variables, figsize=(14, 5)):
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
        if var.min() == var.max():
            print("Variable %r is constant. Skipping." % ft)
            # Why does it appear in `suggested_transform`?
            continue
        
        fig, axes = plt.subplots(1, 3, figsize=figsize)
        #axes[i, 0].set_ylabel(ft)
        axes[0].set_ylabel(ft, fontsize='large')

        axes[0].hist(var, bins=nbins, density=True)
        axes[0].plot(*normal_fit(var), '-')
        axes[0].set_title("original")
        _, xmax0 = axes[0].get_xlim()
        _, ymax0 = axes[0].get_ylim()
        
        varskew = var.skew()
        text = "Skew: %g\nKurtosis: %g\n" % (varskew, var.kurt())
        transform_skews[notransform] = varskew

        current_logtransform = None
        if (var < 0).any():
            if (var > 0).any():
                print("Variable %r has negative and positive values. Shifting to positive." % ft)
                text += "Negative and positive values. Shifting to positive.\n"
                var -= var.min()
                current_logtransform = make_logpostransform_inc()
            else:
                print("Variable %r converted to positive values" % ft)
                text += "Converted to positive values.\n"
                var = -var
                current_logtransform = logneg  # What if null values?
        else:
            current_logtransform = log
        
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
        
        if not n_infinite_vals:
            transform_skews[current_logtransform] = logskew
        else:
            make_logtransform = make_logtransform_inc if current_logtransform is log \
                                else make_logpostransform_inc
            
            suggested_increment = logtransformed_var.quantile(0.05)
            print("%s: Nb of not finite values: %d. Suggested increment: 10^(%g) = %g"\
                  % (ft, n_infinite_vals, suggested_increment, 10**suggested_increment))
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
            
            log_inc_varskew = logtransformed_inc_var.skew()
            logtext_inc = ("Skew     (+inc): %g\nKurtosis (+inc): %g\n"
                           % (log_inc_varskew, logtransformed_inc_var.kurt()))
            transform_skews[make_logtransform(10**suggested_increment)] = log_inc_varskew

            log_inc1_varskew = logtransformed_inc1_var.skew()
            logtext_inc1 = ("Skew     (+1): %g\nKurtosis (+1): %g\n"
                            % (log_inc1_varskew, logtransformed_inc1_var.kurt()))
            transform_skews[make_logtransform(1)] = log_inc1_varskew

        xmin1, xmax1 = axes[1].get_xlim()
        _, ymax1 = axes[1].get_ylim()

        axes[1].set_title("log10 transformed")
        
        sqrttransformed_var = np.sqrt(var)
        sqrttransformed_varskew = sqrttransformed_var.skew()
        sqrttext = "Skew: %g\nKurtosis: %g\n" % (sqrttransformed_varskew,
                                                 sqrttransformed_var.kurt())
        transform_skews[sqrt] = sqrttransformed_varskew

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


def randomforest_regression(X, y, cv=10, n_estimators=100, max_depth=5, **kwargs):
    RF = RandomForestRegressor(n_estimators=n_estimators, max_depth=max_depth)
    RFfit = RF.fit(X, y)
    ft_importance = pd.Series(RFfit.feature_importances_,
                              index=X.columns).sort_values(ascending=False)
    print("R² =", RFfit.score(X, y))
    print(ft_importance)

    # a list of R² for each fit without a subsample
    RFcrossval = cross_val_score(RF, X, y, cv=cv, n_jobs=-1)
    print(RFcrossval)

    return RFfit


def print_pca_loadings(ft_pca, features):
    print("### PC loadings")
    print("%-17s:\t%10s\t%10s" % ('Feature', 'coef PC1', 'coef PC2'))
    for ft, coef1, coef2 in sorted(zip(features, ft_pca.components_[0,:],
                                       ft_pca.components_[1,:]),
                                   key=lambda x: (abs(x[1]), abs(x[2])),
                                   reverse=True):
        print("%-17s:\t%10.6f\t%10.6f" % (ft, coef1, coef2))


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
    ax.set_ylabel("Cumulative ratio of variance explained")
    plt.show()

    # Coefficients of the linear combination of each parameter in the resulting components
    print("Components dimensions:", ft_pca.components_.shape)

    # ### PC loadings

    # weights of each feature in the PCs

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
    fig.suptitle("Features in Principal Component space")
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
    #if var not in data.columns and re.search(r'[+-/*]', var):
    #    var = data[var]
    #    correlated_var = data[correlated_var]

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

def lm_summary(lm, features, response, data):
    """Display the summary statistics of the sklearn multiple linear regression."""
    print("R^2       =", lm.score(data[features], data[response]))
    print("Intercept =", lm.intercept_)
    print("\nSlopes\n------")

    features_by_coef = sorted(zip(features, lm.coef_), key=lambda x: np.abs(x[1]), reverse=True)

    for ft, coef in features_by_coef:
        print("%-17s: %10.6f" % (ft, coef))


def sm_ols_summary(olsfit, renames=None):
    """Nicer look for a StatsModels OLS fit (sorted features by slope)."""
    summary = olsfit.summary()
    if summary is not None:
        r_coefs = pd.read_csv(StringIO(olsfit.summary().tables[1].as_csv()),
                              sep=r'\s*,\s*', index_col=0, engine='python')
    else:
        # For example a Lasso fit (OLS().fit_regularized())
        r_coefs = olsfit.params.to_frame('coef')

    # Sort
    r_coefs['abs_coef'] = r_coefs.coef.abs()
    r_coefs.sort_values("abs_coef", ascending=False, inplace=True)
    r_coefs.drop("abs_coef", axis=1, inplace=True)

    renames = {} if renames is None else renames

    param_names = [renames.get(n, n) for n in olsfit.model.exog_names
                   if n not in ('Intercept', 'const')]
    r_coefs_styled = r_coefs.rename(renames).style.bar(
                        subset=pd.IndexSlice[param_names, "coef"],
                        axis=0,
                        align="zero")
    return r_coefs_styled


def drop_eval(features, func):
    """Example: func=lambda x: multicol_test(df[x])"""
    out = pd.Series(0, index=features, dtype=float)
    for i, ft in enumerate(features):
        out[features[i]] = func(features[:i] + features[(i+1):])

    return out


def display_drop_eval(features, func):
    return drop_eval(features, func).to_frame.style.bar()


def loop_drop_eval(features, func, criterion='min', nloops=None):
    if nloops is None:
        nloops = len(features)

    dropped_features = []
    next_features = [ft for ft in features]
    for k in range(nloops):
        dropped_k = drop_eval(next_features, func)
        display_html(dropped_k.to_frame().style.bar())
        
        if criterion == 'min':
            i = dropped_k.values.argmin()  # Numpy argmin -> get an integer index.
        elif criterion == 'max':
            i = dropped_k.values.argmax()  # Numpy argmin -> get an integer index.
        else:
            raise ValueError('Unknown criterion %r (min/max)' % criterion)
        
        print('%d. DROP %s' % (k, next_features[i]))
        dropped_features.append(next_features[i])
        next_features = next_features[:i] + next_features[(i+1):]

    return dropped_features

