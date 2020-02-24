#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import sys
import os.path as op
from pathlib import Path
import numpy as np
import pandas as pd
import scipy.stats as st
from scipy.stats import gamma, \
                        entropy  # Kullback-Leibler divergence
import matplotlib as mpl
mpl.use('Agg')
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import statsmodels.api as sm

from datasci.graphs import scatter_density, brokenAxes
from datasci.routines import test_transforms
from datasci.dataframe_recipees import bounded_background_gradient
from datasci.savior import HtmlReport, css_dark_style, reroute_loggers, generate_slideshow

from duprates.optim_gamma import inverse_scale #, gamma_negloglikelihood
from codeml.analyse.regress_dating_errors import *


float_ = np.longdouble  # Global precision set here.


import logging
logger = logging.getLogger(__name__)

plt.style.use('softer')
plt.style.use('softdark')  # Also see the default 'dark_background'
#TODO: set background of figure to dark. set better color cycle.
import seaborn as sb
#sb.set_palette('deep')
#sb.set_palette('pastel')
#sb.set_palette('Set2')

mpl.rcParams['figure.figsize'] = (12,7)
myslideshow_js = Path.home().joinpath('mydvpt', 'atavistic-doc-tools', 'myslideshow.js') 
myslideshow_css = Path.home().joinpath('mydvpt', 'atavistic-doc-tools', 'myslideshow.css') 


if __name__ == '__main__':
    try:
        from UItools import colorlog
        colorlog.install(logger)
    except ImportError:
        logging.basicConfig()
    logger.setLevel(logging.INFO)


continuous_distribs = [d for dname, d in sorted(vars(st).items())
                       if isinstance(d, st.rv_continuous)]
continuous_positive_distribs = [d for d in continuous_distribs if d.a >= 0 and np.isinf(d.b)]
# There are 50 distributions implemented.


def fmt_params(params, prec=4):
    fmt='%%.%dg' % prec
    return '(' + ', '.join(fmt % x for x in params) + ')'


def fit_distribs(distribs, df, x=None, nbins=100, stream=None):
    if x is None:
        x = df.columns.tolist()
    ncols = len(x)  # controling x is NotImplemented
    descriptions = ('dup', 'dup (>0)', 'loss', 'loss (dup>0)', 'loss (dup=0)')
    datas = (df.duprate.values,
            df.duprate[~nodup].values,
            df.lossrate.values,
            df.lossrate[~nodup].values,
            df.lossrate[nodup].values)
    for distrib in distribs:
        print('###', distrib.name, file=stream)
        fig, ((ax0top, ax1top), (ax0bottom, ax1bottom)) = plt.subplots(2, ncols=ncols)
        # Y axis with a broken scale : ---//---
        fig.subplots_adjust(hspace=0.1)
        axes_dup = brokenAxes(ax0bottom, ax0top)
        (h_dup, b_dup, _), _ = axes_dup.hist(df.duprate, bins=nbins, density=True)
        axes_dup.dobreak(max(h_dup[1:]))

        axes_loss = brokenAxes(ax1bottom, ax1top)
        ((h_loss_nodup,h_loss), b_loss, _), _ = axes_loss.hist((df.lossrate[nodup],
                                                           df.lossrate[~nodup]),
                                                label=['No dup', 'dup > 0'],
                                                bins=nbins, density=True, stacked=True)
        ymax = h_loss.max()
        axes_loss.dobreak(max(h_loss[h_loss<ymax]))

        fits = [distrib.fit(data) for data in datas]

        # X coordinates of the middle of the bins:
        xdup = b_dup[:-1] + (b_dup[1] - b_dup[0])/2     #, dtype=float_
        xloss = b_loss[:-1] + (b_loss[1] - b_loss[0])/2 #, dtype=float_
        dup_density = distrib.pdf(xdup, *fits[0])
        dup_density_nonzero = distrib.pdf(xdup, *fits[1])
        axes_dup.plot(xdup, dup_density,
                      label=fmt_params(fits[0]))
        axes_dup.plot(xdup, dup_density_nonzero,
                      label='dup>0 %s' % fmt_params(fits[1]))
        ax0top.legend()
        ax0bottom.set_ylabel('% of gene trees')
        ax0top.set_title('Duplication rates')

        loss_density = distrib.pdf(xloss, *fits[2])
        loss_density_nodup = distrib.pdf(xloss, *fits[4])
        loss_density_dup = distrib.pdf(xloss, *fits[3])

        # Kullback-Leibler divergence from the theoretical distribution.
        kl = [entropy(h_dup, dup_density),
              entropy(h_dup[1:], dup_density_nonzero[1:]),
              entropy(h_loss, loss_density),
              entropy(h_loss - h_loss_nodup, loss_density_dup),
              entropy(h_loss_nodup, loss_density_nodup)]
        # negative log-likelihood:
        nlL = [distrib.nnlf(fit, data) for fit, data in zip(fits, datas)]
        aic = [2*len(fit) + 2*nll for fit, nll in zip(fits, nlL)]

        for info in zip(descriptions, kl, nlL, aic):
            print('%-14s:  KL = %8g;  -lL = %8g;  AIC = %g' % info, file=stream)

        axes_loss.plot(xloss, loss_density,
                       label=fmt_params(fits[2]))
        axes_loss.plot(xloss, loss_density_nodup,
                       label='No dup %s' % fmt_params(fits[4]))
        axes_loss.plot(xloss, loss_density_dup,
                       label='dup>0 %s' % fmt_params(fits[3]))
        ax1top.legend()
        ax1top.set_title('Loss rates')
        fig.suptitle('Fitting %s' % distrib.name)
        logger.info('Fitted %s' % (distrib.name))

        yield fig, kl, nlL, aic


def load_generax_familyrates():
    #workdir = op.expanduser('~/ws7/DUPLI_data93/alignments_analysis/duprates')
    #infile = op.join(workdir, 'familyrates.txt')
    workdir = Path.home().joinpath('ws7', 'DUPLI_data93', 'alignments_analysis', 'duprates')
    infile = workdir / 'familyrates.tsv'

    df = pd.read_csv(str(infile), sep='\t', index_col=0, header=None,
                     names=['nb', 'subgenetree', 'duprate', 'lossrate',
                            'phy_lL', 'rec_lL',
                            'S', 'SL', 'D', 'T', 'TL', 'Leaf', 'Invalid'],
                     dtype={'nb': str, 'subgenetree': str,
                            'duprate': float, 'lossrate': float,
                            'phy_lL': float, 'rec_lL': float,
                            'S': int, 'SL': int, 'D': int, 'T': int, 'TL': int,
                            'Leaf': int, 'Invalid': int})
    logger.info('Loaded data (shape %s)', df.shape)
    return df


def filter_invalid_generax(df, out=None):
    invalid = (df.duprate < 0) | (df.lossrate < 0)
    print('%d Invalid dup/loss rates' % invalid.sum(), file=out)
    df = df[~invalid].copy()
    print('%d NA values' % df.isnull().any(axis=1).sum(), file=out)
    # There's a minimum value returned by generax, and many gene trees have this value:
    min_dup = df.duprate.min()
    nodup = (df.duprate == min_dup)
    print('Min Dup rate = %g  (%d gene trees)' % (min_dup, nodup.sum()), file=out)
    min_loss = df.lossrate.min()
    print('Min loss rate = %g  (%d gene trees)' % (min_loss, (df.lossrate == min_loss).sum()), file=out)
    return df


def analysis_1_generax_withgamma():
    # First analyse output from GeneRax (1by1)
    with HtmlReport(str(workdir / 'familyrates.html'), style=css_dark_style) as hr:
        # Invalid values:
        df = filter_invalid_generax(load_generax_familyrates(), out=hr)

        fig, ((ax0top, ax1top), (ax0bottom, ax1bottom)) = plt.subplots(2, ncols=2)
        # Y axis with a broken scale : ---//---
        fig.subplots_adjust(hspace=0.1)
        axes_dup = brokenAxes(ax0bottom, ax0top)
        (heights, _, _), _ = axes_dup.hist(df.duprate, bins=100, density=True)
        axes_dup.dobreak(max(heights[1:]))

        logger.info('Plotted dup rates')

        axes_loss = brokenAxes(ax1bottom, ax1top)
        ((_,heights), _, _), _ = axes_loss.hist((df.lossrate[nodup], df.lossrate[~nodup]),
                                                label=['No dup', 'dup > 0'],
                                                bins=100, density=True, stacked=True)
        ymax = heights.max()
        axes_loss.dobreak(max(heights[heights<ymax]))
        logger.info('Plotted loss rates')

        dup_gamma = gamma.fit(df.duprate.values)
        dup_gamma_nonzero = gamma.fit(df.duprate[~nodup].values)
        loss_gamma = gamma.fit(df.lossrate.values)
        loss_gamma_nodup = gamma.fit(df.lossrate[nodup].values)
        loss_gamma_dup = gamma.fit(df.lossrate[~nodup].values)

        xdup = np.linspace(*ax0bottom.get_xlim(), num=100)
        xloss = np.linspace(*ax1bottom.get_xlim(), num=100)
        axes_dup.plot(xdup, gamma.pdf(xdup, *dup_gamma),
                      label='Gamma(%g, %g, %g)' % dup_gamma)
        #ax0.annotate('Gamma(%g, %g, %g)' % dup_gamma, (1,1), (-2,-2),
        #             xycoords='axes fraction', textcoords='offset points',
        #             va='top', ha='right')
        axes_dup.plot(xdup, gamma.pdf(xdup, *dup_gamma_nonzero),
                      label='dup>0 Gamma(%g, %g, %g)' % dup_gamma_nonzero)
        ax0top.legend()
        ax0bottom.set_ylabel('% of gene trees')
        ax0top.set_title('Duplication rates')

        axes_loss.plot(xloss, gamma.pdf(xloss, *loss_gamma),
                       label='Gamma(%g, %g, %g)' % loss_gamma)
        axes_loss.plot(xloss, gamma.pdf(xloss, *loss_gamma_nodup),
                       label='No dup Gamma(%g, %g, %g)' % loss_gamma_nodup)
        axes_loss.plot(xloss, gamma.pdf(xloss, *loss_gamma_dup),
                       label='dup>0 Gamma(%g, %g, %g)' % loss_gamma_dup)
        #ax1.annotate('Gamma(%g, %g, %g)' % loss_gamma, (1,1), (-2,-2),
        #             xycoords='axes fraction', textcoords='offset points',
        #             va='top', ha='right')
        ax1top.legend()
        ax1top.set_title('Loss rates')
        logger.info('Fitted gamma distribs')

        fig.savefig(str(workdir / 'hist_duploss.pdf'), bbox_inches='tight')
        hr.show()

        for dataset, data_params, data in zip(
            ('dup', 'dup>0', 'loss', 'nodup loss', 'dup>0 loss'),
            (dup_gamma, dup_gamma_nonzero, loss_gamma, loss_gamma_nodup, loss_gamma_dup),
            (df.duprate, df.duprate[~nodup], df.lossrate, df.lossrate[nodup], df.lossrate[~nodup])):
            hr.mkd(r'\\[ - lL(%s\\ params) = %g \\]' % (
                      dataset,
                      gamma.nnlf(inverse_scale(data_params), data.values)))

        nonzerodup_df = df.loc[~nodup]
        suggested_transforms = test_transforms(df, ['duprate', 'lossrate'],
                                               widget=False, out=hr)
        suggested_transforms_nonzerodup = test_transforms(nonzerodup_df, ['duprate', 'lossrate'],
                                               widget=False, out=hr)

        fig, axes = plt.subplots(ncols=2, squeeze=False)
        scatter_density('duprate', 'lossrate', data=df, alpha=0.4, ax=axes[0,0])
        axes[0,0].set_xlabel('dup rate')
        axes[0,0].set_ylabel('loss rate')
        axes[0,0].set_title('Complete dataset', usetex=True)
        scatter_density('duprate', 'lossrate', data=nonzerodup_df, alpha=0.4, ax=axes[0,1])
        axes[0,1].set_xlabel('dup rate')
        axes[0,1].set_ylabel('loss rate')
        axes[0,1].set_title('Dataset with $\delta > 0$', usetex=True)
        hr.show()
        fig.savefig(str(workdir/'scatter_duploss.pdf'), bbox_inches='tight')
        logger.info('Plotted scatter of dup/loss (dup>0, untransformed)')

        hr.print('\n# Regress untransformed data (dup>0)')
        fit = sm.OLS(nonzerodup_df[['lossrate']],
                     sm.add_constant(nonzerodup_df[['duprate']])).fit()
        hr.output(fit.summary())

        hr.print('\n# Regress transformed data (dup>0):')
        for ft, transform in suggested_transforms_nonzerodup.items():
            hr.print('- %s : %s' %(ft, transform.__name__))

        transformed_df = nonzerodup_df.transform(suggested_transforms_nonzerodup)
        fig, ax = plt.subplots()
        scatter_density('duprate', 'lossrate', data=transformed_df, alpha=0.4, ax=ax)
        ax.set_ylabel('%s(loss rate)' % suggested_transforms_nonzerodup['lossrate'].__name__)
        ax.set_xlabel('%s(dup rates)' % suggested_transforms_nonzerodup['duprate'].__name__)
        hr.show()
        fig.savefig(str(workdir/'scatter_transformed-duploss.pdf'), bbox_inches='tight')
        logger.info('Plotted scatter of dup/loss (dup>0, transformed:%s,%s)',
                    suggested_transforms_nonzerodup['duprate'].__name__,
                    suggested_transforms_nonzerodup['lossrate'].__name__
                    )

        fit = sm.OLS(transformed_df[['lossrate']],
                     sm.add_constant(transformed_df[['duprate']])).fit()
        hr.output(fit.summary())


def analysis_2_alldistribs():
    with HtmlReport(str(workdir / 'familyrates_fitdistribs.html'), style=css_dark_style) as hr:
        # Invalid values:
        df = filter_invalid_generax(load_generax_familyrates(), out=hr)
        hr.mkd('## Fitting other distributions')

        descriptions = ('dup', 'dup (>0)', 'loss', 'loss (dup>0)', 'loss (dup=0)')
        #tested_distribs = (gamma, st.expon) #st.dgamma, st.exponweib),
        tested_distribs = continuous_positive_distribs
        #divide by zero: betaprime, exponweib, genpareto, invweibull, loglaplace, mielke, nakagami, ncf
        #invalid value: alpha
        # RuntimeWarning non-integer value: Erlang.
        # Problems with Levy. (Inf KL)

        # Interesting distribs: exp, hyper exponential,
        # beta prime?
        # log-logistique (=Fisk)
        # Gompertz (age-dependent proba of dying)

        all_kl, all_nlL, all_aic = [], [], []
        with PdfPages('distrib_fits_duploss.pdf') as pdfdoc:
            for fig, kl, nlL, aic in fit_distribs(tested_distribs,
                                             df, x=['duprate', 'lossrate'], stream=hr):
                hr.show(fig)
                pdfdoc.savefig(fig, bbox_inches='tight', facecolor='k', transparent=True)
                plt.close()
                all_kl.append(kl)
                all_aic.append(aic)
                hr.flush()
        #logger.info('Saved all histograms.')

        kl_df = pd.DataFrame(all_kl, columns=descriptions, index=[d.name for d in tested_distribs])
        aic_df = pd.DataFrame(all_aic, columns=descriptions, index=[d.name for d in tested_distribs])

        fit_stats = pd.concat((kl_df, aic_df), axis=1, keys=['KL', 'AIC'],
                              verify_integrity=True)
        #fit_stats.format({'KL': '{:.5f}', 'AIC': '{:.2f}'})
        hr.html(fit_stats.style.apply(bounded_background_gradient, cmap='bone_r',
                                      subset=['KL'], bad='#707099', under='k', over='k')\
                               .apply(bounded_background_gradient, cmap='gray_r',
                                      subset=['AIC'], bad='#707099', under='k', over='k'))

    # CCL: Best ones:
    # - dup:    + chi, gengamma, pareto, powerlognorm, recipinvgauss (good KL, badAIC), ~expon
    #           - gompertz looks good but has high KL and high AIC. Nakagami looks good.
    #           - gamma is quite bad!
    # - dup >0:
    # NOTE:    
    # - do we really want truncexpon...?
    # - gengamma > {weibull, gamma}

    fit_stats.to_csv('distrib_fits_duploss.csv', sep='\t')
    #TODO: analyse the eventCounts as well.

def analysis_3_alldistribs():
    """2020/02/24"""
    with HtmlReport(str(workdir / 'familyrates_correlates.html'),
                    style=css_dark_style, css=[myslideshow_css],
                    scripts=[myslideshow_js]) as hr:
        df = filter_invalid_generax(load_generax_familyrates(), out=hr)
        stattypes = ('tree', 'al', 'codeml', 'cleaning', 'alfsa',
                     'codemlfsa', 'cleaningfsa', 'beastS', 'alfsahmmc', 'codemlfsahmmc')
        ss = load_subtree_stats(
                "../subtrees_stats/subtreesGoodQualO2_{stattype}stats-Simiiformes.tsv",
                stattypes=stattypes)
        check_subtree_stats(list(zip(stattypes, ss)))
        hr.show()
        plt.close()

        reg = fullRegression()





if __name__ == '__main__':
    #analysis_1_generax_withgamma()
    analysis_2_alldistribs()
