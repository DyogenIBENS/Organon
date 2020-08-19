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
from datasci.dataframe_recipees import bounded_background_gradient, \
                                       matplotlib_background_gradient
from datasci.savior import HtmlReport, css_dark_style, reroute_loggers, \
                           generate_slideshow, slideshow_generator

from genchron.analyse.regress_dating_errors import *
import genchron.analyse.regress_dating_errors as aregr

from LibsDyogen import myPhylTree


float_ = np.longdouble  # Global precision set here.


import logging
logger = logging.getLogger(__name__)

plt.style.use('softer')
#mpl.rcParams['text.usetex'] = True  # matplotlib can render some math without calling latex.
#mpl.rcParams['mathtext.fontset'] = 'bch'
#mpl.rcParams['mathtext.rm'] = 'Bitstream Charter Sans'
#mpl.rcParams['mathtext.it'] = 'Bitstream Charter Sans:italic'
#mpl.rcParams['mathtext.bf'] = 'Bitstream Charter Sans:bold'
#mpl.rcParams['font.family'] = 'Bitstream Charter'

# Check available fonts with:
#sorted([f.name for f in mpl.font_manager.fontManager.ttflist])
# Install ttf font files in .local/lib/python3.5/site-packages/matplotlib/mpl-data/fonts/ttf/

#TODO: set background of figure to dark. set better color cycle.
import seaborn as sb
#sb.set_palette('deep')
#sb.set_palette('pastel')
#sb.set_palette('Set2')

def inverse_scale(params):
    """Switch between scale and inverse_scale parametrization:
        return (params[0], params[1], 1/params[2])"""
    return (params[0], params[1], 1/params[2])

# For page width 18.5 cm, fontsize 9, fig ratio 0.6
thesisfigsize = (9.71, 5.85)
mpl.rcParams['figure.figsize'] = thesisfigsize
#myslideshow_js = Path.home().joinpath('mydvpt', 'atavistic-doc-tools', 'myslideshow.js') 
#myslideshow_css = Path.home().joinpath('mydvpt', 'atavistic-doc-tools', 'myslideshow.css') 

# Run from '~/ws7/DUPLI_data93/alignments_analysis/duprates/'

#workdir = op.expanduser('~/ws7/DUPLI_data93/alignments_analysis/duprates')
workdir = Path.home().joinpath('ws7', 'DUPLI_data93', 'alignments_analysis', 'duprates')

# Important: keep relative paths here.
myslideshow_js = Path('../../../../').joinpath('mydvpt', 'atavistic-doc-tools', 'myslideshow.js') 
myslideshow_css = Path('../../../../').joinpath('mydvpt', 'atavistic-doc-tools', 'myslideshow.css') 

continuous_distribs = [d for dname, d in sorted(vars(st).items())
                       if isinstance(d, st.rv_continuous)]
continuous_positive_distribs = [d for d in continuous_distribs if d.a >= 0 and np.isinf(d.b)]
# There are 49 or 50 distributions implemented.

restricted_distribs = [
st.betaprime,  # default loc = 0
st.chi,
st.chi2,
st.expon,
st.exponpow,
st.fisk,
st.foldcauchy,
st.foldnorm,
#st.halfcauchy,  # Special case of foldcauchy
#st.halfnorm,    # Special case of foldnorm
st.gompertz,
st.gamma,
st.gengamma,     # > weibull, gamma, expon, halfnorm...
st.invgamma,
st.invgauss,
st.pareto,
st.genpareto,
st.weibull_min,
st.exponweib,
st.invweibull    # = Fréchet
]
# Check all default loc parameter values (expect zero):
#getloc_reg = re.compile(r'^\s+pdf\(.*(loc=[^,]+).*\)', re.M)
#for d in restricted_distribs:
#    print(d.name, getloc_reg.search(d.__doc__).group(1))



# Interesting distribs: exp, hyper exponential (kind of mixture of expon),
# beta prime?
# log-logistique (=Fisk)
# Gompertz (age-dependent proba of dying)

distrib_realnames = {'betaprime': 'Beta prime',
    'chi': 'Chi',
    'chi2': 'Chi²',
    'expon': 'Exponential',
    'exponpow': 'Exponential power',  # (Generalised Gaussian)
    'fisk': 'Log-logistic (Fisk)',
    'foldcauchy': 'Folded Cauchy',
    'foldnorm': 'Folded normal',
    'halfcauchy': 'Half Cauchy',
    'halfnorm': 'Half normal',
    'gompertz': 'Gompertz',
    'gamma': 'Gamma',
    'gengamma': 'Generalised Gamma',
    'invgamma': 'Inverted Gamma',
    'invgauss': 'Inverse Gaussian',
    'pareto': 'Pareto',
    'genpareto': 'Generalised Pareto',
    'weibull_min': 'Weibull',
    'exponweib': 'Exponential Weibull',
    'invweibull': 'Inverted Weibull'
    }
distrib_realnames_fr = {'betaprime': 'Beta prime',
    'chi': 'Chi',
    'chi2': 'Chi²',
    'expon': 'Exponentielle',
    'exponpow': 'Exponentielle puissance',  # (Generalised Gaussian)
    'fisk': 'Log-logistique (Fisk)',
    'foldcauchy': 'Cauchy repliée',
    'foldnorm': 'Normale repliée',
    'halfcauchy': 'Demi-Cauchy',
    'halfnorm': 'Demi-normale',
    'gompertz': 'Gompertz',
    'gamma': 'Gamma',
    'gengamma': 'Gamma généralisée',
    'invgamma': 'Gamma inverse',
    'invgauss': 'Gaussienne inverse',
    'pareto': 'Pareto',
    'genpareto': 'Pareto généralisée',
    'weibull_min': 'Weibull',
    'exponweib': 'Weibull exponentielle',
    'invweibull': 'Fréchet'
    }

distrib_symbols = {'betaprime': r"$\beta'$",
    'chi': r'$\chi$',
    'chi2': r'$\chi^2$',
    'expon': r'$\mathcal{E}$',  # Might require special latex package
    'exponpow': 'Expon. power',
    'fisk': 'Fisk',
    'foldcauchy': 'Folded Cauchy',
    'foldnorm': 'Folded normal',
    'halfcauchy': 'Half Cauchy',
    'halfnorm': 'Half normal',
    'gompertz': 'Gompertz',
    'gamma': r'$\Gamma$',
    'gengamma': r'Generalised $\Gamma$',
    'invgamma': r'$\mathrm{Inv-}\Gamma$',
    'invgauss': 'Inv-Gaussian',
    'pareto': 'Pareto',
    'genpareto': 'Generalised Pareto',
    'weibull_min': 'Weib',
    'exponweib': 'Exp Weib',
    'invweibull': 'Inv-Weib'
    }

#always_estimate_loc = set(('pareto',))  # Loc param > 0
strictly_positive_distribs = set((d for d in continuous_positive_distribs
                                  if d.a > 0))  # Includes only Pareto as of Scipy v1.4.1


def fmt_params(distrib, params, prec=4, scalename='s', floc=None, sep=', '):
    fmt='%%s=%%.%dg' % prec
    if distrib.numargs:
        param_names = distrib.shapes.split(', ') + ['loc', scalename]
    else:
        param_names = ['loc', scalename]
    if floc is not None:
        # Do not display the loc parameter since it is fixed.
        params = params[:-2] + params[-1:]
        del param_names[-2]
    return sep.join(fmt % px for px in zip(param_names, params))


def fmt_params_gamma(params, prec=4, floc=None, sep=', '):
    fmt='%%.%dg' % prec
    template = sep.join([r'$\alpha\!=\!{0}$', r'$\beta\!=\!{2}$'])
    if floc is None:
        template += sep + r'$\mathrm{{loc}}\!=\!{1}$'
    return template.format(*(fmt % p for p in inverse_scale(params)))


def fit_distribs(distribs, df, nodup, x=None, nbins=100, stream=None, lang_fr=True, dark=False,
                 fix_loc=0):
    """
    floc: fix the location parameter when fitting.
    """
    if dark: plt.style.use('softdark')  # Also see the default 'dark_background'
    dist_realnames = distrib_realnames_fr if lang_fr else distrib_realnames
    if x is None:
        x = df.columns.tolist()
    ncols = len(x)  # controling x is NotImplemented
    descriptions = ('dup', 'dup (>0)', 'loss', 'loss (dup>0)', 'loss (dup=0)')
    datas = (df.duprate.values,
            df.duprate[~nodup].values,
            df.lossrate.values,
            df.lossrate[~nodup].values,
            df.lossrate[nodup].values)
    hist_kwargs = dict(density=True, edgecolor='none', alpha=0.7, rwidth=1.05)  # stacked=True
    for i, distrib in enumerate(distribs, start=1):
        floc = fix_loc #floc = None if distrib in strictly_positive_distribs else fix_loc
        # TODO: fix loc for Pareto for datas[1]
        fit_kw = {} if floc is None else {'floc': floc}
        fmt_fit = (partial(fmt_params_gamma, floc=floc) if distrib.name == 'gamma'
                   else partial(fmt_params, distrib, floc=floc))
        distname = dist_realnames.get(distrib.name)
        distname = '%s (%s)' % (distname, distrib.name) if distname else distrib.name.capitalize()
        distrib_k = distrib.numargs+2 if floc is None else distrib.numargs+1
        # Followed by nb of params: numargs(shapes) + loc + scale.
        #if distrib.name == 'pareto':
        #    logger.debug('Pareto: fix_loc=%s floc=%s in_always_est=%s fit_kw=%s',
        #            fix_loc, floc, distrib.name in always_estimate_loc, fit_kw)
        
        print('### %s k=%d' % (distname, distrib_k), file=stream)
        fig, ((ax0top, ax1top), (ax0bottom, ax1bottom)) = plt.subplots(2, ncols=ncols)
        # Y axis with a broken scale : ---//---
        fig.subplots_adjust(hspace=0.1)
        axes_dup = brokenAxes(ax0bottom, ax0top)
        (h_dup, b_dup, _), _ = axes_dup.hist(df.duprate, bins=nbins,
                                             label=r'$\delta$ observ' + ('é' if lang_fr else 'ed'),
                                             **hist_kwargs)
        # Move the color cycler one-step forward:
        axes_dup.plot([])
        # overfitted break position: 3rd max height:
        axes_dup.dobreak(sorted(h_dup)[-3]*1.02)

        axes_loss = brokenAxes(ax1bottom, ax1top)
        ((h_loss,h_loss_nodup), b_loss, _), _ = axes_loss.hist(
                                                    (df.lossrate[~nodup],
                                                     df.lossrate[nodup]),
                                                #label=['No dup', 'dup > 0'],
                                    label=[
                                        r'$\lambda_{\delta\! >\! 0}$ observ' + ('é' if lang_fr else 'ed'),
                                        r'$\lambda_{\delta\! =\! 0}$ observ' + ('é' if lang_fr else 'ed')],
                                    bins=nbins, stacked=True,
                                    **hist_kwargs)
        ymax = h_loss.max()
        axes_loss.dobreak(max(h_loss[h_loss<ymax]))

        fits = [distrib.fit(data, **fit_kw) for data in datas]
        if len(fits[0]) != distrib.numargs + 2:
            logger.error('Unexpected difference in nb params: theo %d VS fitted %d',
                         distrib.numargs+2, len(fits[0]))

        # X coordinates of the middle of the bins:
        xdup = b_dup[:-1] + (b_dup[1] - b_dup[0])/2     #, dtype=float_
        xloss = b_loss[:-1] + (b_loss[1] - b_loss[0])/2 #, dtype=float_
        dup_density = distrib.pdf(xdup, *fits[0])
        dup_density_nonzero = distrib.pdf(xdup, *fits[1])
        axes_dup.plot(xdup, dup_density,
                      label='%s(%s)' % (distrib_symbols.get(distrib.name, distrib.name.capitalize()),
                                        fmt_fit(fits[0])))
        axes_dup.plot(xdup, dup_density_nonzero,
                      label=r'%s(%s);  $\delta\! >\! 0$' % (
                                distrib_symbols.get(distrib.name, distrib.name),
                                fmt_fit(fits[1])))
        # Overfitted legend position
        ax0top.legend(loc='upper left', bbox_to_anchor=(0.1, 1), fontsize='small')
        ax0bottom.set_ylabel("% d'arbres de gènes" if lang_fr else '% of gene trees')
        ax0top.set_title(('Taux de duplication' if lang_fr else 'Duplication rates') + r' $\delta$')

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
        aic = [2*distrib_k + 2*nll for nll in nlL]
        ks, ks_p = zip(*[st.kstest(data, distrib(*fit).cdf, N=nbins) for fit, data
                         in zip(fits, datas)])

        for info in zip(descriptions, kl, nlL, aic, ks, ks_p):
            print('%-14s: KL = %8.6f;  -lL = %8.1f;  AIC = %7.0f;  KS = %8.6g p=%-8g' % info, file=stream)

        axes_loss.plot(xloss, loss_density,
                       label=fmt_fit(fits[2]))
        axes_loss.plot(xloss, loss_density_dup,
                       label=r'%s;  $\delta\!>\!0$' % fmt_fit(fits[3]))
        axes_loss.plot(xloss, loss_density_nodup,
                       label=r'%s;  $\delta\!=\!0$' % fmt_fit(fits[4]))
        ax1top.legend(fontsize='small')
        ax1top.set_title(('Taux de perte' if lang_fr else 'Loss rates') + r' $\lambda$')
        fig.suptitle('Fitting %s' % distrib.name)
        logger.info('Fitted %-16s [%2d/%d]' % (distrib.name, i, len(distribs)))

        yield fig, kl, nlL, aic, ks, ks_p


def load_generax_familyrates():
    #infile = op.join(workdir, 'familyrates.txt')
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
    df['duprate_nonzero'] = (df.duprate.dropna()>1e-7)
    df['generax_robust'] = ~(df[['SL', 'D', 'T', 'TL']].dropna().any(axis=1))
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
    return df, nodup


def analysis_1_generax_withgamma(lang_fr=True, dark=False):
    # First analyse output from GeneRax (1by1)
    with HtmlReport(str(workdir / 'familyrates.html'),
                    style=(css_dark_style if dark else None)) as hr:
        # Invalid values:
        df, nodup = filter_invalid_generax(load_generax_familyrates(), out=hr)

        fig, ((ax0top, ax1top), (ax0bottom, ax1bottom)) = plt.subplots(2, ncols=2)
        # Y axis with a broken scale : ---//---
        fig.subplots_adjust(hspace=0.1)
        axes_dup = brokenAxes(ax0bottom, ax0top)
        (heights, _, _), _ = axes_dup.hist(df.duprate, bins=100, density=True)
        axes_dup.dobreak(max(heights[1:]))

        logger.info('Plotted dup rates')

        axes_loss = brokenAxes(ax1bottom, ax1top)
        ((_,heights), _, _), _ = axes_loss.hist((df.lossrate[nodup], df.lossrate[~nodup]),
                                                #label=['No dup', 'dup > 0'],
                                                label=[r'$\delta=0$', r'$\delta>0$'],
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
                      label=r'$\Gamma(%g, %g, %g)$' % inverse_scale(dup_gamma))
        #ax0.annotate('Gamma(%g, %g, %g)' % dup_gamma, (1,1), (-2,-2),
        #             xycoords='axes fraction', textcoords='offset points',
        #             va='top', ha='right')
        axes_dup.plot(xdup, gamma.pdf(xdup, *dup_gamma_nonzero),
                      #label='dup>0 Gamma(%g, %g, %g)' % dup_gamma_nonzero)
                      label=r'$\delta>0$ $\Gamma(%g, %g, %g)$' % inverse_scale(dup_gamma_nonzero))
        ax0top.legend()
        ax0bottom.set_ylabel("% d'arbres de gènes" if lang_fr else '% of gene trees')
        ax0top.set_title(('Taux de duplication' if lang_fr else 'Duplication rates') + r' $\delta$')

        axes_loss.plot(xloss, gamma.pdf(xloss, *loss_gamma),
                       label=r'$\Gamma(%g, %g, %g)$' % loss_gamma)
        axes_loss.plot(xloss, gamma.pdf(xloss, *loss_gamma_nodup),
                       label=r'$\delta=0$ $\Gamma(%g, %g, %g)$' % inverse_scale(loss_gamma_nodup))
        axes_loss.plot(xloss, gamma.pdf(xloss, *loss_gamma_dup),
                       label=r'$\delta>0$ $\Gamma(%g, %g, %g)$' % inverse_scale(loss_gamma_dup))
        #ax1.annotate('Gamma(%g, %g, %g)' % loss_gamma, (1,1), (-2,-2),
        #             xycoords='axes fraction', textcoords='offset points',
        #             va='top', ha='right')
        ax1top.legend()
        ax1top.set_title(('Taux de perte' if lang_fr else 'Loss rates') + r' $\lambda$')
        logger.info('Fitted gamma distribs')

        fig.savefig(str(workdir / ('fit_duploss_gamma%s.pdf' % ('_darkbg' if dark else ''))),
                    bbox_inches='tight')
        hr.show()

        for dataset, data_params, data in zip(
            #('dup', 'dup>0', 'loss', 'nodup loss', 'dup>0 loss'),
            (r'$\delta$', r'$\delta>0$', r'$\lambda$', r'$\lambda_{\delta=0}$', r'$\lambda_{\delta>0}$'),
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
        axes[0,0].set_xlabel(r'$\delta$')
        axes[0,0].set_ylabel(r'$\lambda$')
        axes[0,0].set_title('Données complètes' if lang_fr else 'Complete dataset')
        scatter_density('duprate', 'lossrate', data=nonzerodup_df, alpha=0.4, ax=axes[0,1])
        axes[0,1].set_xlabel(r'$\delta$')
        axes[0,1].set_ylabel(r'$\lambda$')
        axes[0,1].set_title(r'Données avec $\delta>0$' if lang_fr else r'Dataset with $\delta > 0$')
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
        ax.set_ylabel(r'$\mathrm{%s}(\lambda)$' % suggested_transforms_nonzerodup['lossrate'].__name__)
        ax.set_xlabel(r'$\mathrm{%s}(\delta)$' % suggested_transforms_nonzerodup['duprate'].__name__)
        hr.show()
        fig.savefig(str(workdir/'scatter_transformed-duploss.pdf'), bbox_inches='tight')
        logger.info('Plotted scatter of dup/loss (dup>0, transformed:%s,%s)',
                    suggested_transforms_nonzerodup['duprate'].__name__,
                    suggested_transforms_nonzerodup['lossrate'].__name__
                    )

        fit = sm.OLS(transformed_df[['lossrate']],
                     sm.add_constant(transformed_df[['duprate']])).fit()
        hr.output(fit.summary())


def analysis_2_alldistribs(lang_fr=True, dark=False, restrict_distribs=True,
        fix_loc=0):
    filesuffix = '_restricted' if restrict_distribs else ''
    if fix_loc is None: filesuffix += '_nofloc'
    dist_realnames = distrib_realnames_fr if lang_fr else distrib_realnames
    pd.set_option('display.show_dimensions', True)
    with HtmlReport(str(workdir / 'familyrates_fitdistribs%s.html') % filesuffix,
                    style=(css_dark_style if dark else None),
                    css=[myslideshow_css], scripts=[myslideshow_js]) as hr:
        # Invalid values:
        df, nodup = filter_invalid_generax(load_generax_familyrates(), out=hr)
        hr.mkd('## Fitting other distributions')

        descriptions = ('dup', 'dup (>0)', 'loss', 'loss (dup>0)', 'loss (dup=0)')
        descriptions_tex = (r'$\delta$', r'$\delta > 0$',
                r'$\lambda$', r'$\lambda_{\delta\! >\! 0}$', r'$\lambda_{\delta\! =\! 0}$')
        #tested_distribs = (gamma, st.expon) #st.dgamma, st.exponweib),
        tested_distribs = restricted_distribs if restrict_distribs else continuous_positive_distribs
        #divide by zero: betaprime, exponweib, genpareto, invweibull, loglaplace, mielke, nakagami, ncf
        #invalid value: alpha
        # RuntimeWarning non-integer value: Erlang.
        # Problems with Levy. (Inf KL)

        # Interesting distribs: exp, hyper exponential,
        # beta prime?
        # log-logistique (=Fisk)
        # Gompertz (age-dependent proba of dying)

        all_kl, all_nlL, all_aic, all_ks, all_ks_p = [], [], [], [], []
        with PdfPages('distrib_fits_duploss%s.pdf' % filesuffix) as pdfdoc:
            for _, (fig, kl, nlL, aic, ks, ks_p) in zip(
                    generate_slideshow(hr, len(tested_distribs)),
                    fit_distribs(tested_distribs,
                                         df, nodup, x=['duprate', 'lossrate'],
                                         stream=hr, lang_fr=lang_fr, dark=dark,
                                         fix_loc=fix_loc)):
                hr.show(fig)
                pdfdoc.savefig(fig, bbox_inches='tight', transparent=True,
                               facecolor=('k' if dark else 'none'))
                plt.close()
                all_kl.append(kl)
                all_aic.append(aic)
                all_ks.append(ks)
                all_ks_p.append(ks_p)
                hr.flush()
        #logger.info('Saved all histograms.')

        kl_df = pd.DataFrame(all_kl, columns=descriptions_tex,
                             index=[dist_realnames.get(d.name, d.name.capitalize()) for d in tested_distribs])
        aic_df = pd.DataFrame(all_aic, columns=descriptions_tex,
                              index=[dist_realnames.get(d.name, d.name.capitalize()) for d in tested_distribs])
        ks_df = pd.DataFrame(all_ks, columns=descriptions_tex,
                             index=[dist_realnames.get(d.name, d.name.capitalize()) for d in tested_distribs])
        ks_p_df = pd.DataFrame(all_ks_p, columns=descriptions_tex,
                             index=[dist_realnames.get(d.name, d.name.capitalize()) for d in tested_distribs])

        fit_stats = pd.concat((kl_df, aic_df), axis=1, keys=['KL', 'AIC'],
                              verify_integrity=True)
        fit_stats.rename(dict(zip(descriptions_tex, descriptions)), axis=1, level=1).to_csv('distrib_fits_duploss%s.csv' % filesuffix, sep='\t')

        ks_stats = pd.concat((ks_df, ks_p_df), axis=1, keys=['stat', 'Pval'],
                             verify_integrity=True)
        ks_stats.rename(dict(zip(descriptions_tex, descriptions)), axis=1, level=1).to_csv('distrib_KS_duploss%s.csv' % filesuffix, sep='\t')

        #fit_stats.format({'KL': '{:.5f}', 'AIC': '{:.2f}'})
        logger.info('fit_stats.columns = %s', fit_stats.columns)
        if fit_stats.index.has_duplicates:
            logger.error('fit_stats.index has duplicates!')


        cmapKL = plt.get_cmap('bone_r')
        cmapKL.set_bad('#707099')
        cmapKL.set_under('k')
        cmapKL.set_over('k')
        cmapAIC = plt.get_cmap('gray_r')
        cmapAIC.set_bad('#707099')
        cmapAIC.set_under('k')
        cmapAIC.set_over('k')

        hr.html(fit_stats.rename(lambda n: n.replace('$', '$$'), axis=1, level=1)\
                        .style.apply(bounded_background_gradient, cmap=cmapKL,
                                     subset=['KL'])\
                              .apply(bounded_background_gradient, cmap=cmapAIC,
                                     subset=['AIC']))
        hr.raw('\n<br />\n')
        hr.html(ks_stats.rename(lambda n: n.replace('$', '$$'), axis=1, level=1)\
                        .style.apply(bounded_background_gradient, cmap=cmapKL,
                                     subset=['stat'])\
                              .apply(bounded_background_gradient, cmap=cmapAIC,
                                     subset=['Pval']))
        #fig, (axKL, axAIC) = plt.subplots(ncols=2, hspace=0, sharey=True)
        #print("fit_stats.xs('KL', axis=1) =\n" + str(fit_stats.xs('KL', axis=1)))
        fig, axKL = plt.subplots()
        _, cbar = matplotlib_background_gradient(
                    #fit_stats.transpose().xs('KL'),
                    fit_stats.xs('KL', axis=1),
                    cmapKL, axis=0, ax=axKL)
        #axKL.set_title('Kullback-Leibler divergences from observed distribution')
        fig.set_size_inches((5, 7) if restrict_distribs else (5, 12))
        fig.savefig(str(workdir / 'distrib_KL_duploss%s.pdf') % filesuffix, bbox_inches='tight')
        fig.savefig(str(workdir / 'distrib_KL_duploss%s.png') % filesuffix, bbox_inches='tight')
        hr.mkd('\n![Kullback-Leibler divergences](%s)\n' % ('distrib_KL_duploss%s.png' % filesuffix))

        fig, axKS = plt.subplots()
        _, cbar = matplotlib_background_gradient(
                    #fit_stats.transpose().xs('KL'),
                    ks_stats.xs('stat', axis=1),
                    cmapKL, axis=0, ax=axKS)
        axKS.set_title('Kolmogorov-Smirnov statistic')
        fig.set_size_inches((5, 7) if restrict_distribs else (5, 12))
        fig.savefig(str(workdir / 'distrib_KS_duploss%s.pdf') % filesuffix, bbox_inches='tight')
        fig.savefig(str(workdir / 'distrib_KS_duploss%s.png') % filesuffix, bbox_inches='tight')
        hr.mkd('\n![Kolmogorov-Smirnov statistic](%s)\n' % ('distrib_KS_duploss%s.png' % filesuffix))
    
    # CCL: Best ones:
    # - dup:    + chi, gengamma, pareto, powerlognorm, recipinvgauss (good KL, badAIC), ~expon
    #           - gompertz looks good but has high KL and high AIC. Nakagami looks good.
    #           - gamma is quite bad!
    # - dup >0:
    # NOTE:    
    # - do we really want truncexpon...?
    # - gengamma > {weibull, gamma}
    #
    #TODO: analyse the eventCounts as well.

common_info = ['genetree']

stat_params = dict(
    al=["ingroup_glob_len", "ingroup_mean_GC", "ingroup_mean_N",
        "ingroup_mean_gaps", "ingroup_mean_CpG", "ingroup_std_len",
        "ingroup_std_GC",  "ingroup_std_N",  "ingroup_std_gaps",
        "ingroup_std_CpG", "ingroup_nucl_entropy_mean",
        "ingroup_nucl_entropy_median", "ingroup_nucl_entropy_std",
        "ingroup_nucl_parsimony_mean", "ingroup_nucl_parsimony_median",
        "ingroup_nucl_parsimony_std", "ingroup_codon_entropy_mean",
        "ingroup_codon_entropy_median", "ingroup_codon_entropy_std",
        "ingroup_codon_parsimony_mean", "ingroup_codon_parsimony_median",
        "ingroup_codon_parsimony_std"],
    tree=["root_location", "really_robust", "nodes_robust", "single_child_nodes",
          "aberrant_dists", "rebuilt_topo", "unlike_clock",
          "consecutive_zeros", "sister_zeros", "triplet_zeros",
          "bootstrap_min", "bootstrap_mean"],  #, "root2tip_mean", "root2tip_sd"  (Already in cl_params)
    cl=["ls", "ns", "Nbranches", "NnonsynSites", "NsynSites",
        "kappa", "prop_splitseq", "codemlMP", "convergence_warning",
        "treelen", #"dN_treelen", "dS_treelen",
        "brOmega_mean", "brOmega_std", "brOmega_med", "brOmega_skew",
        "consecutive_zeros_t", "sister_zeros_t", "triplet_zeros_t",
        "consecutive_zeros_dS","sister_zeros_dS", "triplet_zeros_dS",
        "consecutive_zeros_dN", "sister_zeros_dN", "triplet_zeros_dN",
        "lnL", "Niter", "seconds"] \
       + ['r2t_%s_%s' %(m,s) for m in ('t', 'dS', 'dN')
          for s in ('mean', 'std')],
    cleaning=["gb_Nblocks", "gb_percent", "hmmc_propseqs", "hmmc_max"], #, "hmmc_mean_onlycleaned"]
    chronos=['%s_%s' % (method, info)
             for method in ('P1', 'P100', 'P10000', 'R1', 'R100', 'R10000', 'C')
             for info in ('message', 'ploglik', 'loglik', 'PHIIC')]
    )
stat_params['alfsahmmc'] = stat_params['alfsa'] = stat_params['al']
stat_params['codemlfsahmmc'] = stat_params['codemlfsa'] = stat_params['codeml'] = stat_params['cl']
stat_params['cleaningfsa'] = stat_params['cleaning']

dataset_params = ["freq_null_dist", "freq_null_t", "freq_null_dS", "freq_null_dN",
                  "null_dist_before", "null_t_before", "null_dS_before", "null_dN_before",
                  "null_dist_after", "null_t_after", "null_dS_after", "null_dN_after"]
                 # + ["Ndup", "Nspe"] 
rate_params = [statname % m for statname in ('%s_rate', '%s_rate_std') for m in MEASURES]
dataset_params_dS = ['freq_null_dS', 'null_dS_before', 'null_dS_after']
rate_params_dS = ['dS_rate_local', 'dS_rate_std_local', 'dS_rate_nonlocal', 'dS_rate_std_nonlocal']


workspace = Path.home() / 'ws7'
phyltree = myPhylTree.PhylogeneticTree(str(workspace / 'DUPLI_data93/PhylTree.TimeTree201901.Ensembl93-like.goodQual.nwk'))
timetree_ages_CI = pd.read_csv(str(workspace / 'databases/timetree/Primates_conf-int_201909.txt'),
                       sep='\s+', header=None, index_col=0,
                       names=['taxon', 'timetree_CI95_inf', 'timetree_CI95_sup'])#.rename_axis(index='taxon')
dosreis_ages_CI = pd.read_csv(
        str(workspace / 'databases/DosReis2018_Primates-timescale/Primates_dates.tsv'),
        sep='\t', header=0, index_col=0, keep_default_na=False)

ordered_simii_anc = ['Platyrrhini',
                     'Cebidae',
                     'Catarrhini',
                     'Cercopithecidae',
                     'Cercopithecinae',
                     'Papionini',
                     'Macaca',
                     'Hominoidea',
                     'Hominidae',
                     'Homininae',
                     'HomoPan',
                     'Pan']
extra_primates = ['Simiiformes', 'Primates', 'Strepsirrhini', 'Lemuriformes']
calibrations = pd.concat((pd.Series([phyltree.ages[anc] for anc in ordered_simii_anc+extra_primates],
                                    index=ordered_simii_anc+extra_primates,
                                    name='timetree_age'),
                          timetree_ages_CI,
                          dosreis_ages_CI.drop('synonym', 1)
                 ),
                 axis=1, sort=False)

renames = {}
# feature_longnames_fr.tsv
with open(str(workdir.parent / 'subtrees_stats' / 'feature_longnames.tsv')) as f:
    for line in f:
        ft, longname, *_ = line.rstrip().split('\t')
        renames[ft] = longname

def analysis_3_regress_duprates(lang_fr=True, dark=False):
    """2020/02/24"""
    loggers[0].setLevel(logging.DEBUG)
    with reroute_loggers(
            HtmlReport(str(workdir / 'familyrates_correlates.html'),
                       style=(css_dark_style if dark else None), css=[myslideshow_css],
                       scripts=[myslideshow_js], external_content=True),
            loggers) as hr:
        df, nodup = filter_invalid_generax(load_generax_familyrates(), out=hr)

        all_stattypes = ('tree', 'al', 'codeml', 'cleaning', 'alfsa',
                         'codemlfsa', 'cleaningfsa', 'beastS', 'alfsahmmc', 'codemlfsahmmc')
        stattypes = ('tree', 'alfsa', 'codemlfsa', 'cleaning')
        ss = load_subtree_stats(
                "../subtrees_stats/subtreesGoodQualO2_{stattype}stats-Simiiformes.tsv",
                stattypes=stattypes)
        ax = check_subtree_stats(list(zip(stattypes, ss)))
        hr.show(ax.figure); plt.close()

        features = [ft for stype, sdata in zip(stattypes, ss)
                    for ft in stat_params.get(stype,
                                              sdata.drop(common_info, 1).columns.tolist())]

        # I just need Ndup and Nspe
        #for fsa, there are only the robust...
        #ages_file = '../ages/Simiiformes_m1w04_ages.subtreesGoodQualO2-ci-um1.tsv.bz2'
        ages_file = '../ages/Simiiformes_m1w04-fsa_ages.subtreesGoodQualO2-ci-um1.tsv.bz2'
        #ages_treestats, ns = load_prepare_ages(ages_file, ss[0], measures=['dS'])

        print('\n# Loading ages: m1w04-fsa -c%s-%s' %('i', 'um1'), file=hr)
        analysis_fsa = analyse_age_errors(ages_file,
                                        'Simiiformes', phyltree,
                                        calibrations, ss[0], measures=['dS'], #'dN', 't', 'dist'],
                                        control_names=['timetree', 'dosreis'],
                                        out=hr)
        #hr.html(analysis_fsa.control_brlen)
        analysis_fsa.ages_controled_withnonrobust = analysis_fsa.ages_controled_withnonrobust.join(
                pd.read_csv('../ages/Simiiformes_m1w04-fsa-hmmc_star-branch-dS.subtreesGoodQualO2.tsv.bz2',
                            sep='\t', index_col=0), sort=False)
    
        dataset_params_dS = ['freq_null_dS', 'null_dS_before', 'null_dS_after', 'Ndup', 'Nspe']
        alls = pd.concat(([df.set_index('subgenetree')] #, ns[dataset_params_dS]]
                          + [sdata[stat_params.get(stype,
                                          sdata.columns.difference(common_info))]
                             for stype, sdata in zip(stattypes, ss)]),
                         axis=1, join='outer', sort=False, verify_integrity=True)
        
        features.remove('really_robust')
        features.remove('nodes_robust')
        features.remove('single_child_nodes')
        features.remove('ns')
        features.remove('Nbranches')
        # Also: Niter correlates with sitelnL and Ringroup_nucl_parsimony_std
        hr.html(alls.groupby(['duprate_nonzero', 'really_robust'])\
                .size().unstack()) #.style.caption('Comparison between generax and treebest robustness.'))
        hr.html(alls\
                .assign(Ndup_nonzero=(analysis_fsa.ns.Ndup.dropna()>0))\
                .groupby(['duprate_nonzero', 'Ndup_nonzero'])\
                .size().unstack()) #.style.caption('Comparison between generax and treebest duplications.'))

        hr.mkd('# Regression')
        try:
            #reg = fullRegression(
            #        alls,
            #        ['duprate', 'lossrate'],
            #        features,
            #        ref_suggested_transform=None,
            #        impose_transform=aregr._must_transform,
            #        #must_drop_features=aregr._must_drop_features.union(
            #        #    ('single_child_nodes',)),
            #        to_decorr=aregr._must_decorr,
            #        protected_features=aregr._protected_features,
            #        must_drop_data=aregr._must_drop_data,  # TODO: without dropping data.
            #        out=hr,
            #        logger=loggers[0],
            #        widget=slideshow_generator(hr)
            #        )
            #TODO:
            reg = full_dating_regression(analysis_fsa, alls, dataset_params_dS,
                    ['duprate', 'lossrate'], features,
                    ref_suggested_transform=None,
                    impose_transform=aregr._must_transform,
                    must_drop_features=aregr._must_drop_features.union(('Niter',)),
                    to_decorr=aregr._must_decorr,
                    protected_features=aregr._protected_features,
                    #must_drop_data=aregr._must_drop_data,  # TODO: without dropping data.
                    out=hr,
                    logger=loggers[0],
                    widget=slideshow_generator(hr)
                    )
            reg.do_rates(dict(dist_measures=['star_branch_dS'],
                         branchtime='star_median_brlen_dS',
                         mean_condition='type!="dup"',
                         weighted=True))  # Approx method
            print('## Rates computed (spe2spe approx).', file=hr)
            print(reg.cs_rates.head(), file=hr)
            coefs = reg.do()
            fig, axes = reg.plot_coefs(renames=renames)
            outfig = 'familyrates_correlates_coefs.pdf'
            fig.savefig(str(workdir / outfig), bbox_inches='tight')
            plt.close()
            hr.mkd('\n![Tree features coefficients of regression](%s)\n' % outfig)
        finally:
            return reg



if __name__ == '__main__':
    try:
        from UItools import colorlog
        colorlog.ColoredFormatter.install(logger)
    except ImportError:
        logging.basicConfig()
    logger.setLevel(logging.DEBUG)
    logging.getLogger('datasci.dataframe_recipees').setLevel(logging.DEBUG)

    #analysis_1_generax_withgamma()
    #analysis_2_alldistribs()
    # 2020/02/24:
    #reg = analysis_3_regress_duprates()

    # 2020/05/20:
    analysis_2_alldistribs()
    # 2020/05/22 11:11
    #analysis_2_alldistribs(restrict_distribs=False, fix_loc=None)
    # 2020/05/22 11:30
    #analysis_2_alldistribs(restrict_distribs=False)
