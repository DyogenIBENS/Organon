#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import statsmodels.api as sm
import numpy as np
import pandas as pd
from statsmodels.stats.outliers_influence import variance_inflation_factor
from datasci.stats import VIF

df = pd.DataFrame(np.array(
    [[0.48271507, 0.99050744, 0.55593631],
     [0.64841363, 0.42579384, 0.17202914],
     [0.42471806, 0.86893611, 0.20153212],
     [0.14127473, 0.51206951, 0.16471164],
     [0.30567444, 0.11130688, 0.34763897]]), columns=list('abc'))

def test_vif():
    model = sm.OLS(df.a, sm.add_constant(df[['b', 'c']]))
    fit = model.fit()
    vif_const = VIF(fit, 'const')
    vif_b = VIF(fit, 'b')
    vif_c = VIF(fit, 'c')
    ref_vif_const = variance_inflation_factor(model.exog, 0)
    ref_vif_b = variance_inflation_factor(model.exog, 1)
    ref_vif_c = variance_inflation_factor(model.exog, 2)
    assert np.isclose(vif_b, ref_vif_b)
    assert np.isclose(vif_c, ref_vif_c)
    assert np.isclose(vif_const, ref_vif_const)
