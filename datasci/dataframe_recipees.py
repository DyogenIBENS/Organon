#!/usr/bin/env python3

import matplotlib as mpl


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
    c = [mpl.colors.rgb2hex(x) for x in mpl.cm.get_cmap(cmap)(normed)]
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


def tolong(df, x=None, subset=None, varname=None, idname=None):
    """Stack operation while keeping one column as "x".
    
    `x`:      original column(s) to keep in the long dataframe.
    `subset`: columns to be stacked.
    `varname`: name of the new stacked column.
    `idname`:  name of the new column identifying the original column name.

    NOTE: to keep a variable telling the row of origin of the observation, 
    do a `df.reset_index()` first (and add its name to `x`).
    """
    Nx = 1
    if x is not None:
        df = df.set_index(x)
        Nx = 1 if isinstance(x, str) else len(x)
        #else: the original index will become the "level_0" column.
    if subset is not None:
        df = df[subset]
    
    renames = {}
    if varname is not None:
        renames[0] = varname
    if idname is not None:
        renames['level_%d' % Nx] = idname

    return df.stack().reset_index().rename(columns=renames)
