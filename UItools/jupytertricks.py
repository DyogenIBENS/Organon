#!/usr/bin/env python3
# -*- coding: utf-8 -*-


## Capture cell output (stored in variable)

# %%capture captured --no-stderr
#captured.show()
#print(captured.stdout)

## Or with a context manager:
from IPython.utils.io import capture_output
import ipywidgets as widgets

#with capture_output() as captured:
#    MyFunction()

# Create tabs for multiple outputs
#
#out1 = widgets.Output()
#out2 = widgets.Output()
#
#tab = widgets.Tab(children = [out1, out2])
#tab.set_title(0, 'First')
#tab.set_title(1, 'Second')
#display(tab)
#
#with out1:
#    fig1, axes1 = plt.subplots()
#    data1.hist(ax = axes1)
#    plt.show(fig1)
#
#with out2:
#    fig2, axes2 = plt.subplots()
#    data2.hist(ax = axes2)
#    plt.show(fig2)


def make_tabs(iterator, **kwargs):
    outputs = []
    for item in iterator:
        out = widgets.Output(**kwargs)
        outputs.append(out)
        with out:
            try:
                yield item
                # Do stuff with item in your loop.
            except BaseException as err:
                display(err)
                break
    tab = widgets.Tab(children=outputs)
    display(tab)

