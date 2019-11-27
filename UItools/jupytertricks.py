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


def basic_make_tabs(iterator, make_title=repr, **kwargs):
    if make_title is None:
        make_title = lambda item: None
    tab = widgets.Tab()
    display(tab)
    for i, item in enumerate(iterator):
        out = widgets.Output(**kwargs)
        tab.children += (out,)
        tab.set_title(i, make_title(item))
        with out:
            yield item
            # Do stuff with item in your loop.
            #except BaseException as err:
            #    display(err)
            #    break


class make_tabs(object):
    """Take an iterable, and create a new html tab for displaying the output
    of each iteration.
    
    The instance makes the attributes available:
        - `tab`: ipywidget.Tab instance;
        - `current_out`: ipywidget.Output instance of the current iteration;
        - `current_index`: index of the `current_out` in `tab.children`.
    """
    def __init__(self, iterable, make_title=repr, dynamic=False, **output_kw):
        self.iterable = iterable
        self.make_title = (lambda item: None) if (make_title is None) else make_title
        self.dynamic = dynamic
        self.output_kw = output_kw
        self.tab = widgets.Tab()
        display(self.tab)
        self.current_out = None
        self.current_index = -1

    def __iter__(self):
        for item in self.iterable:
            self.current_out = widgets.Output(**self.output_kw)
            self.tab.children += (self.current_out,)
            self.current_index += 1
            if self.dynamic:
                self.tab.selected_index = self.current_index
            self.tab.set_title(self.current_index, self.make_title(item))
            with self.current_out:
                yield item

