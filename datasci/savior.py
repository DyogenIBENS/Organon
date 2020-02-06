#!/usr/bin/env python3
# -*- coding: utf-8 -*-


from contextlib import contextmanager
from io import StringIO, BytesIO
from pprint import pprint, pformat
import base64
import matplotlib.pyplot as plt
import pandas as pd
#from pandas.io.formats.style import Styler
from markdown import markdown
from IPython.display import display_html
import logging
logger = logging.getLogger(__name__)


def format_p(txt):
    return '<p>' + txt + '</p>\n'

def format_pre(txt):
    """Display unformatted text (fixed-width)"""
    return '<pre>\n' + txt + '\n</pre>\n'

def format_mkd(txt):
    pass


def format_fig(fig, format='svg', **kwargs):
    with StringIO() as out:
        fig.savefig(out, format=format, **kwargs)

def format_embed_fig(fig, format='png', bbox_inches='tight', **kwargs):
    open_out = StringIO if format=='svg' else BytesIO
    with open_out() as out:
        fig.savefig(out, format=format, bbox_inches=bbox_inches, **kwargs)
        imgdata = out.getvalue()

    # image/svg+xml
    if format == 'svg':
        #TODO: the <?xml> opening tag should be stripped out.
        #          and the '<!DOCTYPE svg PUBLIC' also.
        #          We should use an XML parser for that, but we will assume
        #          That the 3 first lines are to be removed.
        imglines = imgdata.split('\n', 10)  # Should be within the 10 first lines?
        for i, line in imglines:
            if line.startswith('<!-- ') or line.startswith('<svg '):
                break
            assert line.startswith('<?xml ') or line.startswith('<!DOCTYPE ') or line.startswith(' ')

        return '<figure>\n' + '\n'.join(imglines[i:]) + '\n</figure>'
        #<object type="image/svg+xml" data="image.svg"></object>

    imgdata = base64.encodebytes(imgdata).decode()
    return '<figure>\n<img\n src="data:image/%s;base64,%s">\n</figure>\n' % (format, imgdata)
    # alt=
    # style="width:500px;height:600px;" || width="500" height="600"
    # <figcaption>

# See help of IPython.display.display:
#  - `_repr_html_`: return raw HTML as a string, or a tuple (see below).
#  - `_repr_json_`: return a JSONable dict, or a tuple (see below).
#  - `_repr_jpeg_`: return raw JPEG data, or a tuple (see below).
#  - `_repr_png_`: return raw PNG data, or a tuple (see below).
#  - `_repr_svg_`: return raw SVG data as a string, or a tuple (see below).
#  - `_repr_latex_`: return LaTeX commands in a string surrounded by "$",
#                    or a tuple (see below).

from abc import ABC, abstractmethod

class Report(ABC):
    """Base class for HtmlReport and CellReport. Do not instantiate directly."""
    def __enter__(self):
        return self
    
    def __exit__(self, exc_type, exc_value, traceback):
        logger.debug('Exited context: %r %r %r', exc_type, exc_value, traceback)
    
    @abstractmethod
    def write(self, txt):
        """Write raw string"""

    @abstractmethod
    def print(self, *args, sep=' ', end='\n'):
        """Write unformatted text"""

    @abstractmethod
    def pprint(self, *args):
        """pretty print"""

    @abstractmethod
    def show(self, fig=None, **kwargs):
        """Show figure"""
        #return plt.gcf() if fig is None else fig
    
    @abstractmethod
    def html(self, obj):
        """Write the html representation of an object"""
        #return obj.to_html()

    @abstractmethod
    def mkd(self, *strings):
        pass

    def output(self, *args):
        """Automatically dispatch figures, dataframes and text (as printed)."""
        for arg in args:
            if isinstance(arg, (plt.Figure)):  # TODO: seaborn.axisgrid.FacetGrid
                self.show(arg)
            elif isinstance(arg, (pd.DataFrame, pd.io.formats.style.Styler)):
                self.html(arg)
            else:
                self.print(arg)

    def display(self, *args):
        """like output(), but fall back to markdown text."""
        for arg in args:
            if isinstance(arg, (plt.Figure)):  # TODO: seaborn.axisgrid.FacetGrid
                self.show(arg)
            elif isinstance(arg, (pd.DataFrame, pd.io.formats.style.Styler)):
                self.html(arg)
            else:
                self.mkd(arg)


# For default Jupyter Notebook styles,
# see /static/style/style.min.css?v=29c09309dd70e7fe93378815e5f022ae
# or /Users/grant/Sites/jupyter/notebook/notebook/static/notebook/less/renderedhtml.less
class HtmlReport(Report):
    def __init__(self, filename=None, mode='w', figformat='png', title=None, css=None,
                 scripts=None, mathjax=True, style=None, metas=None):
        self.filename = filename
        self.mode = mode
        self.figformat = figformat
        self.printing = False  # Whether a <pre> is already opened
        if css is None: css = []
        if scripts is None: scripts = []
        if metas is None: metas = ['charset="UTF-8"']
        styles = ['''
            .column{0} {{
              float: left;
              width: {0}%;
            }}'''.format(p) for p in range(10,100,10)]
        styles.append('''
            /* Clear floats after the columns */
            .columns:after {
              content: "";
              display: table;
              clear: both;
            }
            pre {
              white-space: pre-wrap;
            }
            body > * {
              max-width:940px;
              margin-right:auto;
              margin-left:auto;
            }
            table .dataframe {
              margin: 0 auto;
              border: none;
              border-collapse: collapse;
              border-spacing: 0;
              color: black;
              table-layout: fixed;
            }
            .dataframe thead {
              border-bottom: 1px solid black;
              vertical-align: bottom;
            }
            .dataframe tr, .dataframe th, .dataframe td {
              text-align: right;
              vertical-align: middle;
              padding: 0.5em 0.5em;
              line-height: normal;
              white-space: normal;
              max-width: none;
              border: none;
            }
            .dataframe tbody tr:nth-child(2n+1) {
              background: #f5f5f5;
            }
            .dataframe tbody tr th:only-of-type {
              vertical-align: middle;
            }
''')
        #TODO: put the notebook CSS code for dataframes (.rendered_html table)
        if style:
            styles.append(style)
        self.head = '\n'.join(
            ['<!DOCTYPE html>', '<head>',
             '<title>%s</title>' if title else ''] +
            ['<meta %s>' % m for m in metas] +
            ['<link rel="stylesheet" type="text/css" href=%s/>' % c for c in css] +
            ['<script type="text/javascript" src=%s/>' % s for s in scripts] +
            ['<script type="text/javascript" async src="https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.5/latest.js?config=TeX-MML-AM_CHTML"></script>' if mathjax else '',
            '<style>\n%s\n</style>' % '\n'.join(styles),
            '</head>',
            '<body>\n'])
        self.closed = True
        self.begun = []  # opened html tags inside body.
    
    def open(self, mode):
        self.handle = StringIO() if self.filename is None else open(self.filename, mode)
        self.closed = False
    
    def close(self):
        logger.debug('Closing </body>, then file handle.')
        self.handle.write('\n</body>\n')
        if self.filename is None:
            display_html(self.handle.getvalue(), raw=True)
            self.handle.close()
        self.closed = True

    def __enter__(self):
        self.open(self.mode)
        if self.mode == 'w':
            self.handle.write(self.head)
        return self
    
    def __exit__(self, exc_type, exc_value, traceback):
        self.close()
        super().__exit__(exc_type, exc_value, traceback)
    
    def write(self, txt):
        """Write raw string to file."""
        self.handle.write(txt)

    def print(self, *args, sep=' ', end='\n'):
        """Write unformatted text (<pre> environment)"""
        if not self.printing:
            self.handle.write('<pre>\n')
            self.printing = True
        self.handle.write(sep.join(str(a) for a in args) + '\n')

    def pprint(self, *args):
        """pretty print (in <pre> environment)"""
        self.print(*(pformat(a) for a in args), sep='\n')

    def show(self, fig=None, **kwargs):
        """Show figure"""
        if self.printing:
            self.handle.write('</pre>\n')
            self.printing = False
        if fig is None: fig = plt.gcf()
        self.handle.write(format_embed_fig(fig, self.figformat, **kwargs))
    
    def html(self, obj):
        if self.printing:
            self.handle.write('</pre>\n')
            self.printing = False
        obj_html = obj._repr_html_()
        if obj_html is None:
            raise ValueError('method _repr_html_() returned None.')
        self.handle.write(obj_html + '\n')
    
    def mkd(self, *strings):
        if self.printing:
            self.handle.write('</pre>\n')
            self.printing = False
        self.handle.write('\n'.join(markdown(s) for s in strings) + '\n')
    
    def begin(self, tag, class_=None, id_=None, attrs=None, print_off=True):
        """Open an html tag"""
        if print_off and self.printing:
            self.handle.write('</pre>\n')
            self.printing = False
        tagstr = tag
        if class_: 
            tagstr += (' class="%s"' % class_)
        if id_: 
            tagstr += (' id="%s"' % id_)
        if attrs:
            tagstr += ' '.join('%s="%s"' for k,v in attrs.items())
        self.handle.write('<%s>\n' % tagstr)
        self.begun.append((tag, class_, id_, attrs))

    def end(self, tag=None, class_=None, id_=None, print_off=True):
        if print_off and self.printing:
            self.handle.write('</pre>\n')
            self.printing = False
        if not self.begun:
            raise ValueError('Trying to close tag but none were opened.')
        for spec_name, spec, prev_spec in zip(('tag', 'class', 'id'),
                                              (tag, class_, id_),
                                              self.begun[-1][:3]):
            if spec is not None and spec != prev_spec:
                raise ValueError('Mismatching closing %s %r (last opened %r)'
                                 % (spec_name, spec, prev_spec))
        self.handle.write('\n</%s>\n' % tag)
        self.begun.pop()

    def ends(self, *tags):
        """Close html tags, if matching opened ones. Otherwise raise ValueError."""
        for tag in tags:
            self.end(tag)

    def end_all(self):
        while self.begun:
            self.handle.write('\n</%s>' % self.begun.pop()[0])

    def columns(self):
        self.begin('div', class_='columns')

    def column(self, p=50):
        self.begin('div', class_='column%d' % p)

    def end_column(self, p=None):
        if p is None:
            begun, begun_class = self.begun[-1][:2]
            if not begun_class.startswith('column') or begun_class == 'columns':
                raise ValueError('Mismatching closing column (last opened "%s .%s")'
                                 % (begun, begun_class))
        else:
            begun_class = 'column%2d' % p
        self.end('div', begun_class)
        self.write('<!-- end %s -->\n' % begun_class)
    
    def end_columns(self):
        self.end('div', 'columns')
        self.write('<!-- end columns -->\n<br/>\n')

    # Haven't found elegant way to use `cols` with or without `with ...` using @contextmanager
    @contextmanager
    def col(self, p=50):
        self.column(p)
        yield
        self.end_column(p)
    
    @contextmanager
    def cols(self):
        self.columns()
        yield
        self.end_columns()


class CellReport(Report):
    #def __init__(self):
    #    pass
    
    def write(self, txt):
        """Write raw string"""
        print(txt, end='')

    def print(self, *args, sep=' ', end='\n'):
        """Write unformatted text"""
        #TODO: turn off pretty printing
        print(*args, sep=sep, end=end)

    def pprint(self, *args):
        """pretty print"""
        for arg in args:
            pprint(arg, sep=sep, end=end)

    def show(self, fig=None, **kwargs):
        """Show figure"""
        if fig is None:
            plt.show(**kwargs)
        else:
            fig.show()  # Not sure this works as expected.

    def html(self, obj):
        display_html(obj.to_html(), raw=True)  # Using the attribute `.to_html()` is necessary for the output method.

    def mkd(self, *strings):
        display_markdown(*strings)


class DoubleReport(HtmlReport, CellReport):
    """Double displaying: in an HTML document, and in the cell output."""
    # __init__ inherited from HtmlReport
    def write(self, txt):
        """Write raw string"""
        HtmlReport.write(self, txt)  # Will call both parents only if HtmlReport uses a super() call...
        CellReport.write(self, txt)

    def print(self, *args, sep=' ', end='\n'):
        """Write unformatted text"""
        HtmlReport.print(self, *args, sep=sep, end=end)
        CellReport.print(self, *args, sep=sep, end=end)

    def pprint(self, *args):
        """pretty print"""
        HtmlReport.pprint(self, *args)
        CellReport.pprint(self, *args)

    def show(self, fig=None, **kwargs):
        """Show figure"""
        HtmlReport.show(self, fig, **kwargs)
        CellReport.show(self, fig, **kwargs)
    
    def html(self, obj):
        HtmlReport.html(self, obj)
        CellReport.html(self, obj)

    def mkd(self, *strings):
        HtmlReport.mkd(self, *strings)
        CellReport.html(self, *strings)


