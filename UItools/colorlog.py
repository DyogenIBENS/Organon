# -*- coding: utf-8 -*-


"""
Add ANSI color to the levelname, for nicer printing of log records.

Taken from: https://stackoverflow.com/a/384125
"""


import string
from collections import namedtuple
import logging
from copy import copy

colortuple = namedtuple('colorcode',
                        'BLACK RED GREEN YELLOW BLUE MAGENTA CYAN WHITE')


#The background is set with 40 plus the number of the color, and the foreground with 30

#These are the sequences needed to get colored ouput
RESET_SEQ = "\033[0m"
COLOR_SEQ = "\033[0;%dm"
BRIGHT_SEQ = "\033[1;%dm"  # Or is bold?
BOLD_SEQ = "\033[1m"


COLORS_I = colortuple(*range(8))

LVL_I = {
    'WARNING':  COLORS_I.YELLOW,
    'INFO':     COLORS_I.WHITE,
    'DEBUG':    COLORS_I.BLUE,
    'CRITICAL': COLORS_I.YELLOW,
    'ERROR':    COLORS_I.RED}


BG_COLOR = {'bg'+colname.lower(): COLOR_SEQ % (40+i)
            for colname,i in COLORS_I._asdict().items()}
BG_COLOR.update({'bg'+colname: BRIGHT_SEQ % (30+i)
            for colname,i in COLORS_I._asdict().items()})
COLOR =    {colname.lower(): COLOR_SEQ % (30+i)
            for colname,i in COLORS_I._asdict().items()}
COLOR.update({colname: BRIGHT_SEQ % (30+i)
            for colname,i in COLORS_I._asdict().items()})
BG_LVLCOLOR = {lvl: COLOR_SEQ % (40+i) for lvl,i in LVL_I.items()}
LVLCOLOR    = {lvl: COLOR_SEQ % (30+i) for lvl,i in LVL_I.items()}

LVLCOLOR.update(ERROR=COLOR['RED'], CRITICAL=BG_COLOR['bgred'])

BASIC_FORMAT = '$LVL%(levelname)s$RESET:$BLACK%(name)s$RESET:%(message)s'


class ColoredFormatter(logging.Formatter):

    def __init__(self, fmt=None, datefmt=None, style='%'):
        if fmt is not None:
            colored_template = string.Template(fmt)
            fmt = colored_template.substitute(
                            BGLVL=('%(bglvlcol)s' if style=='%' else '{bglvlcol}'),
                            LVL=('%(lvlcol)s' if style=='%' else '{lvlcol}'),
                            RESET=RESET_SEQ,
                            BOLD=BOLD_SEQ,
                            **COLOR, **BG_COLOR)
        super(ColoredFormatter, self).__init__(fmt, datefmt, style)


    def format(self, record):
        levelname = record.levelname
        if levelname in LVLCOLOR:
            record = copy(record)
            record.__dict__.update(bglvlcol=BG_LVLCOLOR[levelname],
                                   lvlcol=LVLCOLOR[levelname])
        return super(ColoredFormatter, self).format(record)

# Pb: if there were width specifications, they aren't visually respected anymore (because ANSI escape codes consume width internally, but not visibly)...
# Better: using newstyle string formatting, then:
#import string
#formatparse = string.Formatter().parse
#
#class newColoredFormatter(logging.Formatter):
#    def __init__(self, fmt=None, datefmt=None, bg=True):
#        self.base = 40 if bg else 30
#        self.colors = {level: code+self.base for level,code in COLORS.items()}
#        
#        for txt, fieldname, format_spec, converter in formatparse(fmt):
#            if fieldname == 'levelname':
#                break
#        else:
#            logging.Formatter.__init__(self, fmt, datefmt, style='{')
#            return
#
#        levelname_fmtstr = '{' + field_name
#        if converter: levelname_fmtstr += '!' + converter
#        if format_spec: levelname_fmtstr += ':' + format_spec
#        levelname_fmtstr += '}'
#
#        levelname_fmtstr = .. + levelname_fmtstr + ..
#        fmt = 
#        logging.Formatter.__init__(self, fmt, datefmt, style='{')
#
#    def format(self, record):
#        levelname = record.levelname
#        if levelname in COLORS:
#            record = copy(record)
#            record.levelname = (COLOR_SEQ % (self.base + COLORS[levelname])
#                               + levelname + RESET_SEQ)
#        return logging.Formatter.format(self, record)
