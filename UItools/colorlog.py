# -*- coding: utf-8 -*-


"""
Add ANSI color to the levelname, for nicer printing of log records.

Taken from: https://stackoverflow.com/a/384125
"""


import logging
from copy import copy


BLACK, RED, GREEN, YELLOW, BLUE, MAGENTA, CYAN, WHITE = range(8)

#The background is set with 40 plus the number of the color, and the foreground with 30

#These are the sequences need to get colored ouput
RESET_SEQ = "\033[0m"
COLOR_SEQ = "\033[1;%dm"
BOLD_SEQ = "\033[1m"


COLORS = {
    'WARNING': YELLOW,
    'INFO': WHITE,
    'DEBUG': BLUE,
    'CRITICAL': YELLOW,
    'ERROR': RED
}


class ColoredFormatter(logging.Formatter):
    def __init__(self, fmt=None, datefmt=None, style='%', bg=True):
        logging.Formatter.__init__(self, fmt, datefmt, style)
        self.base = 40 if bg else 30

    def format(self, record):
        levelname = record.levelname
        if levelname in COLORS:
            record = copy(record)
            record.levelname = (COLOR_SEQ % (self.base + COLORS[levelname])
                               + levelname + RESET_SEQ)
        return logging.Formatter.format(self, record)

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
