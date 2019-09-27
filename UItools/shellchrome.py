#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""ANSI color codes for the terminal"""


from collections import namedtuple


colortuple = namedtuple('colorcode',
                        'BLACK RED GREEN YELLOW BLUE MAGENTA CYAN WHITE')

COLORS_I = colortuple(*range(8))

#These are the sequences needed to get colored ouput
RESET_SEQ = "\033[0m"
COLOR_SEQ = "\033[0;%dm"
BRIGHT_SEQ = "\033[1;%dm"  # Or is bold?
BOLD_SEQ = "\033[1m"

#The background is set with 40 plus the number of the color, and the foreground with 30
BG_COLOR = {'bg'+colname.lower(): COLOR_SEQ % (40+i)
            for colname,i in COLORS_I._asdict().items()}
BG_COLOR.update({'bg'+colname: BRIGHT_SEQ % (30+i)
            for colname,i in COLORS_I._asdict().items()})
COLOR =    {colname.lower(): COLOR_SEQ % (30+i)
            for colname,i in COLORS_I._asdict().items()}
COLOR.update({colname: BRIGHT_SEQ % (30+i)
            for colname,i in COLORS_I._asdict().items()})
