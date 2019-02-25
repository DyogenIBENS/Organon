#!/usr/bin/env python3
# -*- coding: utf-8 -*-


### Test colorlog


import sys
import logging
from UItools import colorlog


print('Checking display of colors: what do you see:')

for i, colname in zip(range(8),
                      ('black', 'red', 'green', 'yellow', 'blue', 'magenta',
                       'cyan', 'white')):

    print(colorlog.COLOR_SEQ % (30 + i) + 'fg:' + colname + '?' + colorlog.RESET_SEQ, end='\t')
    print(colorlog.COLOR_SEQ % (40 + i) + 'bg:' + colname + '?' + colorlog.RESET_SEQ)


logger = logging.getLogger('UItools.colorlog._tests')
logger.setLevel(logging.DEBUG)

#coloredfmt = '$BGLVL%(levelname)s$RESET:$BLACK$BOLD%(name)s$RESET:$LVL%(message)s$RESET'
coloredfmt = colorlog.BASIC_FORMAT

fmtter = colorlog.ColoredFormatter(coloredfmt + ' (Should be colored).')
sh = logging.StreamHandler(sys.stdout)
sh.setFormatter(fmtter)
fhc = logging.FileHandler('colorlog_test_withcolor.txt', 'w')
fhc.setFormatter(fmtter)
fh = logging.FileHandler('colorlog_test_nocolor.txt', 'w')
fh.setFormatter(logging.Formatter(logging.BASIC_FORMAT))

logger.addHandler(sh)
logger.addHandler(fhc)
logger.addHandler(fh)

#logging.basicConfig()


for level in ('DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'):
    emit = getattr(logger, level.lower())
    emit('Level is: %s', level)

# Uninstall all handlers
logger.handlers = logger.handlers[:-3]
