#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#

"""Print an alignment to stdout, with colors."""

import os.path
import argparse

from Bio import AlignIO

NORMAL       = ""
RESET        = "\033[m"
BOLD         = "\033[1m"
RED          = "\033[31m"
GREEN        = "\033[32m"
YELLOW       = "\033[33m"
BLUE         = "\033[34m"
MAGENTA      = "\033[35m"
CYAN         = "\033[36m"
BOLD_RED     = "\033[1;31m"
BOLD_GREEN   = "\033[1;32m"
BOLD_YELLOW  = "\033[1;33m"
BOLD_BLUE    = "\033[1;34m"
BOLD_MAGENTA = "\033[1;35m"
BOLD_CYAN    = "\033[1;36m"
BG_RED       = "\033[41m"
BG_GREEN     = "\033[42m"
BG_YELLOW    = "\033[43m"
BG_BLUE      = "\033[44m"
BG_MAGENTA   = "\033[45m"
BG_CYAN      = "\033[46m"
BOLD_BG_RED     = "\033[1;41m"
BOLD_BG_GREEN   = "\033[1;42m"
BOLD_BG_YELLOW  = "\033[1;43m"
BOLD_BG_BLUE    = "\033[1;44m"
BOLD_BG_MAGENTA = "\033[1;45m"
BOLD_BG_CYAN    = "\033[1;46m"

COL    = "\033[38;5;%dm"
BG_COL = "\033[48;5;%dm"


nucl2col = {'A': BG_RED,
            'T': BG_BLUE,
            'G': BG_YELLOW,
            'C': BG_GREEN}

ext2fmt = {'.fa':    'fasta',
           '.fasta': 'fasta',
           '.mfa':   'fasta', 
           '.phy':   'phylip'}


def filename2format(filename):
    _, ext = os.path.splitext(filename)
    return ext2fmt[ext]


def makeRGBcolorwheel(levels=5, mix=0):
    wheel = [(levels - i, i, mix) for i in range(levels)] + \
            [(mix, levels - i, i) for i in range(levels)] + \
            [(i, mix, levels - i) for i in range(levels)]
    return wheel


def makeRGBcolorwheel2(levels=5, mix=5):
    wheel = [(mix, levels - i, i) for i in range(levels)] + \
            [(levels - i, i, mix) for i in range(levels)] + \
            [(i, mix, levels - i) for i in range(levels)]
    return wheel

def makeRGBpalette(n=21, offset=0.5):
    wheel = makeRGBcolorwheel()
    step = len(wheel) // n
    first = int(offset * step)
    wheel = wheel[first:] + wheel[:first]
    return [wheel[i*step] for i in range(n)]

def RGB2term(rgb):
    return 16 + 36*rgb[0] + 6*rgb[1] + rgb[2]

def maketermpalette(n=21, offset=0.5):
    return [RGB2term(rgb) for rgb in makeRGBpalette(n, offset)]

def printwheels():
    #for L in range(1, 6):
    L = 5
    #for mix in range(6):
    #    termwheel = [BG_COL % RGB2term(rgb) for rgb in makeRGBcolorwheel(L, mix)]
    #    print(' '.join(termwheel) + ' ')

    termwheel = [BG_COL % RGB2term(rgb) for rgb in makeRGBcolorwheel2(L)]
    print(' '.join(termwheel) + ' ')


def pos2tickmark(pos):
    pass


def makeruler(length, base=1):
    ticks = ['.'] * length
    ticks[0] = str(base)
    for i in range(5, length, 5):
        ticks[i-base] = '|'
    for i in range(10, length, 10):
        count = str(i)
        nchars = len(count)
        for char_i, char in enumerate(count):
            ticks[i-base - (nchars-1-char_i)] = char

    return ''.join(ticks)


def colorizerecord(record):
    return ''.join(nucl2col.get(nucl, RESET)+nucl for nucl in record.seq)

#def printblock(records, namefmt, pad):

def printal(infile, wrap=False):
    ### TODO: wrap to column width
    pad = 4*' '
    unit_delim = '.'
    five_delim = '|'

    with open(infile) as al:
        align = AlignIO.read(al, format=filename2format(infile))

    length = align.get_alignment_length()
    name_len = max(len(record.id) for record in align)

    ruler = makeruler(length)
    namefmt = '%%%ds' % name_len

    try:
        if wrap:
            from subprocess import check_output
            ncols = int(check_output(['tput', 'cols']))
            block_width = ncols - name_len - len(pad)
            assert block_width > 0, \
                "Can't wrap on %d columns because sequence names use %d columns" %\
                (ncols, name_len + pad)
            #print(ncols, name_len)

            for block in range(length // block_width + 1):
                start, stop = (block*block_width, (block+1)*block_width)
                #print(start, stop)
                print(' '*name_len + pad + ruler[start:stop])
                for record in align:

                    print(namefmt % record.id + pad + \
                            colorizerecord(record[start:stop]) + RESET)
                print('')

        else:
            print(' '*name_len + pad + ruler)
            for record in align:
                print(namefmt % record.id + pad + colorizerecord(record) + RESET)
    except BrokenPipeError as err:
        #from os import devnull
        #with open(devnull, 'w') as dn:
        #    print(err, file=dn)
        pass


if __name__ == '__main__':
    #printwheels()

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('infile')
    parser.add_argument('-w', '--wrap', action='store_true', 
                        help='Wrap output to terminal width')
    
    args = parser.parse_args()
    printal(**vars(args))
