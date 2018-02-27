#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#

"""Print an alignment to stdout, with colors."""

import sys
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
GREY         = "\033[0;37m"
DGREY        = "\033[0;90m"
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

# 256 color codes
COL    = "\033[38;5;%dm"
BG_COL = "\033[48;5;%dm"

# RGB color codes (format with % (r, g, b))
RGB_ESCAPE = "\033[38;2;%d;%d;%dm"

nucl2col = {'A': BG_RED,
            'T': BG_BLUE,
            'G': BG_YELLOW,
            'C': BG_GREEN,
            'N': GREY,
            '-': DGREY}

# tuples of (bg, fg) codes
CODON_TO_256 = {
    # Stop
    'TAA': (15,16), 'TAG': (15,16), 'TGA': (15,16),
    # Unknown
    #'NNN': (),
    # Methionine
    'ATG': (16,),
    # Phenylalanine
    'TTT': (17,) , 'TTC': (18,),
    # Serine
    'TCT': (46,16), 'TCC': (47,16), 'TCG': (48,16),
    'TCA': (82,16), 'AGT': (83,16), 'AGC': (84,16),
    # Tyrosine
    'TAT': (52,),   'TAC': (88,),
    # Cysteine
    'TGT': (53,),   'TGC': (89,),
    # Tryptophane
    'TGG': (197,),
    # Leucine
    'TTA': (139,16), 'TTG': (140,16), 'CTT': (141,16),
    'CTC': (175,16), 'CTA': (176,16), 'CTG': (177,16),
    # Proline
    'CCT': (24,),   'CCC': (25,), 'CCA': (26,), 'CCG': (27,),
    # Histidine
    'CAT': (58,),   'CAC': (94,),
    # Glutamine
    'CAA': (130,),  'CAG': (166,),
    # Arginine
    'CGT': (38,16), 'CGC': (74,16), 'CGA': (110,16),
    'CGG': (39,16), 'AGA': (75,16), 'AGG': (111,16),
    # Isoleucine
    'ATT': (23,),   'ATC': (59,),   'ATA': (95,),
    # Threonine
    'ACT': (60,),   'ACC': (62,),   'ACA': (62,), 'ACG': (63,),
    # Asparagine
    'AAT': (167,),  'AAC': (203,),
    # Lysine
    'AAA': (134,),  'AAG': (135,),
    # Valine
    'GTT': (142,16),'GTC': (143,16), 'GTA': (144,16), 'GTG': (145,16),
    # Alanine
    'GCT': (179,16),'GCC': (180,16), 'GCA': (215,16), 'GCG': (216,16),
    # Aspartic acid
    'GAT': (214,16),'GAC': (178,16),
    # Glutamic acid
    'GAA': (220,16),'GAG': (221,16),
    # Glycine
    'GGT': (236,),  'GGC': (239,), 'GGA': (242,), 'GGG': (245,)
    }

CODON2COL = {codon: ((BG_COL + COL) % code if len(code)>1 else BG_COL % code) \
                for codon, code in CODON_TO_256.items()}
CODON2COL.update({'---': DGREY})

ext2fmt = {'.fa':    'fasta',
           '.fasta': 'fasta',
           '.mfa':   'fasta',
           '':       'fasta',
           '.phy':   'phylip-relaxed'}


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


def makeruler(length, base=1, stepwidth=1):
    """Set stepwidth=3 for codons"""
    nsteps = length // stepwidth
    minortick='.'
    majortick='|'
    ticks = list(minortick + ' '*(stepwidth-1)) * nsteps
    ticks[0] = str(base)
    for i in range(5, nsteps, 5):
        ticks[(i-base)*stepwidth] = majortick
    for i in range(10, nsteps, 10):
        # update the character at the tick, by taking into account the length
        # of the number.
        count = str(i)
        nchars = len(count)
        for char_i, char in enumerate(count):
            ticks[(i-base)*stepwidth - (nchars-1-char_i)] = char

    return ''.join(ticks)


def colorizerecord(record):
    return ''.join(nucl2col.get(nucl, '')+nucl+RESET for nucl in record.seq)


def iter_codons(seq):
    Nnucl = len(seq)
    assert Nnucl % 3 == 0
    #N = Nnucl // 3
    for i in range(0, Nnucl, 3):
        yield str(seq[i:(i+3)])
    

def codoncolorizerecord(record):
    colorized=''
    unknown_codons = set()
    for codon in iter_codons(record.seq):
        try:
            codoncol = CODON2COL[codon]
        except KeyError:
            unknown_codons.add(codon)
            codoncol = RED
        colorized += codoncol + codon + RESET

    if unknown_codons:
        print("WARNING: unknown codons: %s" % ' '.join(unknown_codons),
                file=sys.stderr)

    return colorized

#def printblock(records, namefmt, pad):

def printal(infile, wrap=False, format=None, slice=None, codon=False):
    ### TODO: wrap to column width
    pad = 4*' '
    #unit_delim = '.'
    #five_delim = '|'

    #with open(infile) as al:
    align = AlignIO.read(infile, format=(format or filename2format(infile.name)))

    length = align.get_alignment_length()
    name_len = max(len(record.id) for record in align)

    if codon:
        stepwidth = 3
        colorize = codoncolorizerecord
    else:
        stepwidth = 1
        colorize = colorizerecord

    ruler = makeruler(length, stepwidth=stepwidth)
    
    namefmt = '%%%ds' % name_len

    try:
        if wrap:
            from subprocess import check_output
            ncols = int(check_output(['tput', 'cols']))
            block_width = ncols - name_len - len(pad)
            if codon:
                block_width -= (block_width % 3)

            assert block_width>0, \
                "Can't wrap on %d columns because sequence names use %d columns" %\
                (ncols, name_len + pad)
            #print(ncols, name_len)

            if slice:
                # -1 because coords are taken in base 1
                if codon:
                    slstart, slend = [(int(pos)-1)*3 for pos in slice.split(':')]
                else:
                    slstart, slend = [int(pos)-1 for pos in slice.split(':')]

                length = slend - slstart
            else:
                slstart, slend = 0, length

            nblocks = length // block_width + 1
            for block in range(nblocks):
                start, stop = (block*block_width, (block+1)*block_width)
                start += slstart
                stop = min(stop + slstart, slend)

                #print(start, stop)
                print(' '*name_len + pad + ruler[start:stop])
                for record in align:

                    print(namefmt % record.id + pad + \
                            colorize(record[start:stop]) + RESET)
                if block < nblocks-1:
                    print('')

        else:
            print(' '*name_len + pad + ruler)
            for record in align:
                print(namefmt % record.id + pad + colorize(record) + RESET)
    except BrokenPipeError as err:
        #from os import devnull
        #with open(devnull, 'w') as dn:
        #    print(err, file=dn)
        pass


if __name__ == '__main__':
    #printwheels()

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('infile', nargs='?', default=sys.stdin,
                        type=argparse.FileType('r'))
    parser.add_argument('-L', '--nowrap', action='store_false', dest='wrap',
                        help='Do not wrap output to terminal width')
    #parser.add_argument('-w', '--wrap', action='store_true', 
    #                    help='Wrap output to terminal width')
    parser.add_argument('-f', '--format', help='Force format usage.' \
                        ' Can be any format accepted by Bio.alignIO')
    parser.add_argument('-s', '--slice',
                        help='select positions (start:end). 1-based, end excluded')
    parser.add_argument('-c', '--codon', action='store_true', 
                        help='Colorize and index alignment by codons.')
    
    args = parser.parse_args()
    printal(**vars(args))
