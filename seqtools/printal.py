#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Copyright © Guillaume LOUVEL
# e-mail : guillaume.louvel@ens.fr
# License LGPL v3

"""Print an alignment to stdout, with colors."""

from sys import stdin
import os.path
import re
import argparse
import logging
logging.basicConfig(format='%(levelname)s:%(funcName)s:%(message)s')


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
            'U': BG_CYAN,
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

# tuples of (bg, fg) codes
AA_TO_256 = {
    '*': (15,16),  # Stop
    #'X': (),      # Unknown
    'M': (16,),    # Methionine
    'F': (17,),    # Phenylalanine
    'S': (46,16),  # Serine
    'Y': (52,),    # Tyrosine
    'C': (53,),    # Cysteine
    'W': (197,),   # Tryptophane
    'L': (139,16), # Leucine
    'P': (24,),    # Proline
    'H': (58,),    # Histidine
    'Q': (130,),   # Glutamine
    'R': (38,16),  # Arginine
    'I': (23,),    # Isoleucine
    'T': (60,),    # Threonine
    'N': (167,),   # Asparagine
    'K': (134,),   # Lysine
    'V': (142,16), # Valine
    'A': (179,16), # Alanine
    'D': (214,16), # Aspartic acid
    'E': (220,16), # Glutamic acid
    'G': (236,)    # Glycine
    #'B': , # Aspartic acid/Asparagine
    #'Z':   # Glutamic acid/Glutamine
    }

AA3_TO_256 = {
    'TERM': (15,16), # Stop
    #'Xaa': (),      # Unknown
    'Met': (16,),    # Methionine
    'Phe': (17,),    # Phenylalanine
    'Ser': (46,16),  # Serine
    'Tyr': (52,),    # Tyrosine
    'Cys': (53,),    # Cysteine
    'Trp': (197,),   # Tryptophane
    'Leu': (139,16), # Leucine
    'Pro': (24,),    # Proline
    'His': (58,),    # Histidine
    'Gln': (130,),   # Glutamine
    'Arg': (38,16),  # Arginine
    'Ile': (23,),    # Isoleucine
    'Thr': (60,),    # Threonine
    'Asn': (167,),   # Asparagine
    'Lys': (134,),   # Lysine
    'Val': (142,16), # Valine
    'Ala': (179,16), # Alanine
    'Asp': (214,16), # Aspartic acid
    'Glu': (220,16), # Glutamic acid
    'Gly': (236,)    # Glycine
    #'Asx': , # Aspartic acid/Asparagine
    #'Glx':   # Glutamic acid/Glutamine
    }

CODON2COL = {codon: ((BG_COL + COL) % code if len(code)>1 else BG_COL % code) \
                for codon, code in CODON_TO_256.items()}
CODON2COL.update({'---': DGREY})

AA2COL = {aa: ((BG_COL + COL) % code if len(code)>1 else BG_COL % code) \
                for aa, code in AA_TO_256.items()}
AA2COL.update({'-': DGREY})

AA32COL = {aa3: ((BG_COL + COL) % code if len(code)>1 else BG_COL % code) \
                for aa3, code in AA3_TO_256.items()}
AA32COL.update({'---': DGREY})

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


def fastcolorizerecord(record, residu2col=nucl2col):
    return ''.join(residu2col.get(r.upper(), '')+r+RESET for r in record.seq)


def iter_step(seq, stepwidth=3):
    Nnucl = len(seq)
    assert Nnucl % stepwidth == 0
    #N = Nnucl // 3
    for i in range(0, Nnucl, stepwidth):
        yield str(seq[i:(i+stepwidth)])
    

def colorizerecord(record, residu2col=CODON2COL, stepwidth=3):
    colorized=''
    unknown_residus = set()
    for residu in iter_step(record.seq, stepwidth):
        try:
            residucol = residu2col[residu.upper()]
        except KeyError:
            unknown_residus.add(residu)
            residucol = RED
        colorized += residucol + residu + RESET

    if unknown_residus:
        logging.warning("Unknown codons: %s", ' '.join(unknown_residus))

    return colorized

#def printblock(records, namefmt, pad):

def printal(infile, wrap=False, format=None, slice=None, alphabet='codon',
            start0=False):
    padlen = 4
    pad = padlen*' '
    #unit_delim = '.'
    #five_delim = '|'
    ruler_end_reg = re.compile(r'\d+$')
    ruler_start_reg = re.compile(r'\d+')

    #with open(infile) as al:
    align = AlignIO.read(infile, format=(format or filename2format(infile.name)))

    length = align.get_alignment_length()
    name_len = max(len(record.id) for record in align)

    if alphabet == 'codon':
        stepwidth = 3
        colorize = colorizerecord
        residu2col = CODON2COL
    elif alphabet == 'aa3':
        stepwidth = 3
        colorize = colorizerecord
        residu2col = AA32COL
    elif alphabet == 'nucl':
        stepwidth = 1
        colorize = fastcolorizerecord
        residu2col = nucl2col
    elif alphabet == 'aa':
        stepwidth = 1
        colorize = fastcolorizerecord
        residu2col = AA2COL


    start1 = int(not start0)
    ruler = makeruler(length, base=start1, stepwidth=stepwidth)
    
    namefmt = '%%%ds' % name_len

    if slice:
        # -1 because coords are taken in base 1
        sliceparts = slice.split(':')
        if not sliceparts[0]: sliceparts[0] = 0
        if not sliceparts[1]: sliceparts[1] = length

        slstart, slend = [(int(pos)-start1)*stepwidth for pos in sliceparts]

        length = slend - slstart
    else:
        slstart, slend = 0, length

    try:
        if wrap:
            from subprocess import check_output
            ncols = int(check_output(['tput', 'cols']))
            block_width = ncols - name_len - padlen
            if stepwidth > 1:
                block_width -= (block_width % 3)

            assert block_width>0, \
                "Can't wrap on %d columns because sequence names use %d columns" %\
                (ncols, name_len + pad)
            #print(ncols, name_len)
        else:
            ncols = length + 1 + name_len + padlen
            block_width = length + 1

        nblocks = length // block_width + 1
        prev_block_number = ''
        for block in range(nblocks):
            start, stop = (block*block_width, (block+1)*block_width)
            start += slstart
            stop = min(stop + slstart, slend)

            blockruler = ruler[start:stop]

            # If the end of the previous column number was split, add it here
            rulerline = ' '*(name_len + padlen - len(prev_block_number)) \
                        + prev_block_number + blockruler

            end_match = ruler_end_reg.search(blockruler)
            continue_match = ruler_start_reg.match(ruler[stop:])
            if end_match and continue_match:
                prev_block_number = end_match.group()
                rulerline = rulerline.rstrip('0123456789')
            else:
                prev_block_number = ''

            print(rulerline)
            #print(blockruler.rjust(endcol))
            
            for record in align:
                print(namefmt % record.id + pad + \
                        colorize(record[start:stop], residu2col) + RESET)
            if block < nblocks-1:
                print('')


        #else:
        #    print(' '*name_len + pad + ruler[slstart:slend])
        #    for record in align[:,slstart:slend]:
        #        print(namefmt % record.id + pad + colorize(record) + RESET)
    
    except BrokenPipeError as err:
        #from os import devnull
        #with open(devnull, 'w') as dn:
        #    print(err, file=dn)
        pass


if __name__ == '__main__':
    #printwheels()

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('infile', nargs='?', default=stdin,
                        type=argparse.FileType('r'))
    parser.add_argument('-L', '--nowrap', action='store_false', dest='wrap',
                        help='Do not wrap output to terminal width')
    #parser.add_argument('-w', '--wrap', action='store_true', 
    #                    help='Wrap output to terminal width')
    parser.add_argument('-f', '--format', help='Force format usage.' \
                        ' Can be any format accepted by Bio.alignIO')
    parser.add_argument('-s', '--slice',
                        help='select positions (start:end). 1-based, end excluded')
    parser.add_argument('-c', '--codon', action='store_const', dest='alphabet',
                        const='codon', default='nucl',
                        help='Colorize and index alignment by codons.')
    parser.add_argument('-a', '--aa', action='store_const', dest='alphabet',
                        const='aa',
                        help='Colorize and index alignment by amino-acids.')
    parser.add_argument('--aa3', action='store_const', dest='alphabet',
                        const='aa3',
                        help='Colorize and index alignment by 3-letters amino-acids.')
    parser.add_argument('-0', '--start0', action='store_true',
                        help='Use 0-based coordinates.')
    
    args = parser.parse_args()
    printal(**vars(args))
