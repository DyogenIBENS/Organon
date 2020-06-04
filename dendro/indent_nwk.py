#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from sys import stdin, stdout, stderr
import argparse
import re
from io import StringIO


# The following regex:
# - matches structural tokens: parentheses and comma;
# - includes white spaces before/after (stripped later);
# - matches delimiters to be jumped over: quotes and comments brackets.
RE_STRUCT = re.compile(r'\s*(?:,|\(|\)|\[\&\&NHX:|\'|\")\s*')
PROTECT = {'[&&NHX:': ']', 
           '"': '"',
           "'": "'"}
INDENT = '  '


def go_in(nindent, indent):
    nindent += 1
    # return new indent size, and string to print before the rest of the text
    newstruct = '(\n' + nindent*indent
    return nindent, newstruct

def go_out(nindent, indent):
    nindent -= 1
    newstruct = '\n' + nindent*indent + ')'
    return nindent, newstruct

def go_same(nindent, indent):
    return nindent, ',\n' + nindent*indent

STRUCT_DO = {',': go_same,
             '(': go_in,
             ')': go_out}


def main(inputtree, outputfile, inplace=False, indent=INDENT):
    all_treetxt = inputtree.read()
    all_trees = all_treetxt.rstrip('; \t\n\r').split(';')
    #print(len(all_trees))
    #print(len(all_trees), all_trees)

    if inplace:
        out = StringIO()  # Buffer the output in-memory, in case of error.
    else:
        out = outputfile
    
    try:
        for treetxt in all_trees:
            # Split text by semantic units: nodes, commas, parentheses.
            nindent = 0
            #prev_pos = 0
            #for m in RE_STRUCT.finditer(treetxt):
            structmatch = RE_STRUCT.search(treetxt)
            while structmatch:
                struct = structmatch.group().strip()
                start,end = structmatch.start(), structmatch.end()
                if struct in STRUCT_DO:
                    nindent, newstruct = STRUCT_DO[struct](nindent, indent)

                    out.write(treetxt[:start].rstrip().lstrip() + newstruct)

                    treetxt = treetxt[end:]
                elif struct in PROTECT:
                    #out.write(treetxt[:end])
                    #treetxt = treetxt[structmatch.end():]
                    closing_char = PROTECT[struct]
                    closing_pos = treetxt.find(closing_char, end)
                    out.write(treetxt[:(closing_pos+1)])
                    treetxt = treetxt[(closing_pos+1):]
                else:
                    raise ValueError("Invalid structure character %r" % struct)

                structmatch = RE_STRUCT.search(treetxt)
            out.write(treetxt + ';\n')
        if inplace:
            inputtree.close()
            with open(inputtree.name, 'w') as realout:
                realout.write(out.getvalue())
    
    except BrokenPipeError:  # IOError in Python2.7
        if out is stdout:
            print('indent_nwk:Caught BrokenPipeError to stdout.', file=stderr)
        else:
            raise
    finally:
        if out is not stdout:
            out.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('inputtree', nargs='?', type=argparse.FileType(),
                        default=stdin)
    parser.add_argument('outputfile', nargs='?', type=argparse.FileType('w'),
                        default=stdout)
    parser.add_argument('-i', '--inplace', action='store_true',
                        help='Edit input file in place')
    
    main(**vars(parser.parse_args()))

