#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from sys import stdin, stdout
#import fileinput
import argparse
import re
#import fileinput

RE_SPACE = re.compile(r'\s+')
RE_STRUCT = re.compile(r',|\(|\)|\[\&\&NHX:|\'|\"') # Exclude contents between quotes and in comments
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

def main(inputtree, outputfile, indent=INDENT):
    #if inputtree is not '-':
    #with open(inputtree) as intree:
    all_treetxt = RE_SPACE.sub('', inputtree.read())
    #else:
    #    all_treetxt = RE_SPACE.sub('', ''.join(line.rstrip() for line in stdin.readlines()))
    all_trees = all_treetxt.rstrip('; \t\n\r').split(';')
    #print(len(all_trees))
    #print(len(all_trees), all_trees)

    out = outputfile
    
    try:
        for treetxt in all_trees:
            # Split text by semantic units: nodes, commas, parentheses.
            nindent = 0
            #prev_pos = 0
            #for m in RE_STRUCT.finditer(treetxt):
            structmatch = RE_STRUCT.search(treetxt)
            while structmatch:
                struct = structmatch.group()
                start,end = structmatch.start(), structmatch.end()
                if struct in STRUCT_DO:
                    nindent, newstruct = STRUCT_DO[struct](nindent, indent)

                    out.write(treetxt[:start] + newstruct)

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
    except:
        if out.name != '<stdout>':
            out.close()
        raise

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('inputtree', nargs='?', type=argparse.FileType(),
                        default=stdin)
    parser.add_argument('outputfile', nargs='?', type=argparse.FileType('w'),
                        default=stdout)
    
    main(**vars(parser.parse_args()))

