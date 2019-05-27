#!/usr/bin/env python3
# -*- coding: utf-8 -*-


"""Display comment fields from a NHX tree as a table."""


from sys import stdin, stdout, stderr
import re
import argparse as ap
from dendro.parsers import read_multinewick
import logging
logger = logging.getLogger(__name__)


NHX_open = re.compile(r'\[\&\&NHX:')
NHX_close = re.compile(r'\]')
one_comma = re.compile(r'^\n*,?\n*')


def tab_nhx(treetxt):
    fieldnames = set('')
    fieldvalues = []

    match_open = NHX_open.search(treetxt)
    while match_open:
        start = match_open.end()
        match_close = NHX_close.search(treetxt)
        end = match_close.start()
        comment = treetxt[start:end]
        logger.debug('COMMENT (%d:%d) = %r' % (start, end, comment))

        fields = dict(field.split('=') for field in comment.split(':'))
        fields[''] = treetxt[:match_open.start()]

        fieldnames.update(fields)
        fieldvalues.append(fields)

        treetxt = treetxt[match_close.end():]
        if treetxt.lstrip('\n')[0] == ',':
            fields[''] += ','

        #treetxt = treetxt.lstrip(',\n')  # Risks failing to display a bare topology (unamed nodes -> several commas)
        treetxt = one_comma.sub('', treetxt, count=1)

        logger.debug('NEXT TEXT: %r' % treetxt[:30])
        match_open = NHX_open.search(treetxt)

    if treetxt:
        fieldvalues.append({'': treetxt})

    # Easier to sort fields alphabetically
    fieldnames = sorted(fieldnames)
    yield fieldnames

    for value in fieldvalues:
        yield [value.get(n, '') for n in fieldnames]


def main(infile, outfile):
    for treetxt in read_multinewick(infile, stripchars=''):
        for row in tab_nhx(treetxt):
            outfile.write('\t'.join(row) + '\n')


if __name__ == '__main__':
    logging.basicConfig()
    #logger.setLevel(logging.DEBUG)
    parser = ap.ArgumentParser(description=__doc__)
    parser.add_argument('infile', nargs='?', default=stdin,
                        type=ap.FileType('r'))
    parser.add_argument('outfile', nargs='?', default=stdout,
                        type=ap.FileType('w'))
    
    args = parser.parse_args()
    main(**vars(args))
