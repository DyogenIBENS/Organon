#!/usr/bin/env python3
# -*- coding: utf-8 -*-


from archiparse import ParseUnit, ParseBloc, floatlist, strlist, OrderedDict, ListOf


Gb_html_parser = ParseBloc('gb.htm',
        [ParseUnit('Nblocks',
            r'Flank positions of the (\d+) selected block\(s\)',
            types=[int]),
         ParseBloc('Flank positions',
             [ParseUnit('rowname', r'^Flanks:'),
              ParseUnit('range', r'(?<!\n)[\t ]+\[(\d+) +(\d+)\]', types=[int, int],
                        repeat=True, optional=True, flags=0)]),
         ParseUnit('filename',
                   r'^New number of positions in (.*):\s*<b>\s*'),
         ParseUnit('positions',
                   r'^(\d+) *</b> *\(([0-9.]+)% of the original (\d+) positions\)[\t ]*$',
                   keys=['new', 'percent', 'old'], types=[int, int, int])
         ])


def parse_gb_html(gbfile):
    with open(gbfile) as s:
        return Gb_html_parser.parse(s.read())[0]
