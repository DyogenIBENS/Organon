#!/usr/bin/env python3
# -*- coding: utf-8 -*-


"""Read notebook and output a new one with the cell defining given parameters"""


from sys import stdin, stdout, stderr
import argparse as ap
#import json
import re


def parse_extra_args(args):
    kwargs = {}
    k = None
    for a in args:
        if a.startswith('--'):
            k = a[2:]
            if not k:
                raise ValueError('Invalid param ""')
        elif k is None:
            raise ValueError('No param name given for the value %r' % a)
        elif k in kwargs:
            raise ValueError('Parameter %s already given.' % k)
        else:
            kwargs[k] = a
    if k not in kwargs:
        raise ValueError('Trailing value missing for parameter %s' %k)
    return kwargs


def main(infile, outfile, **kwargs):
    # TODO: check valid variable names
    param_reg = re.compile(r'\b(' + r'|'.join(kwargs) + r')\b\s*=\s*(.*?)(\s*(?=#|\\n|"))')
    for i, line in enumerate(infile):
        #line = line.rstrip()
        m = param_reg.search(line)
        if m:
            k, default, trailing = m.groups()
            try:
                line = param_reg.sub(r'%s = %s\3' % (k, kwargs.pop(k)), line)
            except KeyError:
                print('Line %d: "%s"\n\tParameter %r already seen before. Ignored.'
                      %(i, line.rstrip(), k),
                      file=stderr)
        outfile.write(line)
    if kwargs:
        print('Parameters not matched:', ' '.join('%s=%r' % kv for kv in kwargs.items()),
                file=stderr)


if __name__=='__main__':
    parser = ap.ArgumentParser(description=__doc__)
    parser.add_argument('infile', nargs='?', type=ap.FileType('r'), default=stdin)
    parser.add_argument('outfile', nargs='?', type=ap.FileType('w'), default=stdout)
    
    args, extra_args = parser.parse_known_args()
    main(**vars(args), **parse_extra_args(extra_args))

