#!/usr/bin/env python3

"""
Takes a table associating a key to a list of values, 
and return a table where the new keys are the set of values, and the new values
are the keys.
"""

import argparse
from datetime import datetime
from collections import defaultdict
try:
    from pandas import Series
except ImportError:
    import sys
    print("Could not import Pandas. `inverse_fromserie` will not be usable",
          file=sys.stderr)

START = datetime.now()

SEP='\t'
VALUE_SEP=' '


def inverse(iterator, has_header=True):
    inversed = defaultdict(set)
    
    new_header = list(reversed(next(iterator))) if has_header else None

    for key, value in iterator:
        inversed[value].add(key)

    return inversed, new_header


def iter_one2many_fromfile(infile, key_column=0, value_column=1, sep=SEP,
                        value_sep=VALUE_SEP, has_header=True):
    with open(infile) as inf:
        if has_header:
            header = next(inf).rstrip().split(sep)
            yield (header[key_column], header[value_column])

        for line in inf:
            fields = line.rstrip().split(sep)
            for elem in fields[value_column].split(value_sep):
                yield (fields[key_column], elem)


def inverse_fromfile(infile, key_column=0, value_column=1, sep=SEP,
                     value_sep=VALUE_SEP, has_header=True):
    links_iterator = iter_one2many_fromfile(infile, key_column, value_column, sep,
                                         value_sep, has_header)
    return inverse(links_iterator, has_header)


def inverse_fromserie(serie):
    """Inverse a Pandas Series.
    
    Serie links 'index' -> list of values.
    Output 'value' -> list of indices.
    """
    value_name = serie.name
    index_name = serie.index.name
    iter_serie = ((key, v) for key, values in serie.items() for v in values)
    inversed = Series({val: sorted(keys) for val,keys in
                        inverse(iter_serie, has_header=False)[0].items()},
                      name=index_name)
    inversed.index.name = value_name
    return inversed


def inverse_fromfile_old(infile, key_column=0, value_column=1, sep=SEP,
                      value_sep=VALUE_SEP, has_header=True):
    inversed = defaultdict(set)
    with open(infile) as inf:
        if has_header:
            header = next(inf).rstrip().split(sep)
            new_header = [header[value_column], header[key_column]]
        else:
            new_header = None
        for line in inf:
            fields = line.rstrip().split(sep)
            for elem in fields[value_column].split(value_sep):
                inversed[elem].add(fields[key_column])

    return inversed, new_header


def save(outfile, relation, header=None, sort=True, sep=SEP, value_sep=VALUE_SEP):
    with open(outfile, 'w') as out:
        if header:
            out.write(sep.join(header) + '\n')
        if sort:
            for key in sorted(relation):
                out.write(key + sep + value_sep.join(sorted(relation[key])) + '\n')
        else:
            for key, targets in relation.items():
                out.write(key + sep + value_sep.join(targets) + '\n')


def main(infile, outfile, key_column=0, value_column=1, sort=True, sep=SEP,
         value_sep=VALUE_SEP, has_header=True):
    start = datetime.now()
    print("Started at %s" % (start - START))
    relation, new_header = inverse_fromfile(infile, key_column, value_column,
                                             sep, value_sep, has_header)

    t_elapsed = datetime.now()
    print("Inversion time: %s\nSaving..." % (t_elapsed - start))
    save(outfile, relation, new_header, sort, sep, value_sep)
    print("\rSaving time: %s" % (datetime.now() - t_elapsed))


if __name__=='__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('infile')
    parser.add_argument('outfile')
    parser.add_argument('-k', '--key-column', type=int, default=0)
    parser.add_argument('-v', '--value-column', type=int, default=1)
    parser.add_argument('-S', '--sep', default='\t')
    parser.add_argument('-s', '--value-sep', default=' ')
    parser.add_argument('--no-sort', action='store_false', dest='sort')
    parser.add_argument('--no-header', action='store_false', dest='has_header',
                        help='Do not interpret the first line as header')
    args = parser.parse_args()

    main(**vars(args))

