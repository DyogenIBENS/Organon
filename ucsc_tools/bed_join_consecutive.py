#!/usr/bin/env python3

"""Given a tabulated file (BED-like) with columns for (sequence, start, end,
strand), merge all intervals that are overlapping."""

import argparse


def convert_line(line, int_columns=[1,2]):
    return [int(elem) if i in int_columns else elem for i,elem in enumerate(line)]


def load_table(infile, skip_header=True, int_columns=[1,2]):
    with open(infile) as stream:
        if skip_header:
            header = stream.readline().split()
        table = [convert_line(line.split(), int_columns) for line in stream]
    return table


def sort_coords(table, b=1, e=2):
    """by increasing start and decreasing end"""
    pass
    # Not needed here because I want to keep the order as mapped on the mouse.


def join_consecutive_overlapping(processed, name, begin, end, strand,
                                 forward_coords=False):
    """Merge previous coordinates with current, if and only if they overlap.
    forward_coords: whether coordinates on minus strands refer to the forward strand."""
    try:
        prev_n, prev_b, prev_e, prev_s = processed[-1]
    except IndexError:
        processed.append([name, begin, end, strand])
        return

    if (forward_coords and strand == '-') and \
            name == prev_n and \
            strand == prev_s and \
            end >= prev_b:
        processed[-1][1] = min(prev_b, begin)
        processed[-1][2] = max(prev_e, end)

    elif (not forward_coords or strand == '+') and \
            name == prev_n and \
            strand == prev_s and \
            begin <= prev_e:
        processed[-1][1] = min(prev_b, begin)
        processed[-1][2] = max(prev_e, end)
    else:
        processed.append([name, begin, end, strand])


def join_consecutive(processed, name, begin, end, strand,
                     forward_coords=False):
    """Join consecutive fragments of contigs.
    They can be separated by gaps or they can overlap, the only condition is 
    that the following one does not start before the previous one (in order to
    respect the order relative to te reference sequence).
    forward_coords: whether coordinates on minus strands refer to the forward strand."""
    
    try:
        prev_n, prev_b, prev_e, prev_s = processed[-1]
    except IndexError:
        processed.append([name, begin, end, strand])
        return

    if (forward_coords and strand == '-') and \
            name == prev_n and \
            strand == prev_s and \
            end <= prev_e:
        processed[-1][1] = min(prev_b, begin)

    elif (not forward_coords or strand == '+') and \
            name == prev_n and \
            strand == prev_s and \
            begin >= prev_b:
        processed[-1][2] = max(prev_e, end)
    else:
        processed.append([name, begin, end, strand])


def process(table, n=0, b=1, e=2, s=5, join_func=join_consecutive,
            forward_coords=False):
    """given a list of name, start, end, strand; apply the merging function.
    Assumes the table is sorted by name and increasing start."""
    processed = []
    for row in table:
        join_func(processed, name=row[n], begin=row[b], end=row[e], strand=row[s],
                  forward_coords=forward_coords)

    return processed


def save(table, output):
    with open(output, 'w') as out:
        out.write('\n'.join('\t'.join(str(elem) for elem in row) for row in
                            table) + '\n')

def main(infile, output, skip_header=True, only_overlapping=False,
         forward_coords=False, n=0, b=1, e=2, s=5):
    table = load_table(infile, skip_header, int_columns=[b,e])
    join_func = join_consecutive_overlapping if only_overlapping else join_consecutive
    joined = process(table, n, b, e, s, join_func, forward_coords)
    save(joined, output)


if __name__=='__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('infile')
    parser.add_argument('output')
    parser.add_argument('--noskip-header', action='store_false',
                        dest='skip_header')
    parser.add_argument('-n', type=int, default=0,
                        help='column for sequence name (chromosome, scaffold)'\
                             ' [%(default)s]')
    parser.add_argument('-b', type=int, default=1,
                        help='column for beginning [%(default)s]')
    parser.add_argument('-e', type=int, default=2,
                        help='column for end [%(default)s]')
    parser.add_argument('-s', type=int, default=5,
                        help='column for strand [%(default)s]')
    parser.add_argument('-o', '--only-overlapping', action='store_true', 
                        help='for multiple sequences (consecutive in the ' \
                             'input) from the same contig, only join when ' \
                             'they overlap.')
    parser.add_argument('-f', '--forward-coords', action='store_true',
                        help="Coords on the '-' strand are also relative to " \
                             "the forward ('+') strand")

    args = parser.parse_args()

    main(**vars(args))
