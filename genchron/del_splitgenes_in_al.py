#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import os
import os.path as op
import argparse
import logging
logger = logging.getLogger(__name__)
ch = logging.StreamHandler()
ch.setFormatter(logging.Formatter("%(levelname)s:%(message)s"))
logger.addHandler(ch)


from del_splitgenes_in_subtrees import ANCESTORLIST, \
                                        load_cladeof, \
                                        find_src_files, \
                                        iter_splitgenes_ancgenes

from seqtools.al_merge_splits import main_fromlist as main_merge_splits


def main(SGlistfile, alignments_dir, src_subtreedir='subtreesCleanO2',
            out_subtreedir='subtreesCleanO2mergedSG'):

    cladeof = load_cladeof()

    outputted = []
    count_genesplits = 0
    
    for split_ancgene, split_descendants in iter_splitgenes_ancgenes(SGlistfile):
        src_file = find_src_files(split_ancgene, cladeof, alignments_dir, src_subtreedir,
                                  end='*_genes.fa')
        out_file = src_file.replace(src_subtreedir, out_subtreedir)
        out_dir = op.dirname(out_file)

        if out_file in outputted:
            src_file = out_file
        
        if not op.exists(out_dir):
            os.mkdir(out_dir)
        
        try:
            main_merge_splits(src_file, [split_descendants], out_file)
        
            outputted.append(out_file)
            count_genesplits += 1
        except KeyError as err:
            logger.warning('KeyError in %s: %s', src_file, err.args[0])

    print('Found %d genesplits.' % count_genesplits, file=stderr)
    print('\n'.join(outputted))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('SGlistfile')
    #parser.add_argument('genetreelistfile')
    parser.add_argument('alignments_dir', nargs='?', default='.')
    parser.add_argument('--out_subtreedir', default='subtreesCleanO2mergedSG')

    args = parser.parse_args()
    main(**vars(args))
