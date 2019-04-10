#!/usr/bin/env python3
# -*- coding: utf-8 -*-


"""Given a list of taxonomic names (from file), return the species tree from
opentreeoflife.org

Output:
    {basename}.tsv (convert name to ott id)
    {basename}.nwk (newick tree)"""


from sys import stdout
import os.path as op
import argparse
import requests
import json
import logging
logger = logging.getLogger(__name__)
logging.basicConfig()


ott_api_url = 'https://api.opentreeoflife.org/v3/'
ott_matchnames = 'tnrs/match_names'
ott_isubtree = 'tree_of_life/induced_subtree'

headers = {'Content-type': 'application/json'}


def get_name_to_ott(taxonlist):
    r0 = requests.post(ott_api_url + ott_matchnames,
                       headers=headers,
                       json={'names': taxonlist})

    if not r0.ok:
        raise RuntimeError("Error %d with 'matchnames' request: %s" % (
                            r0.status_code, r0.reason))
    r0.close()

    names_dict = json.loads(r0.text)

    unmatched_names = names_dict['unmatched_names']
    if unmatched_names:
        logger.warning('Unmatched names : ' + ', '.join(unmatched_names))

    name_to_ott = {}

    for result in names_dict['results']:
        name, matches = result['name'], result['matches']
        goodmatch_i = 0
        if len(matches) > 1:
            logger.warning('%d matches for %s, trying to select the good one, but here are all of them:\n%s',
                    len(matches), name, json.dumps(matches, indent=2))
            for goodmatch_i in range(len(matches)):
                goodmatch = matches[goodmatch_i]
                if goodmatch['is_synonym'] is False \
                    and goodmatch['taxon']['name'] == name \
                    and (goodmatch['taxon']['unique_name'] == name or
                         'in domain Eukaryota' in goodmatch['taxon']['unique_name']):
                    #and goodmatch['taxon']['rank'] in ('species', 'subspecies') \
                    break
            else:
                logger.warning('No good match, keep the first one.')
                goodmatch_i = 0

        name_to_ott[name] = matches[goodmatch_i]['taxon']['ott_id']
    return name_to_ott


def get_induced_subtree(ottidlist):
    r1 = requests.post(ott_api_url + ott_isubtree,
                       headers=headers,
                       json={'label_format': 'name',
                             'ott_ids': ottidlist})

    if not r1.ok:
        raise RuntimeError("Error %d with 'induced_subtree' request: %s" % (
                            r1.status_code, r1.reason))
    r1.close()

    return json.loads(r1.text)['newick']


def main(taxonlistfile, outbase=None):

    if op.exists(outbase + '.tsv'):
        raise FileExistsError(outbase + '.tsv')
    if op.exists(outbase + '.nwk'):
        raise FileExistsError(outbase + '.nwk')

    with open(taxonlistfile) as f:
        taxonlist = [line.rstrip() for line in f if not line.startswith('#')]

    name_to_ott = get_name_to_ott(taxonlist)

    out_names = open(outbase + '.tsv', 'w') if outbase else stdout
    for name in sorted(name_to_ott):
        out_names.write('%s\t%s\n' % (name, name_to_ott[name]))
    if outbase: out_names.close()

    tree = get_induced_subtree(list(name_to_ott.values()))

    out_tree = open(outbase + '.nwk', 'w') if outbase else stdout
    out_tree.write(tree + '\n')
    if outbase: out_tree.close()



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('taxonlistfile')
    parser.add_argument('outbase', nargs='?')
    
    args = parser.parse_args()
    main(**vars(args))
    


