#!/usr/bin/env python3


import sys
import argparse
import requests


URL = "http://www.ensembl.org/biomart/"
#ARCHIVE = "http://{}.archive.ensembl.org/biomart/"
ARCHIVES = {87: "http://dec2016.archive.ensembl.org/biomart/",
            86: "http://oct2016.archive.ensembl.org/biomart/",
            85: "http://jul2016.archive.ensembl.org/biomart/",
            84: "http://mar2016.archive.ensembl.org/biomart/"}

FORM="martservice?query="

# example query

QUERY = """<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE Query>
<Query  virtualSchemaName = "default" formatter = "{formatter}" header = "{header}" uniqueRows = "{uniqueRows}" count = "{count}" datasetConfigVersion = "0.6" >
    <Dataset name = "{dataset}" interface = "default" >
        {filters}
        {attributes}
    </Dataset>
</Query>"""

FILTER = '<Filter name = "{}" value = "{}"/>'
ATTRIBUTE = '<Attribute name = "{}" />'


def query_fromfile(queryfile, outfile='-', ensembl_version=None):
    with open(queryfile) as qf:
        query = qf.read()
    do_query(query, outfile, ensembl_version)


def query_fromargs(outfile='-', ensembl_version=None, formatter='TSV', header=0,
                    uniqueRows=0, count='', dataset='hsapiens_gene_ensembl',
                    filters=None, attributes=None):
    if filters:
        filter_lines = '\n'.join(FILTER.format(*f.split('=')) for f in filters)
    else:
        filter_lines = ''

    if attributes:
        attr_lines = '\n'.join(ATTRIBUTE.format(a) for a in attributes)
    else:
        attr_lines = ''

    query = QUERY.format(formatter=formatter, header=header, uniqueRows=uniqueRows,
                         count=count, dataset=dataset, filters=filter_lines,
                         attributes=attr_lines)
    
    do_query(query, outfile, ensembl_version=None)


def do_query(query, outfile='-', ensembl_version=None):
    url = ARCHIVES.get(ensembl_version, URL) + FORM + query
    response = requests.get(url)
    
    out = sys.stdout if outfile == '-' else open(outfile, 'w')

    out.write(response.text)

    if outfile == '-': out.close()


def main(**kwargs):
    if kwargs.get('queryfile'):
        for kw in ('formatter', 'header', 'uniqueRows', 'count', 'dataset', 
                   'filter', 'attribute'):
            kwargs.pop(kw)

        query_fromfile(**kwargs)
    
    else:
        kwargs.pop('queryfile')
        query_fromargs(**kwargs)


if __name__=='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('queryfile', nargs='?')
    parser.add_argument('-o', '--outfile', default='-')
    parser.add_argument('-e', '--ensembl-version')
    parser.add_argument('-F', '--formatter', choices=['TSV', 'FASTA'], default='TSV')
    parser.add_argument('-H', '--header', choices=[0, 1], default=0)
    parser.add_argument('-u', '--uniqueRows', choices=[0, 1], default=0)
    parser.add_argument('-c', '--count', default='')
    parser.add_argument('-d', '--dataset', default='hsapiens_gene_ensembl')
    parser.add_argument('-f', '--filters', action='append',
                        help='syntax: "<name>=<value>"')
    parser.add_argument('-a', '--attributes', nargs='+')

    args = parser.parse_args()
    main(**vars(args))

