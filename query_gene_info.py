#!/usr/bin/env python

import urllib2
import fileinput

xmlbase = """<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE Query>
<Query  virtualSchemaName = "default" formatter = "TSV" header = "0" uniqueRows = "0" count = "" datasetConfigVersion = "0.6" >
    <Dataset name = "{dataset}_gene_ensembl" interface = "default" >
        <Filter name = "biotype" value = "protein_coding"/>
        <Attribute name = "ensembl_gene_id" />
        <Attribute name = "percentage_gc_content" />
        <Attribute name = "start_position" />
        <Attribute name = "end_position" />
        <Attribute name = "chromosome_name" />
        <Attribute name = "{species}_homolog_ensembl_gene" />
        <Attribute name = "{species}_homolog_chromosome" />
        <Attribute name = "{species}_homolog_chrom_start" />
        <Attribute name = "{species}_homolog_chrom_end" />
        <Attribute name = "{species}_homolog_orthology_type" />
        <Attribute name = "{species}_homolog_orthology_confidence" />
        <Attribute name = "{species}_homolog_perc_id" />
        <Attribute name = "{species}_homolog_perc_id_r1" />
        <Attribute name = "{species}_homolog_dn" />
        <Attribute name = "{species}_homolog_ds" />
    </Dataset>
</Query>"""

url="http://www.ensembl.org/biomart/martservice?query="
# Ensembl 85
url="http://mar2016.archive.ensembl.org/biomart/martservice?query="
