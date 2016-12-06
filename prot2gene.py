#!/usr/bin/env python

"""Convert any Ensembl protein ID to gene ID and vice-versa
EXAMPLE:
    ./prot2gene ~/ws2/DUPLI_data85/gene_info/%s_gene_info.tsv <fastafiles>
"""

import re
import sys
import os.path
import argparse
from bz2 import BZ2File

def myopen(filename, *args, **kwargs):
    if filename.endswith('.bz2'):
        return BZ2File(filename, *args, **kwargs)
    else:
        return open(filename, *args, **kwargs)


def convert_prot2species(modernID):
    prot2sp = { #'Y': 'Saccharomyces cerevisiae',  # there are C.elegans prot with 'Y' too
                'Q0': 'Saccharomyces cerevisiae',
                'FB': 'Drosophila melanogaster',
                #'WBGene0': 'Caenorhabditis elegans',  # No consensus
                'ENSCINP': 'Ciona intestinalis',
                'ENSCSAV': 'Ciona savignyi',
                'ENSCSAP': 'Chlorocebus sabaeus',
                'ENSPMAP': 'Petromyzon marinus',
                'ENSXETP': 'Xenopus tropicalis',
                'ENSPSIP': 'Pelodiscus sinensis',
                'ENSGALP': 'Gallus gallus',
                'ENSMGAP': 'Meleagris gallopavo',
                'ENSTGUP': 'Taeniopygia guttata',
                'ENSFALP': 'Ficedula albicollis',
                'ENSAPLP': 'Anas platyrhynchos',
                'ENSACAP': 'Anolis carolinensis',
                'ENSOANP': 'Ornithorhynchus anatinus',
                'ENSMEUP': 'Macropus eugenii',
                'ENSSHAP': 'Sarcophilus harrisii',
                'ENSMODP': 'Monodelphis domestica',
                'ENSLAFP': 'Loxodonta africana',
                'ENSETEP': 'Echinops telfairi',
                'ENSPCAP': 'Procavia capensis',
                'ENSDNOP': 'Dasypus novemcinctus',
                'ENSCHOP': 'Choloepus hoffmanni',
                'ENSSARP': 'Sorex araneus',
                'ENSEEUP': 'Erinaceus europaeus',
                'ENSMLUP': 'Myotis lucifugus',
                'ENSPVAP': 'Pteropus vampyrus',
                'ENSTTRP': 'Tursiops truncatus',
                'ENSBTAP': 'Bos taurus',
                'ENSOARP': 'Ovis aries',
                'ENSVPAP': 'Vicugna pacos',
                'ENSSSCP': 'Sus scrofa',
                'ENSECAP': 'Equus caballus',
                'ENSMPUP': 'Mustela putorius furo',
                'ENSAMEP': 'Ailuropoda melanoleuca',
                'ENSCAFP': 'Canis lupus familiaris',
                'ENSFCAP': 'Felis catus',
                'ENSTBEP': 'Tupaia belangeri',
                'ENSPANP': 'Papio anubis',
                'ENSMMUP': 'Macaca mulatta',
                'ENSPPYP': 'Pongo abelii',
                'ENSGGOP': 'Gorilla gorilla gorilla',
                'ENSPTRP': 'Pan troglodytes',
                'ENSP000':    'Homo sapiens',       # ENSG
                'ENSNLEP': 'Nomascus leucogenys',
                'ENSCJAP': 'Callithrix jacchus',
                'ENSTSYP': 'Tarsius syrichta',
                'ENSOGAP': 'Otolemur garnettii',
                'ENSMICP': 'Microcebus murinus',
                'ENSOPRP': 'Ochotona princeps',
                'ENSOCUP': 'Oryctolagus cuniculus',
                'ENSCPOP': 'Cavia porcellus',
                'ENSRNOP': 'Rattus norvegicus',
                'ENSMUSP': 'Mus musculus',
                'ENSSTOP': 'Ictidomys tridecemlineatus',
                'ENSDORP': 'Dipodomys ordii',
                'ENSLACP': 'Latimeria chalumnae',
                'ENSLOCP': 'Lepisosteus oculatus',
                'ENSGACP': 'Gasterosteus aculeatus',
                'ENSTNIP': 'Tetraodon nigroviridis',
                'ENSTRUP': 'Takifugu rubripes',
                'ENSONIP': 'Oreochromis niloticus',
                'ENSORLP': 'Oryzias latipes',
                'ENSPFOP': 'Poecilia formosa',
                'ENSXMAP': 'Xiphophorus maculatus',
                'ENSGMOP': 'Gadus morhua',
                'ENSAMXP': 'Astyanax mexicanus',
                'ENSDARP': 'Danio rerio'}
    try:
        return prot2sp[modernID[:7]]
    except KeyError:
        try:
            # Saccharomyces cerevisiae (Q0) or Drosophila melanogaster
            return prot2sp[modernID[:2]]
        except KeyError:
            if re.match('Y[A-Z]', modernID):
                return 'Saccharomyces cerevisiae'
            else:
                return 'Caenorhabditis elegans'


def grep_prot(filename, protID, cprot=2, cgene=0):
    with myopen(filename) as IN:
        for line in IN:
            fields = line.rstrip('\r\n').split('\t')
            if fields[cprot] == protID:
                return fields[cgene]


def grep_gene(filename, geneID, cprot=2, cgene=0):
    with myopen(filename) as IN:
        for line in IN:
            fields = line.rstrip('\r\n').split('\t')
            if fields[cgene] == geneID:
                return fields[cprot]


def convert_prot2gene(protID):
    sp = convert_prot2species(protID)
    if shorten_species:
        spsplit = sp.split()
        sp2 = spsplit[0][0].lower() + spsplit[-1]
    else:
        sp2 = sp.replace(' ', '.')
    return grep_prot(gene_info % sp2, protID)


def rewrite_fastafile(fastafile, outputformat="{0}_genes.fa", cprot=2, cgene=0):
    genetree, ext = os.path.splitext(fastafile)
    if ext == '.bz2': genetree, ext = os.path.splitext(genetree)
    genetreedir, genetreefile = os.path.split(genetree)
    #print >>sys.stderr, genetree, genetreedir, genetreefile
    outfile = outputformat.format(genetreefile)
    if os.path.exists(outfile):
        if force_overwrite:
            print >>sys.stderr, "(Overwriting %s)" % outfile
        else:
            print >>sys.stderr, "%s exists. Skipping." % outfile
            return

    # avoid duplicate genes
    found = {}
    unknowns = 0
    with myopen(fastafile) as IN, myopen(outfile, 'w') as OUT:
        for line in IN:
            if line[0] == '>':
                protID = line[1:].split('/')[0]
                geneID = convert_prot2gene(protID)
                if not geneID and protID.startswith('ENSCSAP'):
                    protID = protID.replace('ENSCSAP', 'ENSCSAVP')
                    geneID = convert_prot2gene(protID)
                    print >>sys.stderr, "converting", geneID
                    if geneID:
                        # Fit names in tree
                        geneID = geneID.replace('ENSCSAVG', 'ENSCSAG')
                if not geneID:
                    unknowns += 1
                    geneID = "unknown_gene_%s" % unknowns
                else:
                    found.setdefault(geneID, 0)
                    found[geneID] += 1
                    if found[geneID] > 1:
                        geneID += ".%d" % found[geneID]
                OUT.write('>' + geneID + '\n')
            else:
                OUT.write(line)
        

if __name__=='__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("gene_info", type=str, help=('string with wildcard,'
                        'for example ../gene_info/%%s_gene_info.tsv'))
    parser.add_argument("fastafiles", nargs="+")
    parser.add_argument("-o", "--outputformat", default="{0}_genes.fa",
                        help=("output file: '{0}' will be replaced by the "
                              "basename of the input file. [%(default)r]"))
    parser.add_argument("--shorten-species", action='store_true',
                        help="change 'Mus musculus' to 'mmusculus'?")
    parser.add_argument("-f", "--force-overwrite", action='store_true',
                        help="overwrite already existing files")
    parser.add_argument("--cprot", type=int, default=2, metavar='INT',
                        help="column for protein [%(default)s]")
    parser.add_argument("--cgene", type=int, default=0, metavar='INT',
                        help="column for gene [%(default)s]")

    args = parser.parse_args()
    gene_info = args.gene_info # global variable
    shorten_species = args.shorten_species # global variable
    force_overwrite = args.force_overwrite
    #for protID in sys.argv[2:]:
    for fastafile in args.fastafiles:
        print >>sys.stderr, fastafile
        rewrite_fastafile(fastafile, args.outputformat, args.cprot, args.cgene)


