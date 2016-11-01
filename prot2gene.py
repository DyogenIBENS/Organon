#!/usr/bin/env python

"""Convert any Ensembl protein ID to gene ID and vice-versa
EXAMPLE:
    ./prot2gene ~/ws2/DUPLI_data85/gene_info/%s_gene_info.tsv ENSMODP00000031002
"""

import sys
import os.path
import argparse


def convert_prot2species(modernID):
    prot2sp = { 'Y': 'Saccharomyces cerevisiae',
                'Q': 'Saccharomyces cerevisiae',
                'FBpp': 'Drosophila melanogaster',
                #'WBGene0': 'Caenorhabditis elegans',  # No consensus
                'ENSCINP': 'Ciona intestinalis',
                'ENSCSAP': 'Ciona savignyi',
                'ENSPMAP': 'Petromyzon marinus',
                'ENSXETP': 'Xenopus tropicalis',
                'ENSPSIP': 'Pelodiscus sinensis',
                'ENSGALP': 'Gallus gallus',
                'ENSMGAP': 'Meleagris gallopavo',
                'ENSTGUP': 'Taeniopygia guttata',
                'ENSFALP': 'Ficedulla albicollis',
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
            return prot2sp[modernID[:4]]
        except KeyError:
            try:
                return prot2sp[modernID[0]]
            except KeyError:
                return 'Caenorhabditis elegans'


def grep_prot(filename, protID):
    with open(filename) as IN:
        for line in IN:
            fields = line.split('\t')
            if fields[2] == protID:
                return fields[0]


def grep_gene(filename, geneID):
    with open(filename) as IN:
        for line in IN:
            fields = line.split('\t')
            if fields[0] == geneID:
                return fields[2]


def convert_prot2gene(protID):
    sp = convert_prot2species(protID)
    spsplit = sp.split()
    sp2 = spsplit[0][0].lower() + spsplit[-1]
    return grep_prot(gene_info % sp2, protID)


def rewrite_fastafile(fastafile, outputformat="{0}_genes.fa"):
    # avoid duplicate genes
    found = {}
    genetree, ext = os.path.splitext(fastafile)
    genetreedir, genetreefile = os.path.split(genetree)
    #print >>sys.stderr, genetree, genetreedir, genetreefile
    outfile = outputformat.format(genetreefile)
    unknowns = 0
    with open(fastafile) as IN, open(outfile, 'w') as OUT:
        for line in IN:
            if line[0] == '>':
                protID = line[1:].split('/')[0]
                geneID = convert_prot2gene(protID)
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
    parser.add_argument("gene_info", type=str, help=("string with wildcard,"
                        "for example ../gene_info/%s_gene_info.tsv"))
    parser.add_argument("fastafiles", nargs="+")
    parser.add_argument("-o", "--outputformat", default="{0}_genes.fa",
                        help=("output file: '{0}' will be replaced by the "
                              "basename of the input file"))
    args = parser.parse_args()
    gene_info = args.gene_info # global variable
    #for protID in sys.argv[2:]:
    for fastafile in args.fastafiles:
        print >>sys.stderr, fastafile
        rewrite_fastafile(fastafile, args.outputformat)


