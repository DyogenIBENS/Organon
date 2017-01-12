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
from multiprocessing import Pool

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
    #print cprot, cgene
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


def convert_prot2gene(protID, gene_info, cprot=2, cgene=0, shorten_species=False):
    sp = convert_prot2species(protID)
    if shorten_species:
        spsplit = sp.split()
        sp2 = spsplit[0][0].lower() + spsplit[-1]
    else:
        sp2 = sp.replace(' ', '.')
    return grep_prot(gene_info % sp2, protID, cprot=cprot, cgene=cgene)


def rewrite_fastafile(fastafile, gene_info, outputformat="{0}_genes.fa", cprot=2,
                      cgene=0, shorten_species=False, force_overwrite=False,
                      verbose=1, strict=False):
    if verbose:
        print fastafile
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
                geneID = convert_prot2gene(protID, gene_info, cprot, cgene,
                                            shorten_species)
                #if not geneID and protID.startswith('ENSCSAP'):
                #    protID = protID.replace('ENSCSAP', 'ENSCSAVP')
                #    geneID = convert_prot2gene(protID)
                #    print >>sys.stderr, "converting", geneID
                #    if geneID:
                #        # Fit names in tree
                #        geneID = geneID.replace('ENSCSAVG', 'ENSCSAG')
                if not geneID:
                    if strict:
                        raise RuntimeError("protein ID %s could not be converted"\
                                            % protID)
                    unknowns += 1
                    geneID = "unknown_gene_%s" % unknowns
                else:
                    found.setdefault(geneID, 0)
                    found[geneID] += 1
                    if found[geneID] > 1:
                        geneID += ".%d" % found[geneID]
                if verbose > 1:
                    print "%s -> %s" % (protID, geneID)
                OUT.write('>' + geneID + '\n')
            else:
                OUT.write(line)

def rewrite_fasta_process(arglist):
    rewrite_fastafile(*arglist)


if __name__=='__main__':
    #fasta_rewriter = FastaRewriter()
    #fasta_rewriter.run()
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("gene_info", type=str, help=('string with wildcard,'
                        'for example ../gene_info/%%s_gene_info.tsv'))
    parser.add_argument("fastafiles", nargs="+")
    parser.add_argument("--fromfile", action='store_true',
                        help=("if True, the positional argument <fastafiles> "
                              "is a file containing one fastafile per line"))
    parser.add_argument("--cores", type=int, default=1, 
                        help="number of cores for parallelization")
    parser.add_argument("-q", "--quiet", action='store_const', const=0,
                        dest='verbose', default=1,
                        help="do not print each fasta file name")
    parser.add_argument("-v", "--verbose", action='store_const', const=2,
                        dest='verbose', default=1,
                        help="print each conversion")
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
    ##TODO: argument to trow error if conversion not found
    parser.add_argument("--strict", action='store_true',
                        help="Exit at first failed conversion")

    args = parser.parse_args()
    #for protID in sys.argv[2:]:
    #for fastafile in args.fastafiles:
    #    print >>sys.stderr, fastafile
    #    rewrite_fastafile(fastafile, args.outputformat, args.cprot, args.cgene)
    pool = Pool(processes=args.cores)
    if args.fromfile:
        if len(args.fastafiles) > 1:
            print >>sys.stderr, "Error: only one 'fastafiles' allowed with "\
                    "--fromfile. See help"
            sys.exit(1)
        else:
            with open(args.fastafiles[0]) as ff:
                fastafiles = [line.rstrip() for line in ff]
    else:
        fastafiles = args.fastafiles

#def _run_process(self, fastafile):
#    rewrite_fastafile(fastafile, **self.args)#fastafile, args.gene_infoargs.outputformat, args.cprot, args.cgene, verbose=True)

    generate_args = ((f,
                        args.gene_info,
                        args.outputformat,
                        args.cprot,
                        args.cgene,
                        args.shorten_species,
                        args.force_overwrite,
                        args.verbose,
                        args.strict) for f in fastafiles)
    pool.map(rewrite_fasta_process, generate_args)

