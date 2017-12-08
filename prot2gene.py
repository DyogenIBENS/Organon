#!/usr/bin/env python3

"""Convert any Ensembl protein ID to gene ID and vice-versa
EXAMPLE:
    ./prot2gene ~/ws2/DUPLI_data85/gene_info/%s_gene_info.tsv <fastafiles>
"""

from __future__ import print_function

import re
import sys
import os.path
import argparse
from copy import deepcopy
from bz2 import BZ2File
from multiprocessing import Pool
from glob import glob

def myopen(filename, *args, **kwargs):
    if filename.endswith('.bz2'):
        return BZ2File(filename, *args, **kwargs)
    else:
        return open(filename, *args, **kwargs)

ENSEMBL_VERSION = 85

CELEGANS_REG = re.compile(r'^([234]R|A[CH]|B[0E]|C[0-5C-E]|cT|[DFHKMRTWY][0-9CHYT]|E[0_EG]|JC|LL|SS|Z[CK]|V[BCFHMWYZ]|P[AD]).*')

PROT2SP = {85:
           { #'Y': 'Saccharomyces cerevisiae',  # there are C.elegans prot with 'Y' too
            'Q0': 'Saccharomyces cerevisiae',
            'FB': 'Drosophila melanogaster',
            #'WBGene0': 'Caenorhabditis elegans',  # No consensus
            '': 'Caenorhabditis elegans', # just so that the values reflect all species.
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
           }

PROT2SP[86] = deepcopy(PROT2SP[85])
PROT2SP[86].update(**{'MGP_SPR': 'Mus spretus'})  # TODO: check the protein id

PROT2SP[87] = deepcopy(PROT2SP[86])
PROT2SP[87].update(**{'ENSTSYP': 'Carlito syrichta'})

PROT2SP[90] = PROT2SP[88] = PROT2SP[87]

PROT2SP[90].update({'ENSMEUP': 'Notamacropus eugenii', # Update of Macropus e.
                    'ENSCAPP': 'Cavia aperea',
                    'ENSCLAP': 'Chinchilla lanigera',
                    'ENSCGRP00001': 'Cricetulus griseus CHOK1GS',
                    'ENSCGRP00000': 'Cricetulus griseus Crigri', # !!!!
                    'ENSFDAP':      'Fukomys damarensis',
                    'ENSHGLP00000': 'Heterocephalus glaber female',
                    'ENSHGLP00100': 'Heterocephalus glaber male',
                    'ENSJJAP': 'Jaculus jaculus',
                    'ENSMAUP': 'Mesocricetus auratus',
                    'ENSMOCP': 'Microtus ochrogaster',
                    'MGP_CAR': 'Mus caroli',
                    'MGP_Pah': 'Mus pahari',
                    'ENSNGAP': 'Nannospalax galili',
                    'ENSODEP': 'Octodon degus',
                    'ENSPEMP': 'Peromyscus maniculatus bairdii'})

def convert_prot2species(modernID, ensembl_version=ENSEMBL_VERSION, default=None):
    if ensembl_version >= 90:
        try:
            return PROT2SP[ensembl_version][modernID[:12]]
        except KeyError:
            pass
    try:
        return PROT2SP[ensembl_version][modernID[:7]]
    except KeyError:
        try:
            # Saccharomyces cerevisiae (Q0) or Drosophila melanogaster
            return PROT2SP[ensembl_version][modernID[:2]]
        except KeyError as err:
            if re.match('Y[A-Z]', modernID):
                return 'Saccharomyces cerevisiae'
            elif CELEGANS_REG.match(modernID):
                return 'Caenorhabditis elegans'
            elif default is not None:
                return default
            else:
                err.args = (err.args[0] + \
                            ' (protein: %s, Ens.%d)' % (modernID, ensembl_version),)
                raise


# My first ever unit-test!
def test_convert_prot2species(ensembl_version, default, gene_info, cprot=2):
    """test the above function for every modernID"""
    # Check rejection of wrong strings
    for wrong in ('xululul', '0000000', 'ENSXXXP', 'ENSG000'):
        predicted_sp = convert_prot2species(wrong, ensembl_version, False)
        assert predicted_sp is False, "%r predicted %r" % (wrong, predicted_sp)

    # Check every valid prot ID in the given files
    splist_module = set(PROT2SP[ensembl_version].values())
    sp_from_filename = re.compile(gene_info.replace('%s', '([A-Za-z0-9.]+)'))
    gene_info_files = glob(gene_info.replace('%s', '*'))
    assert gene_info_files, "No files found, check your path."
    splist_file = set(sp_from_filename.match(fn).group(1).replace('.', ' ') \
                        for fn in gene_info_files)

    if not splist_module == splist_file:
        raise AssertionError('Differences in lists of species:\n' +
                             'module (%d): %s\n' % (len(splist_module),
                                                    splist_module - splist_file) +
                             'files  (%d): %s' % (len(splist_file),
                                                  splist_file - splist_module))
    for sp in splist_file:
        filename = gene_info % sp.replace(' ', '.')
        print("Checking %s in %r" % (sp, os.path.basename(filename)), file=sys.stderr)
        # Check that each species protein return the correct species.
        with myopen(filename) as IN:
            for line in IN:
                prot = line.rstrip('\r\n').split('\t')[cprot]
                try:
                    predicted_sp = convert_prot2species(prot, ensembl_version, default)
                    assert sp == predicted_sp, "%s: %r ≠ %r" % (prot, sp, predicted_sp)
                except KeyError as err:
                    err.args = err.args[:-1] + \
                               (err.args[-1] + ' '.join((sp, prot, "Not found")),)
                    raise
    return True


def grep_prot(filename, protID, cprot=2, cgene=0):
    #print(cprot, cgene)
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


def convert_prot2gene(protID, gene_info, cprot=2, cgene=0, shorten_species=False,
                      ensembl_version=ENSEMBL_VERSION):
    sp = convert_prot2species(protID, ensembl_version)
    if shorten_species:
        spsplit = sp.split()
        sp2 = spsplit[0][0].lower() + spsplit[-1]
    else:
        sp2 = sp.replace(' ', '.')
    return grep_prot(gene_info % sp2, protID, cprot=cprot, cgene=cgene)


def rewrite_fastafile(fastafile, gene_info, outputformat="{0}_genes.fa", cprot=2,
                      cgene=0, shorten_species=False,
                      ensembl_version=ENSEMBL_VERSION, force_overwrite=False,
                      verbose=1, strict=False):
    if verbose:
        print(fastafile)
    genetree, ext = os.path.splitext(fastafile)
    if ext == '.bz2': genetree, ext = os.path.splitext(genetree)
    genetreedir, genetreefile = os.path.split(genetree)
    #print >>sys.stderr, genetree, genetreedir, genetreefile
    outfile = outputformat.format(genetreefile)
    if os.path.exists(outfile):
        if force_overwrite:
            print("(Overwriting %s)" % outfile, file=sys.stderr)
        else:
            print("%s exists. Skipping." % outfile, file=sys.stderr)
            return

    # avoid duplicate genes
    found = {}
    unknowns = 0
    with myopen(fastafile) as IN, myopen(outfile, 'w') as OUT:
        for line in IN:
            if line[0] == '>':
                protID = line[1:].split('/')[0]
                geneID = convert_prot2gene(protID, gene_info, cprot, cgene,
                                           shorten_species, ensembl_version)
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
                    print("%s -> %s" % (protID, geneID))
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
    parser.add_argument('-e', '--ensembl-version', type=int,
                        default=ENSEMBL_VERSION, help='[%(default)s]')
    
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
            print("Error: only one 'fastafiles' allowed with --fromfile. See help", file=sys.stderr)
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
                        args.ensembl_version,
                        args.force_overwrite,
                        args.verbose,
                        args.strict) for f in fastafiles)
    pool.map(rewrite_fasta_process, generate_args)

