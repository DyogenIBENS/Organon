#!/usr/bin/env python3


"""Module to help converting identifiers for proteins/genes/species from Genomicus/Ensembl.
"""

from __future__ import print_function

from sys import stderr
import os.path as op
import re
from copy import deepcopy
from bz2 import BZ2File
import warnings


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
            '': 'Caenorhabditis elegans',  # just so that the values reflect all species.
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

PROT2SP[90].update({'ENSMEUP': 'Notamacropus eugenii',  # Update of Macropus e.
                    'ENSCAPP': 'Cavia aperea',
                    'ENSCLAP': 'Chinchilla lanigera',
                    'ENSCGRP00001': 'Cricetulus griseus CHOK1GS',
                    'ENSCGRP00000': 'Cricetulus griseus Crigri',  # !!!!
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

GENE2SP = {}

GENE2SP[85] = {'Y': 'Saccharomyces cerevisiae',
               'Q': 'Saccharomyces cerevisiae',
               'FBgn': 'Drosophila melanogaster',
               'WBGene0': 'Caenorhabditis elegans',
               'ENSCING': 'Ciona intestinalis',
               'ENSCSAV': 'Ciona savignyi',
               'ENSCSAG': 'Chlorocebus sabaeus',
               'ENSPMAG': 'Petromyzon marinus',
               'ENSXETG': 'Xenopus tropicalis',
               'ENSPSIG': 'Pelodiscus sinensis',
               'ENSGALG': 'Gallus gallus',
               'ENSMGAG': 'Meleagris gallopavo',
               'ENSTGUG': 'Taeniopygia guttata',
               'ENSFALG': 'Ficedula albicollis',
               'ENSAPLG': 'Anas platyrhynchos',
               'ENSACAG': 'Anolis carolinensis',
               'ENSOANG': 'Ornithorhynchus anatinus',
               'ENSMEUG': 'Macropus eugenii',
               'ENSSHAG': 'Sarcophilus harrisii',
               'ENSMODG': 'Monodelphis domestica',
               'ENSLAFG': 'Loxodonta africana',
               'ENSETEG': 'Echinops telfairi',
               'ENSPCAG': 'Procavia capensis',
               'ENSDNOG': 'Dasypus novemcinctus',
               'ENSCHOG': 'Choloepus hoffmanni',
               'ENSSARG': 'Sorex araneus',
               'ENSEEUG': 'Erinaceus europaeus',
               'ENSMLUG': 'Myotis lucifugus',
               'ENSPVAG': 'Pteropus vampyrus',
               'ENSTTRG': 'Tursiops truncatus',
               'ENSBTAG': 'Bos taurus',
               'ENSOARG': 'Ovis aries',
               'ENSVPAG': 'Vicugna pacos',
               'ENSSSCG': 'Sus scrofa',
               'ENSECAG': 'Equus caballus',
               'ENSMPUG': 'Mustela putorius furo',
               'ENSAMEG': 'Ailuropoda melanoleuca',
               'ENSCAFG': 'Canis lupus familiaris',
               'ENSFCAG': 'Felis catus',
               'ENSTBEG': 'Tupaia belangeri',
               'ENSPANG': 'Papio anubis',
               'ENSMMUG': 'Macaca mulatta',
               'ENSPPYG': 'Pongo abelii',
               'ENSGGOG': 'Gorilla gorilla gorilla',
               'ENSPTRG': 'Pan troglodytes',
               'ENSG000':    'Homo sapiens',       # ENSG
               'ENSNLEG': 'Nomascus leucogenys',
               'ENSCJAG': 'Callithrix jacchus',
               'ENSTSYG': 'Tarsius syrichta',
               'ENSOGAG': 'Otolemur garnettii',
               'ENSMICG': 'Microcebus murinus',
               'ENSOPRG': 'Ochotona princeps',
               'ENSOCUG': 'Oryctolagus cuniculus',
               'ENSCPOG': 'Cavia porcellus',
               'ENSRNOG': 'Rattus norvegicus',
               'ENSMUSG': 'Mus musculus',
               'ENSSTOG': 'Ictidomys tridecemlineatus',
               'ENSDORG': 'Dipodomys ordii',
               'ENSLACG': 'Latimeria chalumnae',
               'ENSLOCG': 'Lepisosteus oculatus',
               'ENSGACG': 'Gasterosteus aculeatus',
               'ENSTNIG': 'Tetraodon nigroviridis',
               'ENSTRUG': 'Takifugu rubripes',
               'ENSONIG': 'Oreochromis niloticus',
               'ENSORLG': 'Oryzias latipes',
               'ENSPFOG': 'Poecilia formosa',
               'ENSXMAG': 'Xiphophorus maculatus',
               'ENSGMOG': 'Gadus morhua',
               'ENSAMXG': 'Astyanax mexicanus',
               'ENSDARG': 'Danio rerio'}

GENE2SP[86] = deepcopy(GENE2SP[85])
GENE2SP[86].update(**{'MGP_SPR': 'Mus spretus'})

GENE2SP[87] = deepcopy(GENE2SP[86])
GENE2SP[87].update({'ENSTSYG': 'Carlito syrichta'})

GENE2SP[90] = GENE2SP[88] = GENE2SP[87]

GENE2SP[90].update({'ENSMEUG': 'Notamacropus eugenii',  # Update of Macropus e.
                    'ENSCAPG': 'Cavia aperea',
                    'ENSCLAG': 'Chinchilla lanigera',
                    'ENSCGRG00001': 'Cricetulus griseus CHOK1GS',
                    'ENSCGRG00000': 'Cricetulus griseus Crigri',  # !!!!
                    'ENSFDAG':      'Fukomys damarensis',
                    'ENSHGLG00000': 'Heterocephalus glaber female',
                    'ENSHGLG00100': 'Heterocephalus glaber male',
                    'ENSJJAG': 'Jaculus jaculus',
                    'ENSMAUG': 'Mesocricetus auratus',
                    'ENSMOCG': 'Microtus ochrogaster',
                    'MGP_CAR': 'Mus caroli',
                    'MGP_Pah': 'Mus pahari',
                    'ENSNGAG': 'Nannospalax galili',
                    'ENSODEG': 'Octodon degus',
                    'ENSPEMG': 'Peromyscus maniculatus bairdii'})

SP2GENEID = {}
for version, conversion in GENE2SP.items():
    SP2GENEID[version] = {}
    for geneid, species in conversion.items():
        if species in SP2GENEID[version]:
            SP2GENEID[version][species] += '|' + geneid
        else:
            SP2GENEID[version][species] = geneid
#from tabletools import inverse


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
                raise KeyError(err.args[0] + ' (protein: %s, Ens.%d)' % \
                                (modernID, ensembl_version))


def convert_gene2species(modernID, ensembl_version=ENSEMBL_VERSION):
    """Give the species of an **extant** sequence id (from ENSEMBL)."""
    gene2sp = GENE2SP[ensembl_version]
    if ensembl_version >= 90:
        try:
            return gene2sp[modernID[:12]]
        except KeyError:
            pass
    try:
        return gene2sp[modernID[:7]]
    except KeyError:
        try:
            return gene2sp[modernID[:4]]
        except KeyError:
            try:
                return gene2sp[modernID[0]]
            except KeyError:
                raise KeyError("%s can't be assigned to a species" % modernID)


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


def load_assembly2species(filename="~/ws2/UCSC_genome_releases_full.tsv",
                          fromcol=2, tocol=0):
    conversion = {}
    with open(op.expanduser(filename)) as stream:
        header = next(stream)
        for line in stream:
            fields = line.rstrip().split('\t')
            conversion[fields[fromcol]] = fields[tocol]
    return conversion


ucsc_conv_filename = '~/ws2/UCSC_genome_releases_full.tsv'
try:
    assembly2species = load_assembly2species(ucsc_conv_filename)
except FileNotFoundError:
    warnings.warn("Conversion file not found: %r" % ucsc_conv_filename)
    assembly2species = {}


def ultimate_seq2sp(seqname, ensembl_version=ENSEMBL_VERSION):
    """From an extant sequence name, find the corresponding species.
    Recognizes Ensembl gene IDs, Ensembl protein IDs, and also UCSC assembly
    names such as 'loxAfr3'"""
    try:
        sp = convert_gene2species(seqname, ensembl_version)
    except KeyError:
        try:
            sp = convert_prot2species(seqname, ensembl_version)
        except KeyError:
            assembly = re.match('[A-Za-z0-9]+', seqname).group()
            sp = assembly2species[assembly]
    return sp


# My first ever unit-test!
def test_convert_prot2species(ensembl_version, default, gene_info, cprot=2):
    """test the `convert_prot2species` function for every modernID"""
    # Check rejection of wrong strings
    from glob import glob
    for wrong in ('xululul', '0000000', 'ENSXXXP', 'ENSG000'):
        predicted_sp = convert_prot2species(wrong, ensembl_version, False)
        assert predicted_sp is False, "%r predicted %r" % (wrong, predicted_sp)

    # Check every valid prot ID in the given files
    splist_module = set(PROT2SP[ensembl_version].values())
    sp_from_filename = re.compile(gene_info.replace('%s', '([A-Za-z0-9.]+)'))
    gene_info_files = glob(gene_info.replace('%s', '*'))
    assert gene_info_files, "No files found, check your path."
    splist_file = set(sp_from_filename.match(fn).group(1).replace('.', ' ')
                      for fn in gene_info_files)

    if not splist_module == splist_file:
        raise AssertionError('Differences in lists of species:\n' +
                             'module (%d): %s\n' % (len(splist_module),
                                                    splist_module - splist_file) +
                             'files  (%d): %s' % (len(splist_file),
                                                  splist_file - splist_module))
    for sp in splist_file:
        filename = gene_info % sp.replace(' ', '.')
        print("Checking %s in %r" % (sp, op.basename(filename)), file=stderr)
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


def test_convert_gene2species(ensembl_version, gene_info, cgene=1):
    """test the above function for every modernID"""
    # Check rejection of wrong strings
    from glob import glob
    for wrong in ('xululul', '0000000', 'ENSXXXG', 'ENSP000'):
        try:
            predicted_sp = convert_gene2species(wrong, ensembl_version)
        except KeyError:
            predicted_sp = False

        assert predicted_sp is False, "%r predicted %r" % (wrong, predicted_sp)

    # Check every valid gene ID in the given files
    splist_module = set(GENE2SP[ensembl_version].values())
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
        print("Checking %s in %r" % (sp, op.basename(filename)), file=stderr)
        # Check that each species gene return the correct species.
        with open(filename) as IN:
            for line in IN:
                gene = line.rstrip('\r\n').split('\t')[cgene]
                try:
                    predicted_sp = convert_gene2species(gene, ensembl_version)
                except KeyError as err:
                    err.args = err.args[:-1] + \
                               (err.args[-1] + ' '.join((sp, gene, "Not found")),)
                    raise
                assert sp == predicted_sp, "%s: %r ≠ %r" % (gene, sp, predicted_sp)
    return True
