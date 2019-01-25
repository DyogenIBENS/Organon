#!/usr/bin/env python3

"""Convert any Ensembl protein ID to gene ID and vice-versa
EXAMPLE:
    ./prot2gene ~/ws2/DUPLI_data85/gene_info/%s_gene_info.tsv <fastafiles>
"""

from __future__ import print_function

from sys import version_info, stderr, exit
import os.path as op
import argparse
from bz2 import BZ2File
from multiprocessing import Pool
import logging
logger = logging.getLogger(__name__)
logging.basicConfig(format='%(levelname)s:%(funcName)s:%(message)s')


#from genomicustools.identify import convert_prot2gene
from genomicustools.identify import convert_prot2species


ENSEMBL_VERSION = 85


class noop_output(object):
    """Object that behaves like an output file object but does absolutely nothing"""
    def __init__(self, *args, **kwargs):
        pass
    def __enter__(self):
        return self
    def __exit__(self, type, value, traceback):
        pass
    def write(self, *args):
        pass
    def writelines(self, *args):
        pass


def myopen(filename, *args, **kwargs):
    if filename.endswith('.bz2'):
        return BZ2File(filename, *args, **kwargs)
    else:
        return open(filename, *args, **kwargs)


def load_prot2gene(filename, cprot, cgene):
    conversions = {}
    with open(filename) as f:
        for line in f:
            fields = line.rstrip().split('\t')
            conversions[fields[cprot]] = fields[cgene]
    return conversions


def rewrite_fastafile(fastafile, gene_info, outputformat="{0}_genes.fa", cprot=2,
                      cgene=0, shorten_species=False,
                      ensembl_version=ENSEMBL_VERSION, force_overwrite=False,
                      verbose=1, strict=False, dryrun=False):
    if verbose:
        print(fastafile)
    genetree, ext = op.splitext(fastafile)
    if ext == '.bz2': genetree, ext = op.splitext(genetree)
    genetreedir, genetreefile = op.split(genetree)
    #print(genetree, genetreedir, genetreefile, file=stderr)
    outfile = outputformat.format(genetreefile)
    if op.exists(outfile):
        if force_overwrite:
            if not dryrun:
                logger.warning("(Overwriting %s)", outfile)
            else:
                logger.warning("(Dry-run: real run would overwrite %s)", outfile)
        else:
            logger.warning("%s exists. Skipping.", outfile)
            return

    # avoid duplicate genes
    found = {}
    unknowns = 0

    if version_info.major == 3 and fastafile.endswith('.bz2'):
        iter_lines = lambda F: (line.decode() for line in F)
    else:
        iter_lines = lambda F: F
    
    prot2gene = load_prot2gene(gene_info, cprot, cgene)

    with myopen(fastafile) as IN, \
            (noop_output() if dryrun else myopen(outfile, 'w')) as OUT:
        for line in iter_lines(IN):
            if line[0] == '>':
                protID = line[1:].split('/')[0]
                geneID = prot2gene.get(protID)
                         #convert_prot2gene(protID, gene_info, cprot, cgene,
                         #                  shorten_species, ensembl_version)
                #if not geneID and protID.startswith('ENSCSAP'):
                #    protID = protID.replace('ENSCSAP', 'ENSCSAVP')
                #    geneID = convert_prot2gene(protID)
                #    print("converting", geneID, file=stderr)
                #    if geneID:
                #        # Fit names in tree
                #        geneID = geneID.replace('ENSCSAVG', 'ENSCSAG')
                if not geneID:
                    notfound_msg = "Protein ID %s could not be converted" % protID
                    if strict:
                        raise LookupError(notfound_msg)
                    else:
                        logger.error(notfound_msg)
                    unknowns += 1
                    species = convert_prot2species(protID, ensembl_version, 'unknown')
                    geneID = "%s_gene_%d" % (species, unknowns)
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
    parser.add_argument('-n', '--dryrun', action='store_true')
    parser.add_argument("-o", "--outputformat", default="{0}_genes.fa",
                        help=("output file: '{0}' will be replaced by the "
                              "basename of the input file. [%(default)r]"))
    parser.add_argument("-f", "--force-overwrite", action='store_true',
                        help="overwrite already existing files")
    parser.add_argument("--cprot", type=int, default=2, metavar='INT',
                        help="column for protein [%(default)s]")
    parser.add_argument("--cgene", type=int, default=0, metavar='INT',
                        help="column for gene [%(default)s]")
    parser.add_argument("--strict", action='store_true',
                        help="Exit at first failed conversion")
    parser.add_argument('-e', '--ensembl-version', type=int,
                        default=ENSEMBL_VERSION, help='[%(default)s]')
    parser.add_argument("--shorten-species", action='store_true',
                        help="DEPRECATED. Change 'Mus musculus' to 'mmusculus'?")

    args = parser.parse_args()
    #for protID in argv[2:]:
    #for fastafile in args.fastafiles:
    #    print(fastafile, file=stderr)
    #    rewrite_fastafile(fastafile, args.outputformat, args.cprot, args.cgene)
    pool = Pool(processes=args.cores)
    if args.fromfile:
        if len(args.fastafiles) > 1:
            logger.error("Only one 'fastafiles' allowed with --fromfile. See help")
            exit(1)
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
                      args.strict,
                      args.dryrun) for f in fastafiles)
    pool.map(rewrite_fasta_process, generate_args)

