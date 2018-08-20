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

from genomicustools.identify import convert_prot2gene


ENSEMBL_VERSION = 85


def myopen(filename, *args, **kwargs):
    if filename.endswith('.bz2'):
        return BZ2File(filename, *args, **kwargs)
    else:
        return open(filename, *args, **kwargs)


def rewrite_fastafile(fastafile, gene_info, outputformat="{0}_genes.fa", cprot=2,
                      cgene=0, shorten_species=False,
                      ensembl_version=ENSEMBL_VERSION, force_overwrite=False,
                      verbose=1, strict=False):
    if verbose:
        print(fastafile)
    genetree, ext = op.splitext(fastafile)
    if ext == '.bz2': genetree, ext = op.splitext(genetree)
    genetreedir, genetreefile = op.split(genetree)
    #print(genetree, genetreedir, genetreefile, file=stderr)
    outfile = outputformat.format(genetreefile)
    if op.exists(outfile):
        if force_overwrite:
            print("(Overwriting %s)" % outfile, file=stderr)
        else:
            print("%s exists. Skipping." % outfile, file=stderr)
            return

    # avoid duplicate genes
    found = {}
    unknowns = 0

    if version_info.major == 3 and fastafile.endswith('.bz2'):
        iter_lines = lambda F: (line.decode() for line in F)
    else:
        iter_lines = lambda F: F

    with myopen(fastafile) as IN, myopen(outfile, 'w') as OUT:
        for line in iter_lines(IN):
            if line[0] == '>':
                protID = line[1:].split('/')[0]
                geneID = convert_prot2gene(protID, gene_info, cprot, cgene,
                                           shorten_species, ensembl_version)
                #if not geneID and protID.startswith('ENSCSAP'):
                #    protID = protID.replace('ENSCSAP', 'ENSCSAVP')
                #    geneID = convert_prot2gene(protID)
                #    print("converting", geneID, file=stderr)
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
    #for protID in argv[2:]:
    #for fastafile in args.fastafiles:
    #    print(fastafile, file=stderr)
    #    rewrite_fastafile(fastafile, args.outputformat, args.cprot, args.cgene)
    pool = Pool(processes=args.cores)
    if args.fromfile:
        if len(args.fastafiles) > 1:
            print("Error: only one 'fastafiles' allowed with --fromfile. See help", file=stderr)
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
                      args.strict) for f in fastafiles)
    pool.map(rewrite_fasta_process, generate_args)

