#! /usr/bin/perl

use strict;
use warnings;

use Bio::AlignIO;

use Bio::EnsEMBL::Registry;

my $reg = 'Bio::EnsEMBL::Registry';
$reg->load_registry_from_db(
    -host=>'ensembldb.ensembl.org',
    -user=>'anonymous', 
    -db_version=>85,
);
my $tree_adaptor = $reg->get_adaptor('Multi', 'compara', 'GeneTree'); 


#open INFILE, "<all_gene_trees.txt" ;
#my $GeneTree=$ARGV[0];
foreach (my $GeneTree = <STDIN>) {
	chomp $GeneTree ;
	my $outfile = "$GeneTree.fa" ;
	if ( (! -f $outfile) && (! -f "$outfile.bz2") ) {
		print STDERR "$GeneTree\n" ;
		
		my $tree = $tree_adaptor->fetch_by_stable_id($GeneTree);
		$tree->preload;

		$tree = $tree->root;
		#my $prot_align = $tree->get_AlignedMemberSet()->get_SimpleAlign();
		my $cds_align  = $tree->get_AlignedMemberSet()->get_SimpleAlign(-SEQ_TYPE => 'cds');

		#print "ALIGNMENT PROT\n";
		#print $alignIO $prot_align;

		#print "ALIGNMENT CDS\n";
		#open OUTFILE, ">$outfile" ;
		my $alignIO = Bio::AlignIO->newFh(-file => ">$outfile", -format => "fasta");
		print $alignIO $cds_align;
	} else {
		print STDERR "$outfile already exists!\n" ;
	}
}
#close(INFILE)
