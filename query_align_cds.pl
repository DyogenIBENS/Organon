#! /usr/bin/perl

use strict;
use warnings;

use Bio::AlignIO;
use Bio::EnsEMBL::Registry;

use IO::Compress::Bzip2 qw(bzip2 $Bzip2Error) ;


# Set autoflush;
#$| = 1;


my $db_version = $ARGV[0] || 86 ;

$db_version += 0 ; # convert to integer?
print STDERR "Ensembl compara database version: $db_version\n" ;


my $reg = 'Bio::EnsEMBL::Registry';
$reg->load_registry_from_db(
    -host=>'ensembldb.ensembl.org',
    -user=>'anonymous', 
    -db_version=>$db_version,
);
my $tree_adaptor = $reg->get_adaptor('Multi', 'compara', 'GeneTree'); 


#open INFILE, "<all_gene_trees.txt" ;
#my $GeneTree=$ARGV[0];
#my $prev_outfile = '';

my $count = 0;

while (my $GeneTree = <STDIN>) {
	chomp $GeneTree ;
	my $outfile = "$GeneTree.fa" ;
	$count++ ;
	if ( (! -f $outfile) && (! -f "$outfile.bz2") ) {
		print STDERR "$count : Downloading $GeneTree.\n"; # Compressing '$prev_outfile'\n" ;		

		my $tree = $tree_adaptor->fetch_by_stable_id($GeneTree);
		$tree->preload;

		$tree = $tree->root;
		my $cds_align  = $tree->get_AlignedMemberSet()->get_SimpleAlign(-SEQ_TYPE => 'cds');

		#print "ALIGNMENT PROT\n";
		#print $alignIO $prot_align;

		#print "ALIGNMENT CDS\n";
		#open OUTFILE, ">$outfile" ;
		
		# Launch bzipping on previous file in parallel:
		#if ( $prev_outfile ) {
		#	system("bzip2 $prev_outfile &") == 0
		#		or die "bzip2 $prev_outfile failed $?\nReason: $!\n";
		#}

		my $alignIO = Bio::AlignIO->newFh(-file => ">$outfile", -format => "fasta");
		#my $z = new IO::Compress::Bzip2 $alignIO
		#	or die "IO::Compress::Bzip2 failed: $Bzip2Error\n";;
		#print $z $cds_align;
		#my $alignIO = Bio::AlignIO->newFh(-file => $z, -format => "fasta");
		print $alignIO $cds_align;
		bzip2 "$outfile" => "$outfile.bz2" or die "bzip2 failed: $Bzip2Error\n";
		system("rm $outfile") == 0 or die "rm failed ($?): $!\n";

		#$prev_outfile = $outfile ;
	} else {
		print STDERR "$count : $outfile already exists!\n" ;
		#$prev_outfile = '';
	}
}
#if ( $prev_outfile ) {
#	print STDERR "Compressing $prev_outfile\n" ;
#	system("bzip2 $prev_outfile &") == 0
#		or die "bzip2 $prev_outfile failed $?\nReason: $!\n";
#}

#close(INFILE)
