#! /usr/bin/perl

use strict;
use warnings;

use Bio::AlignIO;
use Bio::EnsEMBL::Registry;

use IO::Compress::Bzip2 qw(bzip2 $Bzip2Error) ;

# Set autoflush;
$| = 1;


my $db_version = $ARGV[0] || 86 ;

$db_version += 0 ; # convert to integer?
print STDERR "Ensembl compara database version: $db_version\n" ;


my $reg = 'Bio::EnsEMBL::Registry';
$reg->load_registry_from_db(
    -host=>'ensembldb.ensembl.org',
    -user=>'anonymous', 
    -db_version=>$db_version,
);

sub clean_exit {
	my ($sig) = @_;
	print STDERR "Caught a SIG$sig. Closing connections and exiting.";
	$reg->disconnect_all(); # NOT SURE HERE. CHECK.
	exit(1);
}

sub clean_die {
	my ($msg) = @_;
	print STDERR "Die for reason: $msg.\nClosing connections and exiting.";
	$reg->disconnect_all(); # NOT SURE HERE. CHECK.
	exit(1);
}

#$SIG{'INT'}   = \&clean_exit; # Those signals are caught by sigtrap normal-signals already
#$SIG{'TERM'}  = \&clean_exit;
$SIG{__DIE__} = \&clean_die;
use sigtrap qw(handler clean_exit normal-signals error-signals);
#use sigtrap qw(handler clean_exit error-signals);

# Arg [1]     : name of the species to add the adaptor to in the registry.
# Arg [2]     : name of the group   to add the adaptor to in the registry.
# Arg [3]     : name of the type    to add the adaptor to in the registry.
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
		# fetch_all?
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

		my $alignIO = Bio::AlignIO->new(-file => ">$outfile", -format => "fasta");
		#my $z = new IO::Compress::Bzip2 $alignIO
		#	or die "IO::Compress::Bzip2 failed: $Bzip2Error\n";;
		#print $z $cds_align;
		#my $alignIO = Bio::AlignIO->newFh(-file => $z, -format => "fasta");
		$alignIO->write_aln($cds_align);
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
$reg->disconnect_all(); # NOT SURE HERE. CHECK.
