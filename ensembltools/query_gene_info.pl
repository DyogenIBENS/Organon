#!/usr/bin/env perl

# An example script demonstrating the use of BioMart API.
# This perl API representation is only available for configuration versions >=  0.5 
# Guillaume edit 2016/09/04: useful tuto at [https://www.biostars.org/p/53241/]

use warnings;
use strict;
use BioMart::Initializer;
use BioMart::Query;
use BioMart::QueryRunner;

my $confFile = $ENV{"HOME"} . "/.local/lib/perl/biomart-perl/conf/martURLLocation85.xml";
#
# NB: change action to 'clean' if you wish to start a fresh configuration
# and to 'cached' if you want to skip configuration step on subsequent runs from the same registry
#

my $action='cached';
my $initializer = BioMart::Initializer->new('registryFile'=>$confFile, 'action'=>$action);
my $registry = $initializer->getRegistry;


while (my $dataset = <>) {
	print STDERR $dataset ;
	chomp $dataset ;

	my $query = BioMart::Query->new('registry'=>$registry,'virtualSchemaName'=>'default');
		
	$query->setDataset("${dataset}_gene_ensembl");
	$query->addFilter("biotype", ["protein_coding"]);
	$query->addAttribute("ensembl_gene_id");
	$query->addAttribute("ensembl_transcript_id");
	$query->addAttribute("ensembl_peptide_id");
	$query->addAttribute("chromosome_name");
	$query->addAttribute("start_position");
	$query->addAttribute("end_position");
	#$query->addAttribute("strand");
	#$query->addAttribute("band");
	$query->addAttribute("external_gene_name");
	$query->addAttribute("external_gene_source");
	$query->addAttribute("transcript_count");
    #$query->addAttribute("percentage_gc_content");  #NOT FOUND
	$query->addAttribute("description");
	$query->addAttribute("go_id");
	$query->addAttribute("name_1006");
	$query->addAttribute("definition_1006");
	$query->addAttribute("go_linkage_type");
	$query->addAttribute("namespace_1003");
	$query->addAttribute("goslim_goa_accession");
	$query->addAttribute("hgnc_id");

	$query->formatter("TSV");

	my $query_runner = BioMart::QueryRunner->new();
	############################## GET COUNT ############################
	# $query->count(1);
	# $query_runner->execute($query);
	# print $query_runner->getCount();
	#####################################################################


	############################## GET RESULTS ##########################
	# to obtain unique rows only
	# $query_runner->uniqueRowsOnly(1);

	$query_runner->execute($query);
	open (OUT, ">${dataset}_gene_info.tsv") ;
	$query_runner->printHeader(\*OUT);
	$query_runner->printResults(\*OUT);
	$query_runner->printFooter(\*OUT);
	close OUT ;
	#####################################################################

    # bizarrement, une fois le fichier rempli, la connexion a la base de
    # donnees bloque, et l'iteration suivante n'est executee qu'apres
    # "Problems with the web server: 500 read timeout"
}

